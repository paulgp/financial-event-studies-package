# Many-event batch mode ---------------------------------------------------

#' Estimate many events in one parallelized call
#'
#' Fits [event_study()] once per event (each with its own event date, treated
#' unit(s), and donor pool), then averages ATT paths across events with
#' cross-event standard errors — the design of the index-inclusion and M&A
#' applications. Per-event inference is deliberately off (cross-event
#' variation is the SE; placebo reps per event do not scale to thousands of
#' events).
#'
#' @param data,unit,time,ret As in [event_study()]: one long panel containing
#'   all events' windows.
#' @param events Data frame: `unit`, `event_time`, optional `event` id (rows
#'   sharing an id form one multi-unit event). Extra columns (e.g. decade,
#'   deal type) are carried through to `$events` for grouped summaries.
#' @param method,window,est_window,returns,cumulate,align,factors,match_on,V,solver,lambda,r,force
#'   As in [event_study()], applied to every event.
#' @param exclude_treated Drop from each event's donor pool any unit with its
#'   own event between `est_window[1]` and `window[2]` of this event (in
#'   trading-period offsets), so contemporaneously-treated units never serve
#'   as donors.
#' @param se `"cross"` (default): cross-event mean and SE per horizon
#'   (`sd/sqrt(n)`, or its weighted analog under `weights`); `"none"`. With
#'   `event_se` active, `$se` also carries `within_att` (the pooled
#'   within-event SE component, `sqrt(sum(wn^2 se_i^2))`) and `total_att`
#'   (`pmax(att, within_att)`). The cross SE already embeds each event's
#'   estimation noise when events are independent — the two components are
#'   deliberately not added — but the cross dispersion is itself estimated,
#'   so `total_att` floors it by the known within component
#'   (DerSimonian-Laird method-of-moments logic), which binds when few
#'   events are pooled.
#' @param event_se Per-event inference passed to [event_study()]'s `se`
#'   (default `"none"`, the fast path). When active, each event's SE path
#'   lands in `$ses`, its `att_avg_se` in `$events`, and the pooled
#'   `within_att`/`total_att` fields are filled. `"conformal"` is not
#'   supported here (it stores CIs, not SEs).
#' @param event_reps Per-event `reps` passed to [event_study()] when
#'   `event_se` needs draws (placebo/bootstrap).
#' @param weights Per-event weights for the pooled paths and SEs:
#'   `NULL` (equal, the default), `"n_treated"` (each event's realized
#'   treated count), or the name of a numeric column of `events`. Applied
#'   consistently to `$att`/`$car`, the cross/within SEs, and [summary()].
#' @param cores Events are fit with `parallel::mclapply(mc.cores = cores)`,
#'   in blocks, with each worker pinned to one data.table thread (a full
#'   thread pool per fork oversubscribes the CPU and thrashes memory).
#' @param gsynth_cores Worker cap when `method = "gsynth"`, whose fits carry
#'   several times the memory of the other engines; default
#'   `max(1, cores %/% 3)`. Set to `cores` to disable the derate.
#' @param checkpoint_dir Optional directory for per-event checkpoints
#'   (`event_<i>.rds`). Completed events are reloaded instead of refit on a
#'   rerun, so a crash (or an OOM-killed worker) costs one block of work
#'   instead of the whole run.
#' @param keep_fits Store the full per-event `fes_fit` list (memory-heavy).
#' @return An object of class `fes_batch`: `$att`/`$car` (pooled paths),
#'   `$se`, `$atts`/`$cars` (event x horizon matrices), `$ses` (per-event SE
#'   paths when `event_se` is active), `$events` (input table + n_treated,
#'   n_donors, pre_rmse, att_avg_se, status), `$weights` (per-event weights
#'   when set), `$fits` if kept.
#' @export
event_study_batch <- function(data, unit, time, ret, events,
                              method = c("mean", "did", "market", "factor",
                                         "sc", "ridge", "sdid", "gsynth", "cfm",
                                         "apm"),
                              window = c(0, 10), est_window = c(-250, -11),
                              returns, cumulate = c("auto", "sum", "compound", "log"),
                              align = c("position", "value"),
                              factors = NULL, exclude_treated = TRUE,
                              match_on = c("ret", "cumret"), V = NULL,
                              solver = c("hybrid", "fw", "qp"), lambda = NULL,
                              r = c(0, 5), force = c("unit", "none", "two-way"),
                              se = c("cross", "none"),
                              event_se = "none", event_reps = NULL,
                              weights = NULL, cores = 1L,
                              gsynth_cores = NULL, checkpoint_dir = NULL,
                              keep_fits = FALSE) {
  method <- match.arg(method)
  se <- match.arg(se)
  align <- match.arg(align)
  event_se <- match.arg(event_se, c("none", "auto", "tstat", "placebo",
                                    "bootstrap", "analytic"))
  events <- as.data.frame(events)
  stopifnot(all(c("unit", "event_time") %in% names(events)))
  # [["event"]], not $event: partial matching would hit `event_time` whenever
  # the id column is absent, so the id would never be created and the final
  # merge() on "event" would fail for the documented minimal events table
  if (is.null(events[["event"]]))
    events$event <- seq_len(nrow(events))
  ids <- unique(events$event)

  # offsets between any two event times, for donor exclusion
  times <- sort(unique(data[[time]]))
  ev_pos <- if (align == "value") as.numeric(events$event_time)
            else match(events$event_time, times)
  if (anyNA(ev_pos)) stop("some `event_time`s not found among panel times")

  # the unit universe for donor exclusion is invariant across events: scan
  # the panel once, not once per event (a full-column scan per event is ~27s
  # on the 66.7M-row index-inclusion panel)
  all_units <- if (exclude_treated) unique(as.character(data[[unit]]))

  fit_one <- function(eid) {
    # one data.table thread per worker: `cores` forks x the default thread
    # pool oversubscribes the CPU and thrashes memory on wide panels;
    # restored on exit for the serial path
    old_dt <- data.table::setDTthreads(1L)
    on.exit(data.table::setDTthreads(old_dt), add = TRUE)
    rows <- events$event == eid
    s <- ev_pos[rows][1]
    donors <- NULL
    if (exclude_treated) {
      contam <- ev_pos >= s + est_window[1] & ev_pos <= s + window[2]
      bad <- unique(as.character(events$unit[contam]))
      donors <- setdiff(all_units, bad)
    }
    tryCatch({
      f <- event_study(data, unit, time, ret,
                       treated = as.character(events$unit[rows]),
                       event_time = events$event_time[rows][1],
                       method = method, window = window,
                       est_window = est_window, returns = returns,
                       cumulate = cumulate, align = align, factors = factors,
                       donors = donors, match_on = match_on, V = V,
                       solver = solver, lambda = lambda, r = r, force = force,
                       se = event_se, reps = event_reps, keep_data = FALSE)
      list(event = eid, att = f$att, car = f$car,
           se_att = if (!is.null(f$se)) f$se$att,
           att_avg_se = f$att_avg_se %||% NA_real_,
           n_treated = f$diagnostics$n_treated,
           n_donors = f$diagnostics$n_donors,
           pre_rmse = f$diagnostics$info$pre_rmse %||% NA_real_,
           status = "ok", fit = if (keep_fits) f)
    }, error = function(e)
      list(event = eid, att = NULL, car = NULL, se_att = NULL,
           att_avg_se = NA_real_, n_treated = NA_integer_,
           n_donors = NA_integer_, pre_rmse = NA_real_,
           status = paste0("dropped: ", conditionMessage(e)), fit = NULL))
  }

  eff_cores <- as.integer(cores)
  if (method == "gsynth")
    # a gsynth fit carries several times the memory of the other engines
    # (internal factor-number CV over refits); a full complement of gsynth
    # workers is what OOM-killed the 32 GB index-inclusion batch, so derate
    # unless the caller sets `gsynth_cores` explicitly
    eff_cores <- if (is.null(gsynth_cores)) max(1L, eff_cores %/% 3L)
                 else as.integer(gsynth_cores)

  if (!is.null(checkpoint_dir))
    dir.create(checkpoint_dir, showWarnings = FALSE, recursive = TRUE)
  fit_or_load <- function(i) {
    ck <- if (!is.null(checkpoint_dir))
      file.path(checkpoint_dir, sprintf("event_%d.rds", i))
    if (!is.null(ck) && file.exists(ck)) {
      x <- tryCatch(readRDS(ck), error = function(e) NULL)
      # files are keyed by position; the stored id guards against a reordered
      # events table silently serving another event's cached fit
      if (!is.null(x) && identical(x$event, ids[i])) return(x)
    }
    x <- fit_one(ids[i])
    if (!is.null(ck)) saveRDS(x, ck)
    x
  }
  # evaluate in blocks so a killed worker (OOM) costs at most one block, not
  # the whole run, and checkpoints (if any) land as each block completes
  block <- max(eff_cores, 1L) * 4L
  res <- vector("list", length(ids))
  for (b in split(seq_along(ids), ceiling(seq_along(ids) / block))) {
    res[b] <- if (eff_cores > 1L)
      parallel::mclapply(b, fit_or_load, mc.cores = eff_cores)
    else lapply(b, fit_or_load)
  }
  killed <- vapply(res, is.null, TRUE)
  if (any(killed)) {
    # a fork killed by the OS returns NULL with no condition at all — name
    # the likely cause loudly instead of failing downstream on shape errors
    warning(sum(killed), " event(s) lost to killed worker(s), likely out of ",
            "memory; recorded as dropped. Reduce `cores` (or `gsynth_cores`)",
            if (is.null(checkpoint_dir))
              " and consider `checkpoint_dir` to make progress durable"
            else "; completed events are checkpointed and will be reused")
    res[killed] <- lapply(which(killed), function(i)
      list(event = ids[i], att = NULL, car = NULL, se_att = NULL,
           att_avg_se = NA_real_, n_treated = NA_integer_,
           n_donors = NA_integer_, pre_rmse = NA_real_,
           status = "dropped: worker killed (likely out of memory)",
           fit = NULL))
  }

  ok <- vapply(res, function(x) identical(x$status, "ok"), TRUE)
  if (!any(ok)) stop("all events failed; first error: ", res[[1]]$status)
  horizon <- names(res[ok][[1]]$att)
  # every successful event must expose the same horizon labels: otherwise
  # do.call(rbind, ...) recycles a short/gapped vector into the wrong columns
  # (only a warning, easily lost under mclapply) and mixes effects across
  # different event days. With #1 fixed, a truncated event now fails its fit
  # and is dropped, but assert here so a mismatch is a hard error, not silent.
  bad <- Filter(function(x) !identical(names(x$att), horizon), res[ok])
  if (length(bad))
    stop("events have mismatched horizons; cannot combine: ",
         paste(vapply(bad, function(x) as.character(x$event), ""), collapse = ", "))
  atts <- do.call(rbind, lapply(res[ok], `[[`, "att"))
  cars <- do.call(rbind, lapply(res[ok], `[[`, "car"))
  colnames(atts) <- colnames(cars) <- horizon
  rownames(atts) <- rownames(cars) <- vapply(res[ok], function(x) as.character(x$event), "")

  ev_out <- data.frame(event = vapply(res, `[[`, ids[1], "event"),
                       n_treated = vapply(res, `[[`, NA_integer_, "n_treated"),
                       n_donors = vapply(res, `[[`, NA_integer_, "n_donors"),
                       pre_rmse = vapply(res, `[[`, NA_real_, "pre_rmse"),
                       att_avg_se = vapply(res, `[[`, NA_real_, "att_avg_se"),
                       status = vapply(res, `[[`, "", "status"))
  extra <- events[!duplicated(events$event),
                  setdiff(names(events), c("unit", "event_time")), drop = FALSE]
  ev_out <- merge(ev_out, extra, by = "event", sort = FALSE)

  ses <- NULL
  if (event_se != "none") {
    ses <- do.call(rbind, lapply(res[ok], function(x)
      if (is.null(x$se_att)) rep(NA_real_, length(horizon))
      else unname(x$se_att)))
    dimnames(ses) <- dimnames(atts)
  }

  n_ev <- sum(ok)
  # per-event weights: equal (NULL), the realized treated counts, or an
  # `events` column; equal weights keep the original colMeans / sd/sqrt(n)
  # code paths so long-standing results stay bit-identical, and the weighted
  # branch nests them (n/(n-1) * sum(wn^2 (x - xbar_w)^2) = var(x)/n at
  # wn = 1/n)
  w <- rep(1, n_ev)
  if (!is.null(weights)) {
    w <- if (identical(weights, "n_treated")) {
      as.numeric(vapply(res[ok], `[[`, NA_integer_, "n_treated"))
    } else {
      first <- !duplicated(events$event)
      wcol <- events[[weights]]
      if (is.null(wcol))
        stop("`weights` must be \"n_treated\" or a column of `events`")
      as.numeric(wcol[first][match(rownames(atts),
                                   as.character(events$event[first]))])
    }
    if (anyNA(w) || any(w <= 0))
      stop("event weights must be positive and known for every fitted event")
  }
  wn <- w / sum(w)
  pool <- function(M) if (is.null(weights)) colMeans(M)
          else stats::setNames(as.vector(wn %*% M), colnames(M))
  cross_se <- function(M) {
    if (n_ev < 2L)
      return(stats::setNames(rep(NA_real_, ncol(M)), colnames(M)))
    if (is.null(weights)) return(apply(M, 2, stats::sd) / sqrt(n_ev))
    ctr <- sweep(M, 2, as.vector(wn %*% M))
    stats::setNames(sqrt(colSums(ctr^2 * wn^2) * n_ev / (n_ev - 1)),
                    colnames(M))
  }
  se_out <- if (se == "cross") {
    # within-event component and its floor on the cross SE. The cross-event
    # sd already embeds each event's estimation noise under independence, so
    # adding the two would double count; `total_att` instead floors the
    # cross SE by the known within component (the DerSimonian-Laird
    # method-of-moments logic: total = sqrt(max(s^2 - m, 0) + m)/sqrt(n)
    # = max(cross, within) at equal weights), which binds when few events
    # make the cross dispersion estimate unreliably small.
    within_att <- if (!is.null(ses)) sqrt(colSums(ses^2 * wn^2))
    list(att = cross_se(atts), car = cross_se(cars),
         within_att = within_att,
         total_att = if (!is.null(within_att))
           pmax(cross_se(atts), within_att),
         method = "cross", n_events = n_ev)
  }

  structure(list(
    att = pool(atts), car = pool(cars), se = se_out,
    atts = atts, cars = cars, ses = ses, events = ev_out,
    weights = if (!is.null(weights)) stats::setNames(w, rownames(atts)),
    fits = if (keep_fits) lapply(res, `[[`, "fit"),
    method = method,
    conventions = list(returns = match.arg(returns, c("simple", "log")),
                       window = window, est_window = est_window),
    call = match.call()
  ), class = "fes_batch")
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' @export
print.fes_batch <- function(x, ...) {
  okn <- sum(x$events$status == "ok")
  cat("feventr batch fit: method '", x$method, "', ", okn, "/",
      nrow(x$events), " events\n", sep = "")
  # att names are the post-window offsets; window = c(1, 10) or a gapped
  # value-aligned panel need not contain "0", so pick the first available
  h0 <- if ("0" %in% names(x$att)) "0" else names(x$att)[1]
  cat("ATT (cross-event avg) at horizon ", h0, ": ",
      formatC(x$att[[h0]], digits = 4, format = "f"),
      if (!is.null(x$se)) paste0(" (se ",
                                 formatC(x$se$att[[h0]], digits = 4, format = "f"), ")"),
      "\n", sep = "")
  invisible(x)
}

#' Cross-event summary, optionally grouped
#'
#' @param object An `fes_batch`.
#' @param by Optional column of `object$events` (e.g. `"decade"`): grouped
#'   cross-event means and SEs of the ATT at `horizon`.
#' @param horizon Event time(s) to summarize; default all.
#' @param cumulative Summarize CARs instead of per-period ATTs.
#' @param ... Unused.
#' @export
summary.fes_batch <- function(object, by = NULL, horizon = NULL,
                              cumulative = FALSE, ...) {
  M <- if (cumulative) object$cars else object$atts
  if (!is.null(horizon)) M <- M[, as.character(horizon), drop = FALSE]
  ok <- object$events$status == "ok"
  wtd <- !is.null(object$weights)
  wstat <- function(Mi, wi) {
    if (!wtd)
      return(list(m = colMeans(Mi),
                  se = apply(Mi, 2, stats::sd) / sqrt(nrow(Mi))))
    wni <- wi / sum(wi)
    m <- as.vector(wni %*% Mi)
    se <- if (nrow(Mi) < 2L) rep(NA_real_, ncol(Mi))
          else sqrt(colSums(sweep(Mi, 2, m)^2 * wni^2) *
                      nrow(Mi) / (nrow(Mi) - 1))
    list(m = m, se = se)
  }
  w_all <- if (wtd) object$weights[rownames(M)]
  if (is.null(by)) {
    s <- wstat(M, w_all)
    data.frame(event_time = as.numeric(colnames(M)),
               estimate = s$m, se = s$se, row.names = NULL)
  } else {
    g <- object$events[[by]][ok][match(rownames(M),
                                       as.character(object$events$event[ok]))]
    do.call(rbind, lapply(split(seq_len(nrow(M)), g), function(i) {
      s <- wstat(M[i, , drop = FALSE], w_all[i])
      data.frame(group = g[i][1], event_time = as.numeric(colnames(M)),
                 estimate = s$m, se = s$se,
                 n_events = length(i), row.names = NULL)
    }))
  }
}

#' @export
coef.fes_batch <- function(object, cumulative = FALSE, ...) {
  if (cumulative) object$car else object$att
}
