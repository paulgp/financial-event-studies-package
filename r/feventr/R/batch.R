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
#' @param se `"cross"` (default): cross-event mean and `sd/sqrt(n)` per
#'   horizon; `"none"`.
#' @param cores Events are fit with `parallel::mclapply(mc.cores = cores)`.
#' @param keep_fits Store the full per-event `fes_fit` list (memory-heavy).
#' @return An object of class `fes_batch`: `$att`/`$car` (cross-event mean
#'   paths), `$se`, `$atts`/`$cars` (event x horizon matrices), `$events`
#'   (input table + n_treated, n_donors, pre_rmse, status), `$fits` if kept.
#' @export
event_study_batch <- function(data, unit, time, ret, events,
                              method = c("mean", "did", "market", "factor",
                                         "sc", "ridge", "sdid", "gsynth"),
                              window = c(0, 10), est_window = c(-250, -11),
                              returns, cumulate = c("auto", "sum", "compound", "log"),
                              align = c("position", "value"),
                              factors = NULL, exclude_treated = TRUE,
                              match_on = c("ret", "cumret"), V = NULL,
                              solver = c("hybrid", "fw", "qp"), lambda = NULL,
                              r = c(0, 5), force = c("unit", "none", "two-way"),
                              se = c("cross", "none"), cores = 1L,
                              keep_fits = FALSE) {
  method <- match.arg(method)
  se <- match.arg(se)
  align <- match.arg(align)
  events <- as.data.frame(events)
  stopifnot(all(c("unit", "event_time") %in% names(events)))
  if (is.null(events$event))
    events$event <- seq_len(nrow(events))
  ids <- unique(events$event)

  # offsets between any two event times, for donor exclusion
  times <- sort(unique(data[[time]]))
  ev_pos <- if (align == "value") as.numeric(events$event_time)
            else match(events$event_time, times)
  if (anyNA(ev_pos)) stop("some `event_time`s not found among panel times")

  fit_one <- function(eid) {
    rows <- events$event == eid
    s <- ev_pos[rows][1]
    donors <- NULL
    if (exclude_treated) {
      contam <- ev_pos >= s + est_window[1] & ev_pos <= s + window[2]
      bad <- unique(as.character(events$unit[contam]))
      donors <- setdiff(unique(as.character(data[[unit]])), bad)
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
                       se = "none", keep_data = FALSE)
      list(event = eid, att = f$att, car = f$car,
           n_treated = f$diagnostics$n_treated,
           n_donors = f$diagnostics$n_donors,
           pre_rmse = f$diagnostics$info$pre_rmse %||% NA_real_,
           status = "ok", fit = if (keep_fits) f)
    }, error = function(e)
      list(event = eid, att = NULL, car = NULL, n_treated = NA_integer_,
           n_donors = NA_integer_, pre_rmse = NA_real_,
           status = paste0("dropped: ", conditionMessage(e)), fit = NULL))
  }

  res <- if (cores > 1L) parallel::mclapply(ids, fit_one, mc.cores = cores)
         else lapply(ids, fit_one)

  ok <- vapply(res, function(x) identical(x$status, "ok"), TRUE)
  if (!any(ok)) stop("all events failed; first error: ", res[[1]]$status)
  horizon <- names(res[ok][[1]]$att)
  atts <- do.call(rbind, lapply(res[ok], `[[`, "att"))
  cars <- do.call(rbind, lapply(res[ok], `[[`, "car"))
  rownames(atts) <- rownames(cars) <- vapply(res[ok], function(x) as.character(x$event), "")

  ev_out <- data.frame(event = vapply(res, `[[`, ids[1], "event"),
                       n_treated = vapply(res, `[[`, NA_integer_, "n_treated"),
                       n_donors = vapply(res, `[[`, NA_integer_, "n_donors"),
                       pre_rmse = vapply(res, `[[`, NA_real_, "pre_rmse"),
                       status = vapply(res, `[[`, "", "status"))
  extra <- events[!duplicated(events$event),
                  setdiff(names(events), c("unit", "event_time")), drop = FALSE]
  ev_out <- merge(ev_out, extra, by = "event", sort = FALSE)

  n_ev <- sum(ok)
  se_out <- if (se == "cross")
    list(att = apply(atts, 2, stats::sd) / sqrt(n_ev),
         car = apply(cars, 2, stats::sd) / sqrt(n_ev),
         method = "cross", n_events = n_ev)

  structure(list(
    att = colMeans(atts), car = colMeans(cars), se = se_out,
    atts = atts, cars = cars, events = ev_out,
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
  cat("ATT (cross-event avg) at horizon 0: ",
      formatC(x$att[["0"]], digits = 4, format = "f"),
      if (!is.null(x$se)) paste0(" (se ",
                                 formatC(x$se$att[["0"]], digits = 4, format = "f"), ")"),
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
  if (is.null(by)) {
    data.frame(event_time = as.numeric(colnames(M)),
               estimate = colMeans(M),
               se = apply(M, 2, stats::sd) / sqrt(nrow(M)), row.names = NULL)
  } else {
    g <- object$events[[by]][ok][match(rownames(M),
                                       as.character(object$events$event[ok]))]
    do.call(rbind, lapply(split(seq_len(nrow(M)), g), function(i) {
      data.frame(group = g[i][1], event_time = as.numeric(colnames(M)),
                 estimate = colMeans(M[i, , drop = FALSE]),
                 se = apply(M[i, , drop = FALSE], 2, stats::sd) / sqrt(length(i)),
                 n_events = length(i), row.names = NULL)
    }))
  }
}

#' @export
coef.fes_batch <- function(object, cumulative = FALSE, ...) {
  if (cumulative) object$car else object$att
}
