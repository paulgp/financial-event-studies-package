# M&A full refit — per-deal acquirer CARs with feventr estimators.
#
# Rebuilds every announcement-date event panel from CRSP daily exactly as
# sdc_ma_malmendier_gsynth_batch.do (and table6.R Part C, where the recipe
# was validated against the saved gsynth output): calendar event days
# -280..+250 from ma_sl_cleaned_2023_event_date, donors = permnos with
# complete 531-day CRSP coverage, missing daret zero-filled, treated
# acquirers from sdc_ma_details_sl_m_cleaned_permno_dateindex_2023.
#
# Unlike the original (one joint multi-treated gsynth per date), each deal
# is fit separately — treated = that acquirer, donors exclude every
# acquirer announced the same day — which is what per-deal counterfactuals
# from sc/apm require. Panels are sliced in-memory from the keyed CRSP
# table (no 40GB cohort cache; forked workers share it copy-on-write).
#
# Fit: window c(-30, 250) (treatment onset -30, as the original),
# est_window c(-280, -31), returns "simple", cumulate "log", se "none".
# Per deal we store the published summary quantities: the [-1,+1] and
# [-1,+250] log CARs (car path differences), plus diagnostics.
#
# Checkpointed one CSV per (method, date_index) under ma/ma_refit_out/;
# reruns skip completed cohorts. Methods from the command line
# (default: sc apm).
#
# Run from replication/:
#   nohup Rscript ma/ma_refit_full.R sc apm > output/ma_refit_full.log 2>&1 &
source("config.R")
suppressMessages({library(feventr); library(haven); library(data.table)
                  library(parallel)})
ma_work <- function(...) dind("M&A", "data", "work", ...)

methods <- commandArgs(trailingOnly = TRUE)
if (!length(methods)) methods <- c("sc", "apm")
stopifnot(all(methods %in% c("sc", "apm", "gsynth")))
# MA_REFIT_PLACEBO=1: instead of the real acquirers, fit one date-matched
# placebo donor per deal (deterministic per-cohort draw, real acquirers
# excluded from pools), storing the matched real deal for subsample tags.
# Output goes to ma_refit_out/placebo_<method>/.
# MA_REFIT_PLACEBO=runup: additionally match each placebo to its acquirer
# by nearest estimation-window cumulative return (so placebos inherit the
# acquirers' selection on pre-announcement runup — separates the
# reversion-of-selection bias of intercept estimators from the mechanical
# noise bias the random placebo measures). Deterministic; output goes to
# ma_refit_out/placebo_runup_<method>/.
placebo_env <- Sys.getenv("MA_REFIT_PLACEBO", "")
placebo <- nzchar(placebo_env) && placebo_env != "0"
match_runup <- identical(placebo_env, "runup")
# MA_REFIT_DEMEAN=1: SC-with-intercept (Doudchenko-Imbens / Ferman-Pinto) —
# demean every unit's daret by its estimation-window mean before fitting.
# Kills any stationary unit-level mean gap (noise inflation AND alpha; the
# two are not separable here — see MEMO_longrun_bias.md). att stays a
# simple-return gap; car_log is on demeaned returns (do not interpret its
# level). Output prefix demean_.
demean <- nzchar(Sys.getenv("MA_REFIT_DEMEAN", ""))
# MA_REFIT_ABK=1: Asparouhova-Bessembinder-Kalcheva prior-gross-return
# weighted counterfactual — same fitted SC weights, but
# y0hat_t = sum_j w_j (1+r_{j,t-1}) r_{jt} / sum_j w_j (1+r_{j,t-1}),
# which purges each donor's price-noise inflation to first order while
# keeping the fixed-loading estimand (weights drift one day, not
# cumulatively). att/car_log recomputed from the ABK counterfactual.
# Output prefix abk_.
abk <- nzchar(Sys.getenv("MA_REFIT_ABK", ""))

cat("loading inputs (CRSP daily is 4GB; a few minutes) ...\n")
crsp <- as.data.table(read_dta(ma_work("crsp_daily_raw_2023_cleaned.dta"),
                               col_select = c("permno", "date", "daret")))
crsp[, date := as.Date(date)]
setkey(crsp, date)
cal <- as.data.table(read_dta(ma_work("ma_sl_cleaned_2023_event_date.dta"),
                              col_select = c("ann_tdate", "date", "event_date")))
cal[, `:=`(ann_tdate = as.character(as.Date(ann_tdate)), date = as.Date(date))]
cal <- cal[event_date >= -280 & event_date <= 250]
pdix <- as.data.table(read_dta(
  ma_work("sdc_ma_details_sl_m_cleaned_permno_dateindex_2023.dta")))
di <- as.data.table(read_dta(
  ma_work("sdc_ma_details_sl_m_cleaned_dateindex_2023.dta"),
  col_select = c("date", "date_index")))
di_dates <- unique(di[, .(ann_tdate = as.character(as.Date(date)), date_index)])

cohorts <- merge(unique(pdix[, .(date_index)]), di_dates, by = "date_index")
setorder(cohorts, date_index)
cat("cohorts:", nrow(cohorts), " deals:", nrow(pdix), "\n")

# pilot mode: MA_REFIT_LIMIT=n restricts to the first n-1 cohorts plus the
# first multi-deal cohort; checkpoints land in the real output dirs, so the
# full run reuses them
lim <- Sys.getenv("MA_REFIT_LIMIT", "")
if (nzchar(lim)) {
  multi <- pdix[, .N, by = date_index][N > 1L, date_index]
  keep <- unique(c(head(cohorts$date_index, as.integer(lim) - 1L),
                   head(intersect(cohorts$date_index, multi), 1L)))
  cohorts <- cohorts[date_index %in% keep]
  cat("PILOT: restricted to", nrow(cohorts), "cohorts:",
      paste(cohorts$date_index, collapse = ", "), "\n")
}

fit_cohort <- function(ix, ann, method) {
  days <- cal[ann_tdate == ann]
  if (nrow(days) < 531L)
    return(data.table(date_index = ix, ann_tdate = ann, permno = NA_character_,
                      status = sprintf("skipped: calendar %d < 531 days",
                                       nrow(days))))
  setorder(days, event_date)
  pan <- crsp[.(days$date), nomatch = 0L]
  pan <- pan[, ok := .N == 531L, by = permno][ok == TRUE][, ok := NULL]
  pan <- merge(pan, days[, .(date, event_date)], by = "date")
  pan[is.na(daret), daret := 0]
  pan[, permno := as.character(permno)]
  if (demean)
    pan[, daret := daret -
          mean(daret[event_date >= -280L & event_date <= -31L]),
        by = permno]
  tr_all <- as.character(pdix[date_index == ix, permno])
  tr_in <- intersect(tr_all, unique(pan$permno))
  donors <- setdiff(unique(pan$permno), tr_all)
  if (placebo) {
    set.seed(1000000L + ix)
    cand <- setdiff(unique(pan$permno), tr_all)
    if (match_runup) {
      # one placebo per acquirer WITH complete coverage (its runup must be
      # computable), nearest-neighbor on est-window cumulative return,
      # without replacement within the cohort
      cr <- pan[event_date >= -280L & event_date <= -31L,
                .(cr = sum(daret)), by = permno]
      crv <- setNames(cr$cr, cr$permno)
      fit_units <- character(0)
      matched <- character(0)
      pool <- cand
      for (tr0 in tr_in) {
        if (!length(pool)) break
        pick <- pool[which.min(abs(crv[pool] - crv[[tr0]]))]
        fit_units <- c(fit_units, pick)
        matched <- c(matched, tr0)
        pool <- setdiff(pool, pick)
      }
    } else {
      fit_units <- sample(cand, min(length(tr_all), length(cand)))
      matched <- tr_all[seq_along(fit_units)]
    }
    donors_fit <- setdiff(donors, fit_units)
  } else {
    fit_units <- tr_all
    matched <- tr_all
    donors_fit <- donors
  }
  rows <- list()
  for (k in seq_along(fit_units)) {
    tr <- fit_units[k]
    if (!placebo && !tr %in% tr_in) {
      rows[[k]] <- data.table(date_index = ix, ann_tdate = ann, permno = tr,
                              matched_permno = matched[k],
                              status = "skipped: no complete CRSP coverage")
      next
    }
    res <- tryCatch({
      f <- event_study(pan, unit = "permno", time = "event_date",
                       ret = "daret", treated = tr, event_time = 0,
                       method = method, window = c(-30, 250),
                       est_window = c(-280, -31), returns = "simple",
                       cumulate = "log", se = "none", keep_data = FALSE,
                       donors = donors_fit)
      # full per-deal path (one row per event day, like the original
      # stacked gsynth output): att and the cumulative log CAR from -30;
      # horizon CARs from -1 are path differences vs event_date == -2
      days_v <- as.integer(names(f$car))
      att_v <- as.numeric(f$att)
      car_v <- as.numeric(f$car)
      if (abk) {
        w <- f$weights$omega
        w <- w[w > 1e-10]
        dsub <- pan[permno %in% names(w) & event_date >= -31L,
                    .(permno, event_date, daret)]
        setorder(dsub, permno, event_date)
        dsub[, glag := 1 + data.table::shift(daret), by = permno]
        dsub <- dsub[event_date >= -30L]
        dsub[, wj := w[permno]]
        y0 <- dsub[, .(y0hat = sum(wj * glag * daret) / sum(wj * glag)),
                   by = event_date]
        m <- merge(pan[permno == tr & event_date >= -30L,
                       .(event_date, daret)], y0, by = "event_date")
        setorder(m, event_date)
        stopifnot(identical(as.integer(m$event_date), days_v))
        att_v <- m$daret - m$y0hat
        car_v <- cumsum(log1p(m$daret) - log1p(m$y0hat))
      }
      data.table(date_index = ix, ann_tdate = ann, permno = tr,
                 matched_permno = matched[k],
                 event_date = days_v,
                 att = att_v, car_log = car_v,
                 n_units = length(donors_fit) + 1L,
                 n_treated_date = length(tr_all),
                 r = f$diagnostics$info$r %||% NA_integer_,
                 pre_rmse = f$diagnostics$info$pre_rmse %||% NA_real_,
                 status = "ok")
    }, error = function(e)
      data.table(date_index = ix, ann_tdate = ann, permno = tr,
                 matched_permno = matched[k],
                 status = paste0("failed: ", conditionMessage(e))))
    rows[[k]] <- res
  }
  rbindlist(rows, fill = TRUE)
}

`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a

for (method in methods) {
  outdir <- file.path("ma", "ma_refit_out",
                      paste0(if (match_runup) "placebo_runup_"
                             else if (placebo) "placebo_" else "",
                             if (demean) "demean_" else "",
                             if (abk) "abk_" else "", method))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  # MA_REFIT_MAX_NEW=n caps the number of not-yet-checkpointed cohorts this
  # invocation fits, so a run can be sized to finish inside a time budget
  # and chained; completed cohorts are skipped as usual
  max_new <- Sys.getenv("MA_REFIT_MAX_NEW", "")
  work <- cohorts
  if (nzchar(max_new)) {
    have <- file.exists(file.path(outdir, sprintf("cohort_%d.csv",
                                                  cohorts$date_index)))
    work <- cohorts[!have][seq_len(min(as.integer(max_new), sum(!have)))]
    cat(sprintf("%s: %d cohorts checkpointed, fitting next %d\n",
                method, sum(have), nrow(work)))
  }
  t0 <- proc.time()[3]
  st <- mclapply(seq_len(nrow(work)), function(k) {
    setDTthreads(1L)
    ix <- work$date_index[k]
    fn <- file.path(outdir, sprintf("cohort_%d.csv", ix))
    if (file.exists(fn)) return("skip")
    res <- tryCatch(fit_cohort(ix, work$ann_tdate[k], method),
                    error = function(e) e)
    if (inherits(res, "error")) return(paste0("failed: ",
                                              conditionMessage(res)))
    tmp <- paste0(fn, ".tmp")
    fwrite(res, tmp)
    file.rename(tmp, fn)
    "done"
  }, mc.cores = as.integer(Sys.getenv("MA_REFIT_CORES", "6")),
     mc.preschedule = FALSE)
  nulls <- sum(vapply(st, is.null, TRUE))
  cat(sprintf("%s: %s%s | %.1f min\n", method,
              paste(names(table(unlist(st))), table(unlist(st)),
                    collapse = ", "),
              if (nulls) paste0(", killed-fork ", nulls) else "",
              (proc.time()[3] - t0) / 60))
}
