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
  tr_all <- as.character(pdix[date_index == ix, permno])
  tr_in <- intersect(tr_all, unique(pan$permno))
  donors <- setdiff(unique(pan$permno), tr_all)
  rows <- list()
  for (tr in tr_all) {
    if (!tr %in% tr_in) {
      rows[[tr]] <- data.table(date_index = ix, ann_tdate = ann, permno = tr,
                               status = "skipped: no complete CRSP coverage")
      next
    }
    res <- tryCatch({
      f <- event_study(pan, unit = "permno", time = "event_date",
                       ret = "daret", treated = tr, event_time = 0,
                       method = method, window = c(-30, 250),
                       est_window = c(-280, -31), returns = "simple",
                       cumulate = "log", se = "none", keep_data = FALSE)
      ed <- as.integer(names(f$car))
      car_at <- function(h) unname(f$car[ed == h])
      data.table(date_index = ix, ann_tdate = ann, permno = tr,
                 n_units = length(donors) + 1L,
                 n_treated_date = length(tr_all),
                 r = f$diagnostics$info$r %||% NA_integer_,
                 pre_rmse = f$diagnostics$info$pre_rmse %||% NA_real_,
                 car_log_1 = car_at(1L) - car_at(-2L),
                 car_log_250 = car_at(250L) - car_at(-2L),
                 status = "ok")
    }, error = function(e)
      data.table(date_index = ix, ann_tdate = ann, permno = tr,
                 status = paste0("failed: ", conditionMessage(e))))
    rows[[tr]] <- res
  }
  rbindlist(rows, fill = TRUE)
}

`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a

for (method in methods) {
  outdir <- file.path("ma", "ma_refit_out", method)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  t0 <- proc.time()[3]
  st <- mclapply(seq_len(nrow(cohorts)), function(k) {
    setDTthreads(1L)
    ix <- cohorts$date_index[k]
    fn <- file.path(outdir, sprintf("cohort_%d.csv", ix))
    if (file.exists(fn)) return("skip")
    res <- tryCatch(fit_cohort(ix, cohorts$ann_tdate[k], method),
                    error = function(e) e)
    if (inherits(res, "error")) return(paste0("failed: ",
                                              conditionMessage(res)))
    tmp <- paste0(fn, ".tmp")
    fwrite(res, tmp)
    file.rename(tmp, fn)
    "done"
  }, mc.cores = 6L, mc.preschedule = FALSE)
  nulls <- sum(vapply(st, is.null, TRUE))
  cat(sprintf("%s: %s%s | %.1f min\n", method,
              paste(names(table(unlist(st))), table(unlist(st)),
                    collapse = ", "),
              if (nulls) paste0(", killed-fork ", nulls) else "",
              (proc.time()[3] - t0) / 60))
}
