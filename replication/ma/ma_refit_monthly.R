# M&A monthly refit — 3-year-horizon per-deal CATTs at MONTHLY frequency,
# from the WRDS CRSP monthly pull (wrds_monthly_pull.R). The monthly
# frequency bills the price-noise inflation 12x/year instead of 250x
# (channel 2 of MEMO_longrun_bias.md), at the cost of a short estimation
# window (T0 = 35 months) for the fit.
#
# Design mirrors ma_refit_full.R one level up: cohorts are announcement
# MONTHS (all that month's acquirers treated, same-month acquirers
# excluded from every donor pool), panels are event months -36..+36
# (complete 73-month presence required, NA returns zero-filled), fits are
# per deal with est_window = months -36..-2, onset buffer at -1,
# event_time 0 = the announcement month. CAR paths downstream are summed
# att from month 0, rebased at -1. Deals announced after 2021-12 fall
# out (no full +36m window by 2024-12).
#
# Modes: method from the command line (sc | gsynth). Env knobs as in the
# daily runner: MA_REFIT_PLACEBO=1|runup, MA_REFIT_DEMEAN=1,
# MA_REFIT_CORES, MA_REFIT_MAX_NEW, MA_REFIT_LIMIT. Checkpoints one CSV
# per (mode, announcement month) under ma/ma_refit_out/monthly_<dir>/.
#
# Run from replication/:
#   MA_REFIT_CORES=6 Rscript ma/ma_refit_monthly.R sc
source("config.R")
suppressMessages({library(feventr); library(haven); library(data.table)
                  library(parallel)})
ma_work <- function(...) dind("M&A", "data", "work", ...)

method <- commandArgs(trailingOnly = TRUE)
method <- if (length(method)) method[1] else "sc"
stopifnot(method %in% c("sc", "gsynth"))
placebo_env <- Sys.getenv("MA_REFIT_PLACEBO", "")
placebo <- nzchar(placebo_env) && placebo_env != "0"
match_runup <- identical(placebo_env, "runup")
demean <- nzchar(Sys.getenv("MA_REFIT_DEMEAN", ""))
# MA_REFIT_ABK=1 (method sc only): prior-gross-return weighted
# counterfactual from the fitted weights — purges transitory-noise
# inflation (channel-2 member i) but by construction NOT the real-vol
# convexity wedge (member ii); see MEMO_longrun_bias.md.
abk <- nzchar(Sys.getenv("MA_REFIT_ABK", ""))

cat("loading monthly CRSP ...\n")
mo <- fread(ma_work("crsp_monthly_wrds_1974_2024.csv.gz"))
mo[, permno := as.character(permno)]
months_all <- sort(unique(mo$month))
mo[, m_ix := match(month, months_all)]
setkey(mo, m_ix)

pdix <- as.data.table(read_dta(
  ma_work("sdc_ma_details_sl_m_cleaned_permno_dateindex_2023.dta")))
di <- as.data.table(read_dta(
  ma_work("sdc_ma_details_sl_m_cleaned_dateindex_2023.dta"),
  col_select = c("date", "date_index")))
di_dates <- unique(di[, .(ann_tdate = as.character(as.Date(date)),
                          date_index)])
deals <- merge(pdix[, .(permno = as.character(permno), date_index)],
               di_dates, by = "date_index")
deals[, ann_month := substr(ann_tdate, 1, 7)]
deals[, am_ix := match(ann_month, months_all)]
deals <- deals[!is.na(am_ix) & am_ix - 36L >= 1L &
                 am_ix + 36L <= length(months_all)]
cohorts <- sort(unique(deals$am_ix))
cat("announcement months:", length(cohorts), " deals in horizon:",
    nrow(deals), "of", nrow(pdix), "\n")

lim <- Sys.getenv("MA_REFIT_LIMIT", "")
if (nzchar(lim)) cohorts <- head(cohorts, as.integer(lim))

fit_cohort <- function(aix) {
  am <- months_all[aix]
  win <- seq.int(aix - 36L, aix + 36L)
  pan <- mo[.(win)]
  pan <- pan[, ok := .N == 73L, by = permno][ok == TRUE][, ok := NULL]
  pan[is.na(ret), ret := 0]
  pan[, event_month := m_ix - aix]
  if (demean)
    pan[, ret := ret - mean(ret[event_month >= -36L & event_month <= -2L]),
        by = permno]
  tr_all <- unique(deals[am_ix == aix, permno])
  tr_in <- intersect(tr_all, unique(pan$permno))
  donors <- setdiff(unique(pan$permno), tr_all)
  if (placebo) {
    set.seed(3000000L + aix)
    cand <- donors
    if (match_runup) {
      cr <- pan[event_month >= -36L & event_month <= -2L,
                .(cr = sum(ret)), by = permno]
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
      rows[[k]] <- data.table(ann_month = am, permno = tr,
                              matched_permno = matched[k],
                              status = "skipped: no 73-month coverage")
      next
    }
    res <- tryCatch({
      f <- event_study(pan, unit = "permno", time = "event_month",
                       ret = "ret", treated = tr, event_time = 0,
                       method = method, window = c(-1, 36),
                       est_window = c(-36, -2), returns = "simple",
                       cumulate = "log", se = "none", keep_data = FALSE,
                       donors = donors_fit)
      days_v <- as.integer(names(f$car))
      att_v <- as.numeric(f$att)
      car_v <- as.numeric(f$car)
      if (abk) {
        w <- f$weights$omega
        w <- w[w > 1e-10]
        dsub <- pan[permno %in% names(w) & event_month >= -2L,
                    .(permno, event_month, ret)]
        setorder(dsub, permno, event_month)
        dsub[, glag := 1 + data.table::shift(ret), by = permno]
        dsub <- dsub[event_month >= -1L]
        dsub[, wj := w[permno]]
        y0 <- dsub[, .(y0hat = sum(wj * glag * ret) / sum(wj * glag)),
                   by = event_month]
        m <- merge(pan[permno == tr & event_month >= -1L,
                       .(event_month, ret)], y0, by = "event_month")
        setorder(m, event_month)
        stopifnot(identical(as.integer(m$event_month), days_v))
        att_v <- m$ret - m$y0hat
        car_v <- cumsum(log1p(m$ret) - log1p(m$y0hat))
      }
      data.table(ann_month = am, permno = tr, matched_permno = matched[k],
                 event_month = days_v,
                 att = att_v, car_log = car_v,
                 n_units = length(donors_fit) + 1L,
                 n_treated_month = length(tr_all),
                 r = f$diagnostics$info$r %||% NA_integer_,
                 pre_rmse = f$diagnostics$info$pre_rmse %||% NA_real_,
                 status = "ok")
    }, error = function(e)
      data.table(ann_month = am, permno = tr, matched_permno = matched[k],
                 status = paste("error:", conditionMessage(e))))
    rows[[k]] <- res
  }
  rbindlist(rows, fill = TRUE)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

outdir <- file.path("ma", "ma_refit_out",
                    paste0("monthly_",
                           if (match_runup) "placebo_runup_"
                           else if (placebo) "placebo_" else "",
                           if (demean) "demean_" else "",
                           if (abk) "abk_" else "", method))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
have <- file.exists(file.path(outdir, sprintf("cohort_%s.csv",
                                              months_all[cohorts])))
todo <- cohorts[!have]
maxn <- Sys.getenv("MA_REFIT_MAX_NEW", "")
if (nzchar(maxn)) todo <- head(todo, as.integer(maxn))
cores <- as.integer(Sys.getenv("MA_REFIT_CORES", "6"))
cat("todo:", length(todo), "cohorts on", cores, "cores\n")
t0 <- Sys.time()
invisible(mclapply(todo, function(aix) {
  out <- tryCatch(fit_cohort(aix), error = function(e)
    data.table(ann_month = months_all[aix],
               status = paste("cohort error:", conditionMessage(e))))
  fwrite(out, file.path(outdir, sprintf("cohort_%s.csv",
                                        months_all[aix])))
  NULL
}, mc.cores = cores, mc.preschedule = FALSE))
cat(sprintf("%s: done %d | %.1f min\n", basename(outdir), length(todo),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
