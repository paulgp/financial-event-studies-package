# Hybrid M&A refit — SC weights fit on DAILY pre-period data (T0 = 250
# obs, well-identified loadings), counterfactual evaluated on MONTHLY
# returns out to +36 months (12 price touches/year instead of 250, so
# channel 2 of MEMO_longrun_bias.md is billed ~21x less often). The
# fitting step still selects donors on daily noise (piece A), so the
# per-touch differential survives once per month — the placebo columns
# measure exactly that residual.
#
# Cohorts are announcement DATES (the daily est window -280..-31 is
# date-specific); donors must pass BOTH screens (complete daily presence
# over the est window and complete 73-month presence around the
# announcement month) and exclude every same-MONTH acquirer. Weights via
# feventr::solve_simplex_ls on the daily est window; per-deal monthly
# att path stored for event months -36..+36 (the out-of-sample pre
# months are a fit diagnostic). Placebo knobs as in the daily runner
# (MA_REFIT_PLACEBO=1|runup; runup matches on daily est-window cumret).
#
# Checkpoints: ma/ma_refit_out/hybrid_<prefix>sc/cohort_<date_index>.csv
# Run from replication/: MA_REFIT_CORES=6 Rscript ma/ma_refit_hybrid.R
source("config.R")
suppressMessages({library(feventr); library(haven); library(data.table)
                  library(parallel)})
ma_work <- function(...) dind("M&A", "data", "work", ...)

placebo_env <- Sys.getenv("MA_REFIT_PLACEBO", "")
placebo <- nzchar(placebo_env) && placebo_env != "0"
match_runup <- identical(placebo_env, "runup")

cat("loading inputs (daily CRSP is 4GB; a few minutes) ...\n")
crsp <- as.data.table(read_dta(ma_work("crsp_daily_raw_2023_cleaned.dta"),
                               col_select = c("permno", "date", "daret")))
crsp[, date := as.Date(date)]
setkey(crsp, date)
mo <- fread(ma_work("crsp_monthly_wrds_1974_2024.csv.gz"))
mo[, permno := as.character(permno)]
months_all <- sort(unique(mo$month))
mo[, m_ix := match(month, months_all)]
setkey(mo, m_ix)
cal <- as.data.table(read_dta(ma_work("ma_sl_cleaned_2023_event_date.dta"),
                              col_select = c("ann_tdate", "date",
                                             "event_date")))
cal[, `:=`(ann_tdate = as.character(as.Date(ann_tdate)),
           date = as.Date(date))]
cal <- cal[event_date >= -280 & event_date <= -31]
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
cohorts <- unique(deals[, .(date_index, ann_tdate, am_ix)])
setorder(cohorts, date_index)
cat("cohorts:", nrow(cohorts), " deals:", nrow(deals), "\n")

lim <- Sys.getenv("MA_REFIT_LIMIT", "")
if (nzchar(lim)) cohorts <- head(cohorts, as.integer(lim))

fit_cohort <- function(ix, ann, aix) {
  days <- cal[ann_tdate == ann]
  if (nrow(days) < 250L)
    return(data.table(date_index = ix, ann_month = months_all[aix],
                      permno = NA_character_,
                      status = sprintf("skipped: est calendar %d < 250",
                                       nrow(days))))
  setorder(days, event_date)
  dpan <- crsp[.(days$date), nomatch = 0L]
  dpan <- dpan[, ok := .N == nrow(days), by = permno][ok == TRUE][
    , ok := NULL]
  dpan[is.na(daret), daret := 0]
  dpan[, permno := as.character(permno)]
  win <- seq.int(aix - 36L, aix + 36L)
  mpan <- mo[.(win)]
  mpan <- mpan[, ok := .N == 73L, by = permno][ok == TRUE][, ok := NULL]
  mpan[is.na(ret), ret := 0]
  mpan[, event_month := m_ix - aix]
  both <- intersect(unique(dpan$permno), unique(mpan$permno))
  tr_month <- unique(deals[am_ix == aix, permno])   # same-MONTH acquirers
  tr_all <- unique(deals[date_index == ix, permno])
  tr_in <- intersect(tr_all, both)
  donors <- setdiff(both, tr_month)
  if (placebo) {
    set.seed(4000000L + ix)
    cand <- donors
    if (match_runup) {
      cr <- dpan[permno %in% both, .(cr = sum(daret)), by = permno]
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
  # daily est-window return matrix, donors x days
  dsub <- dpan[permno %in% c(donors_fit, fit_units)]
  D <- dcast(dsub, permno ~ date, value.var = "daret")
  dn <- D$permno
  D <- as.matrix(D[, -1])
  rownames(D) <- dn
  Ad <- t(D[donors_fit, , drop = FALSE])
  # monthly return matrix, units x event months
  msub <- mpan[permno %in% c(donors_fit, fit_units)]
  M <- dcast(msub, permno ~ event_month, value.var = "ret")
  mn <- M$permno
  M <- as.matrix(M[, -1])
  rownames(M) <- mn
  emonths <- as.integer(colnames(M))
  rows <- list()
  for (k in seq_along(fit_units)) {
    tr <- fit_units[k]
    if (!tr %in% rownames(D) || !tr %in% rownames(M)) {
      rows[[k]] <- data.table(date_index = ix,
                              ann_month = months_all[aix], permno = tr,
                              matched_permno = matched[k],
                              status = "skipped: no dual coverage")
      next
    }
    res <- tryCatch({
      sol <- solve_simplex_ls(Ad, D[tr, ])
      w <- sol$w
      y0 <- as.vector(crossprod(M[donors_fit, , drop = FALSE], w))
      att <- M[tr, ] - y0
      data.table(date_index = ix, ann_month = months_all[aix],
                 permno = tr, matched_permno = matched[k],
                 event_month = emonths, att = as.numeric(att),
                 car_log = cumsum(log1p(M[tr, ]) - log1p(y0)),
                 n_units = length(donors_fit) + 1L,
                 n_treated_month = length(tr_month),
                 pre_rmse = sqrt(mean((D[tr, ] -
                                         as.vector(Ad %*% w))^2)),
                 status = "ok")
    }, error = function(e)
      data.table(date_index = ix, ann_month = months_all[aix],
                 permno = tr, matched_permno = matched[k],
                 status = paste("error:", conditionMessage(e))))
    rows[[k]] <- res
  }
  rbindlist(rows, fill = TRUE)
}

outdir <- file.path("ma", "ma_refit_out",
                    paste0("hybrid_",
                           if (match_runup) "placebo_runup_"
                           else if (placebo) "placebo_" else "", "sc"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
have <- file.exists(file.path(outdir, sprintf("cohort_%d.csv",
                                              cohorts$date_index)))
todo <- cohorts[!have]
maxn <- Sys.getenv("MA_REFIT_MAX_NEW", "")
if (nzchar(maxn)) todo <- head(todo, as.integer(maxn))
cores <- as.integer(Sys.getenv("MA_REFIT_CORES", "6"))
cat("todo:", nrow(todo), "cohorts on", cores, "cores\n")
t0 <- Sys.time()
invisible(mclapply(seq_len(nrow(todo)), function(i) {
  out <- tryCatch(
    fit_cohort(todo$date_index[i], todo$ann_tdate[i], todo$am_ix[i]),
    error = function(e)
      data.table(date_index = todo$date_index[i],
                 status = paste("cohort error:", conditionMessage(e))))
  fwrite(out, file.path(outdir, sprintf("cohort_%d.csv",
                                        todo$date_index[i])))
  NULL
}, mc.cores = cores, mc.preschedule = FALSE))
cat(sprintf("%s: done %d | %.1f min\n", basename(outdir), nrow(todo),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
