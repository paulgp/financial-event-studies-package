# Table 6 — M&A acquirer 3-day CARs (Goldsmith-Pinkham & Lyu, p.49)
#
# Published cells = cross-deal means of LOG CARs over event days [-1,+1]:
#   sum_{t=-1..1} [log(1+r_treated,t) - log(1+r_counterfactual,t)]
# market counterfactual = CRSP vwretd; gsynth counterfactual = per-deal Y.ct
# from gsynth(force="unit", r=c(1,10)), pre [-280,-31], treatment onset -30.
# Subsamples: full (14,847), public (3,297), private (7,030), other/Sub.
# (4,520), 100% cash (9,261), 100% stock (5,592). Targets are in PERCENT.
#
# Source pipeline (verified): M&A/code/sdc_ma_malmendier_gsynth.do (~400-573)
# -> output/car_log_sc_summ_stats.tex; market column from
# sdc_ma_malmendier_canonical.do -> output/sl_m_deals_car_1_250_mkt.dta.
#
# Compute strategy (documented deviations, accepted by orchestrator):
# - Raw per-date event panels (event_panel_i.dta) are not in Dropbox and
#   re-fitting 7,052 gsynth panels is multi-day compute. The gsynth column is
#   therefore RECOMPUTED (the [-1,+1] log-CAR aggregation) from the stacked
#   day-level gsynth output sc_ma_1_7052.dta (daret_treated / daret_sc per
#   deal-day) rather than trusting the saved deal-level CARs, and the package
#   gsynth path is validated LIVE on a small sample of announcement dates with
#   panels rebuilt from CRSP daily (Part C below).
# - Daily CRSP vwretd is NOT in Dropbox (only monthly indexes; the RA's $ccm
#   path held the daily index file), so the market column consumes the saved
#   per-deal market CARs (sl_m_deals_car_1_250_mkt.dta). Part B rebuilds the
#   treated-stock side of the market CAR from CRSP daily from first principles
#   and checks the implied 3-day market log-return is constant across deals
#   sharing an announcement date (it must be if the saved file is internally
#   consistent), which validates everything except the index series itself.
#
# Run from replication/:  Rscript ma/table6.R

source("config.R")
suppressMessages({
  library(haven)
  library(data.table)
})

ma_work <- function(...) dind("M&A", "data", "work", ...)
ma_out  <- function(...) dind("M&A", "output", ...)

msg <- function(...) cat(sprintf(...), "\n")

## ---------------------------------------------------------------------------
## Part A — Table 6 from the stacked day-level gsynth output
## ---------------------------------------------------------------------------

msg("Part A: reading stacked gsynth output (sc_ma_1_7052.dta) ...")
sc <- as.data.table(read_dta(
  ma_work("sdc_ma_sl_m_gsynth", "sc_ma_1_7052.dta"),
  col_select = c("permno", "event_date", "daret_treated", "daret_sc",
                 "ann_tdate")))

## 3-day log CAR per gsynth run (permno x ann_tdate), event days -1..+1.
## Stata: cum_log_treat - cum_log_sc evaluated at the 3rd row of the
## event_date >= -1 panel, i.e. the sum over event_date in {-1,0,1}.
car_sc <- sc[event_date %in% c(-1, 0, 1),
             .(car_log_sc_1 = sum(log1p(daret_treated)) - sum(log1p(daret_sc)),
               n_days = .N),
             by = .(permno, ann_tdate)]
stopifnot(all(car_sc$n_days == 3L))
car_sc[, n_days := NULL]
car_sc[, date := as.Date(ann_tdate)]
msg("  gsynth runs (permno x ann_tdate): %d over %d announcement dates",
    nrow(car_sc), uniqueN(car_sc$ann_tdate))

## Expand to deals: merge 1:m on permno + announcement trading day.
di <- as.data.table(read_dta(ma_work("sdc_ma_details_sl_m_cleaned_dateindex_2023.dta")))
di[, date := as.Date(date)]
deals <- merge(car_sc, di, by = c("permno", "date"))

## Market CARs (saved deal-level file; see header for why) by deal group.
mkt <- as.data.table(read_dta(
  ma_out("sl_m_deals_car_1_250_mkt.dta"),
  col_select = c("group", "car_log_vwretd_1")))
deals <- merge(deals, mkt, by = "group")

## Deal characteristics.
det <- as.data.table(read_dta(
  ma_work("sdc_ma_details_sl_m_cleaned_2023.dta"),
  col_select = c("master_deal_no", "permno", "tpublic", "pct_cash", "pct_stk")))
deals <- merge(deals, det, by = c("master_deal_no", "permno"))

msg("  deal-level sample: N = %d (published 14,847)", nrow(deals))
stopifnot(!anyDuplicated(deals$group))

## Cross-check the recomputed gsynth CARs against the saved deal-level file.
saved <- as.data.table(read_dta(
  ma_out("sl_m_deals_car_1_250_gsynth.dta"),
  col_select = c("group", "car_log_sc_1", "car_log_vwretd_1")))
chk <- merge(deals[, .(group, car_log_sc_1, car_log_vwretd_1)],
             saved, by = "group", suffixes = c("_recomputed", "_saved"))
msg("  recomputed vs saved car_log_sc_1: matched %d/%d, max |diff| = %.3g",
    nrow(chk), nrow(deals),
    chk[, max(abs(car_log_sc_1_recomputed - car_log_sc_1_saved))])

## Subsample means (in percent).
cells <- function(d) {
  sub <- list(
    "Full sample"     = rep(TRUE, nrow(d)),
    "Public targets"  = d$tpublic == "Public",
    "Private targets" = d$tpublic == "Priv.",
    "Other targets"   = d$tpublic == "Sub.",
    "Cash merger"     = d$pct_cash == 100,
    "Stock merger"    = d$pct_stk == 100)
  rbindlist(lapply(names(sub), function(s) {
    i <- which(sub[[s]])
    data.table(
      row = c("Market mean", "Gsynth mean"),
      col = s,
      estimate_pct = 100 * c(mean(d$car_log_vwretd_1[i]),
                             mean(d$car_log_sc_1[i])),
      n = length(i))
  }))
}
tab6 <- cells(deals)

write.csv(tab6, out_path("table6.csv"), row.names = FALSE)
msg("  wrote %s", out_path("table6.csv"))

## Comparison against targets.
tgt <- as.data.table(read_target(6))
tgt[, n_target := as.integer(gsub(",", "", n))]
cmp <- merge(tgt[, .(row, col, target_pct = estimate, n_target)],
             tab6, by = c("row", "col"))
cmp[, diff_pp := estimate_pct - target_pct]
## published values rounded to 0.1pp; tolerance +-0.1pp on the point
cmp[, ok_estimate := abs(diff_pp) <= 0.1]
cmp[, ok_n := n == n_target]
setorder(cmp, -row, col)  # Market rows then Gsynth rows
write.csv(cmp[, .(row, col, target_pct, estimate_pct, diff_pp,
                  n_target, n, ok_estimate, ok_n)],
          out_path("table6_comparison.csv"), row.names = FALSE)
msg("  wrote %s", out_path("table6_comparison.csv"))
print(cmp[, .(row, col, target_pct, estimate_pct = round(estimate_pct, 3),
              n_target, n, ok_estimate, ok_n)])
msg("Part A: %d/%d cells within tolerance, %d/%d Ns exact",
    sum(cmp$ok_estimate), nrow(cmp), sum(cmp$ok_n), nrow(cmp))

## ---------------------------------------------------------------------------
## Part B — market column: treated-side rebuild from CRSP daily
## ---------------------------------------------------------------------------
## car_log_vwretd_1 = sum_{-1..1} log(1+daret) - sum_{-1..1} log(1+vwretd)
## (Stata sum() skips missing daret == zero-fill). We rebuild the first term
## from crsp_daily_raw_2023_cleaned.dta and check that the implied market
## term (treated logsum - saved CAR) is constant across deals within an
## announcement date.

msg("Part B: reading CRSP daily (4.3GB; this takes a few minutes) ...")
crsp <- as.data.table(read_dta(
  ma_work("crsp_daily_raw_2023_cleaned.dta"),
  col_select = c("permno", "date", "daret")))
crsp[, date := as.Date(date)]
setkey(crsp, permno, date)

cal <- as.data.table(read_dta(ma_work("ma_sl_cleaned_2023_event_date.dta")))
cal[, `:=`(ann_tdate = as.character(as.Date(ann_tdate)), date = as.Date(date))]

ev3 <- cal[event_date %in% c(-1, 0, 1)]
need <- merge(deals[, .(group, ann_tdate, permno)], ev3,
              by = "ann_tdate", allow.cartesian = TRUE)
need <- merge(need, crsp, by = c("permno", "date"), all.x = TRUE)
treat_side <- need[, .(logsum_treat = sum(log1p(fcoalesce(daret, 0)))),
                   by = .(group, ann_tdate)]
treat_side <- merge(treat_side, deals[, .(group, car_log_vwretd_1)],
                    by = "group")
treat_side[, implied_mkt := logsum_treat - car_log_vwretd_1]
bydate <- treat_side[, .(spread = max(implied_mkt) - min(implied_mkt),
                         n = .N), by = ann_tdate]
msg("  implied 3-day market log-return spread within announcement date:")
msg("    max spread = %.3g across %d dates (%d dates with >1 deal)",
    bydate[, max(spread)], nrow(bydate), bydate[n > 1L, .N])
msg("  treated-side rebuild consistent with saved market CARs: %s",
    ifelse(bydate[, max(spread)] < 1e-6, "YES", "NO -- investigate"))

rm(need, ev3, treat_side); invisible(gc())

## ---------------------------------------------------------------------------
## Part C — live validation of feventr's gsynth path on sample dates
## ---------------------------------------------------------------------------
## Rebuild event panels from CRSP daily exactly as
## sdc_ma_malmendier_gsynth_batch.do (calendar dates event_date -280..+250,
## donors = permnos with complete 531-day CRSP coverage, missing daret
## zero-filled, treated from sdc_ma_details_sl_m_cleaned_permno_dateindex),
## fit feventr::event_study(method = "gsynth", force = "unit", r = c(1,10),
## window = c(-30,250), est_window = c(-280,-31)) and compare the per-day
## counterfactual to the saved daret_sc.

msg("Part C: live gsynth validation via feventr")
suppressMessages(library(feventr))

pdix <- as.data.table(read_dta(
  ma_work("sdc_ma_details_sl_m_cleaned_permno_dateindex_2023.dta")))
di_dates <- unique(di[, .(ann_tdate = as.character(date), date_index)])

## candidate dates: exactly one treated permno, present in the stacked file
n_tr <- pdix[, .(n_treated = uniqueN(permno)), by = date_index]
in_sc <- merge(unique(car_sc[, .(ann_tdate)]), di_dates, by = "ann_tdate")
cand <- merge(in_sc, n_tr[n_treated == 1L], by = "date_index")
setorder(cand, date_index)
pick <- cand[round(quantile(seq_len(.N), c(0.10, 0.35, 0.60, 0.85)))]
msg("  validation dates (date_index / ann_tdate): %s",
    paste(sprintf("%d/%s", pick$date_index, pick$ann_tdate), collapse = ", "))

val <- list()
for (k in seq_len(nrow(pick))) {
  ix <- pick$date_index[k]; ann <- pick$ann_tdate[k]
  days <- cal[ann_tdate == ann & event_date >= -280 & event_date <= 250]
  setorder(days, event_date)
  stopifnot(nrow(days) == 531L)

  pan <- merge(days[, .(date, event_date)], crsp, by = "date")
  pan <- pan[, ok := .N == 531L, by = permno][ok == TRUE][, ok := NULL]
  pan[is.na(daret), daret := 0]
  pan[, permno := as.character(permno)]

  tr <- as.character(pdix[date_index == ix, permno])
  tr <- intersect(tr, unique(pan$permno))
  stopifnot(length(tr) == 1L)
  msg("  [%d/%d] date_index %d (%s): %d units (1 treated), fitting gsynth ...",
      k, nrow(pick), ix, ann, uniqueN(pan$permno))

  t0 <- Sys.time()
  fit <- event_study(pan, unit = "permno", time = "event_date",
                     ret = "daret", treated = tr, event_time = 0,
                     method = "gsynth", force = "unit", r = c(1, 10),
                     window = c(-30, 250), est_window = c(-280, -31),
                     returns = "simple", cumulate = "log",
                     se = "none", keep_data = FALSE)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  ref <- sc[ann_tdate == ann & permno == as.numeric(tr)]
  setorder(ref, event_date)
  stopifnot(identical(as.numeric(names(fit$paths$synthetic)),
                      as.numeric(ref$event_date)))
  yct_new <- unname(fit$paths$synthetic)
  yct_old <- ref$daret_sc
  ev <- ref$event_date %in% c(-1, 0, 1)
  car_new <- sum(log1p(unname(fit$paths$treated)[ev])) - sum(log1p(yct_new[ev]))
  car_old <- sum(log1p(ref$daret_treated[ev])) - sum(log1p(yct_old[ev]))
  val[[k]] <- data.table(
    date_index = ix, ann_tdate = ann, permno = tr,
    n_units = uniqueN(pan$permno), r_cv = fit$diagnostics$info$r,
    cor_path = cor(yct_new, yct_old),
    mae_path = mean(abs(yct_new - yct_old)),
    max_abs_diff = max(abs(yct_new - yct_old)),
    car_log_3d_feventr = car_new, car_log_3d_saved = car_old,
    car_diff = car_new - car_old, fit_seconds = round(secs, 1))
  msg("    r.cv=%d  cor=%.6f  MAE=%.2e  3d CAR: feventr %.6f vs saved %.6f (%.1fs)",
      val[[k]]$r_cv, val[[k]]$cor_path, val[[k]]$mae_path,
      car_new, car_old, secs)
}
val <- rbindlist(val)
write.csv(val, out_path("table6_gsynth_validation.csv"), row.names = FALSE)
msg("  wrote %s", out_path("table6_gsynth_validation.csv"))
msg("Done.")
