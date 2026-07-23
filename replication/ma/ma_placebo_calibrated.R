# Placebo-calibrated long-run CATT — prototype of correction #1 in
# MEMO_longrun_bias.md. The estimator is (pooled treated path) minus
# (pooled date-matched placebo path), with the difference taken INSIDE
# each circular-block-bootstrap replicate (shared month blocks), so the
# 95% band is on the calibrated effect itself. No estimand change, no
# model of the price noise: the placebo carries the realized design bias
# out of the estimate.
#
# Methods: SC (feventr refit treated vs placebo_sc) and Gsynth
# (published per-deal fits vs placebo_gsynth; same seeded placebo units).
# Additive CATT only. Writes output/ma_placebo_calibrated.{csv,png}.
#
# Run from replication/: Rscript ma/ma_placebo_calibrated.R
source("config.R")
suppressMessages({library(haven); library(data.table)})
ma_work <- function(...) dind("M&A", "data", "work", ...)
ma_out  <- function(...) dind("M&A", "output", ...)

B_BOOT <- 1000L
BLOCK_MONTHS <- 18L
set.seed(20260721)

## ---- deal sample -------------------------------------------------------------
di <- as.data.table(read_dta(
  ma_work("sdc_ma_details_sl_m_cleaned_dateindex_2023.dta")))
saved <- as.data.table(read_dta(ma_out("sl_m_deals_car_1_250_gsynth.dta"),
                                col_select = c("group")))
det <- as.data.table(read_dta(
  ma_work("sdc_ma_details_sl_m_cleaned_2023.dta"),
  col_select = c("master_deal_no", "permno", "pct_cash", "pct_stk")))
deals <- merge(di[, .(master_deal_no, permno, group, date_index,
                      ann_tdate = as.character(as.Date(date)))],
               saved, by = "group")
deals <- merge(deals, det, by = c("master_deal_no", "permno"))

## ---- per-deal additive CATT paths -------------------------------------------
read_refit <- function(m, tag_col = "permno") {
  fs <- list.files(file.path("ma", "ma_refit_out", m),
                   pattern = "^cohort_[0-9]+[.]csv$", full.names = TRUE)
  if (!length(fs)) stop("no fits in ma_refit_out/", m)
  d <- rbindlist(lapply(fs, fread), fill = TRUE)[status == "ok"]
  setorder(d, permno, date_index, event_date)
  d[, car_arith := cumsum(fifelse(is.finite(att), att, 0)),
    by = .(permno, date_index)]
  d[, base := car_arith[event_date == -2L], by = .(permno, date_index)]
  d <- d[event_date >= -1L,
         .(tag = as.numeric(get(tag_col)), date_index,
           ann_tdate = as.character(ann_tdate),
           event_date, car = car_arith - base)]
  merge(d, deals[, .(tag = permno, date_index, pct_cash, pct_stk)],
        by = c("tag", "date_index"), allow.cartesian = TRUE)
}
read_gsynth_treated <- function() {
  sc <- as.data.table(read_dta(
    ma_work("sdc_ma_sl_m_gsynth", "sc_ma_1_7052.dta"),
    col_select = c("permno", "event_date", "daret_treated", "daret_sc",
                   "ann_tdate")))
  sc <- sc[event_date >= -1L & event_date <= 250L]
  sc[, ann_tdate := as.character(as.Date(ann_tdate))]
  setorder(sc, permno, ann_tdate, event_date)
  sc[, car := cumsum(fcoalesce(daret_treated, 0)) -
       cumsum(fcoalesce(daret_sc, 0)), by = .(permno, ann_tdate)]
  sc <- sc[, n_days := .N, by = .(permno, ann_tdate)][n_days == 252L]
  merge(sc[, .(permno, ann_tdate, event_date, car)],
        deals[, .(permno, ann_tdate, date_index, pct_cash, pct_stk)],
        by = c("permno", "ann_tdate"), allow.cartesian = TRUE)
}

pairs <- list(
  SC = list(tr = read_refit("sc"), pl = read_refit("placebo_sc",
                                                   "matched_permno")),
  Gsynth = list(tr = read_gsynth_treated(),
                pl = read_refit("placebo_gsynth", "matched_permno")))

SUBS <- function(d) list("Full sample"  = rep(TRUE, nrow(d)),
                         "Cash merger"  = d$pct_cash == 100,
                         "Stock merger" = d$pct_stk == 100)
months_all <- sort(unique(substr(unlist(lapply(pairs, function(p)
  c(p$tr$ann_tdate, p$pl$ann_tdate))), 1, 7)))
nm <- length(months_all)
nb <- ceiling(nm / BLOCK_MONTHS)
starts <- matrix(sample.int(nm, B_BOOT * nb, replace = TRUE), B_BOOT, nb)

grid_SN <- function(ds, horizons) {
  agg <- ds[, .(s_car = sum(car), n = .N), by = .(month, event_date)]
  S <- matrix(0, nm, length(horizons), dimnames = list(months_all, horizons))
  N <- S
  S[cbind(match(agg$month, months_all),
          match(agg$event_date, horizons))] <- agg$s_car
  N[cbind(match(agg$month, months_all),
          match(agg$event_date, horizons))] <- agg$n
  list(S = S, N = N)
}

calibrate <- function(p, m) {
  tr <- p$tr[is.finite(car)][, month := substr(ann_tdate, 1, 7)]
  pl <- p$pl[is.finite(car)][, month := substr(ann_tdate, 1, 7)]
  out <- list()
  for (s in names(SUBS(tr))) {
    ts <- tr[SUBS(tr)[[s]]]
    ps <- pl[SUBS(pl)[[s]]]
    horizons <- sort(intersect(unique(ts$event_date),
                               unique(ps$event_date)))
    gt <- grid_SN(ts, horizons)
    gp <- grid_SN(ps, horizons)
    reps <- matrix(NA_real_, B_BOOT, length(horizons))
    for (b in seq_len(B_BOOT)) {
      sel <- as.vector(vapply(starts[b, ], function(st)
        ((st - 1L + 0:(BLOCK_MONTHS - 1L)) %% nm) + 1L,
        integer(BLOCK_MONTHS)))
      reps[b, ] <- colSums(gt$S[sel, , drop = FALSE]) /
        pmax(colSums(gt$N[sel, , drop = FALSE]), 1) -
        colSums(gp$S[sel, , drop = FALSE]) /
        pmax(colSums(gp$N[sel, , drop = FALSE]), 1)
    }
    out[[s]] <- data.table(
      method = m, subsample = s, event_date = horizons,
      calibrated = colSums(gt$S) / colSums(gt$N) -
        colSums(gp$S) / colSums(gp$N),
      boot_se = apply(reps, 2, sd),
      lo = apply(reps, 2, quantile, 0.025),
      hi = apply(reps, 2, quantile, 0.975),
      n = as.integer(colSums(gt$N)))
  }
  rbindlist(out)
}
pp <- rbindlist(Map(calibrate, pairs, names(pairs)))
write.csv(pp, out_path("ma_placebo_calibrated.csv"), row.names = FALSE)

hz <- c(1, 21, 63, 126, 250)
cat("\nPlacebo-calibrated CATT (%), 95% block-bootstrap CI on the difference:\n")
print(dcast(pp[event_date %in% hz,
               .(method, subsample, event_date,
                 cell = sprintf("%6.2f [%6.2f, %6.2f]", 100 * calibrated,
                                100 * lo, 100 * hi))],
            method + subsample ~ event_date, value.var = "cell"),
      row.names = FALSE)

## ---- figure -------------------------------------------------------------------
cols <- c(Gsynth = "#2a78d6", SC = "#008300")
png(out_path("ma_placebo_calibrated.png"), width = 2600, height = 1000,
    res = 240)
par(mfrow = c(1, 3), mar = c(4, 4.2, 2.5, 2), mgp = c(2.4, 0.7, 0), las = 1)
ylim <- 100 * range(pp$lo, pp$hi, finite = TRUE)
for (s in c("Full sample", "Cash merger", "Stock merger")) {
  plot(NA, xlim = c(-1, 255), ylim = ylim, xlab = "Event day",
       ylab = if (s == "Full sample")
         "Placebo-calibrated CATT from day -1 (%)" else "",
       main = s, cex.main = 0.95, font.main = 1, bty = "n", xaxt = "n")
  axis(1, at = c(0, 50, 100, 150, 200, 250))
  grid(nx = NA, ny = NULL, col = "#d8d8d8", lty = 1, lwd = 0.5)
  abline(h = 0, col = "#9a9a9a", lwd = 0.7)
  for (m in names(cols)) {
    x <- pp[method == m & subsample == s][order(event_date)]
    polygon(c(x$event_date, rev(x$event_date)),
            100 * c(x$lo, rev(x$hi)),
            col = grDevices::adjustcolor(cols[[m]], 0.15), border = NA)
    lines(x$event_date, 100 * x$calibrated, col = cols[[m]], lwd = 2)
  }
  if (s == "Full sample")
    legend("bottomleft",
           legend = c("Gsynth  (treated - placebo)", "SC  (treated - placebo)",
                      "bands: 95% block bootstrap of the difference"),
           col = c(cols, NA), lwd = c(2, 2, NA), bty = "n", cex = 0.75)
}
dev.off()
cat("\nwrote", out_path("ma_placebo_calibrated.csv"), "and",
    out_path("ma_placebo_calibrated.png"), "\n")
