# Runup-selection (Ashenfelter dip) check for the gsynth long-run effect.
# Three pooled additive-CATT paths, same block bootstrap, shared draws:
#   Treated        — the paper's saved per-deal gsynth fits
#   Random placebo — placebo_gsynth (no selection: measures mechanical bias)
#   Runup placebo  — placebo_runup_gsynth (each placebo matched to its
#                    acquirer on estimation-window cumulative return, so it
#                    inherits the selection on pre-announcement runup)
# Decomposition, with replicate-level CIs:
#   reversion  = runup placebo - random placebo  (what the gsynth unit FE
#                wrongly carries forward if runups are transient)
#   corrected  = treated - runup placebo         (long-run effect net of
#                both mechanical bias and runup reversion)
# Writes output/ma_runup_check.{csv,png}.
#
# Run from replication/: Rscript ma/ma_runup_check.R
source("config.R")
suppressMessages({library(haven); library(data.table)})
ma_work <- function(...) dind("M&A", "data", "work", ...)
ma_out  <- function(...) dind("M&A", "output", ...)

B_BOOT <- 1000L
BLOCK_MONTHS <- 18L
set.seed(20260721)

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

read_refit <- function(m, tag_col = "matched_permno") {
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

paths <- list("Treated"        = read_gsynth_treated(),
              "Random placebo" = read_refit("placebo_gsynth"),
              "Runup placebo"  = read_refit("placebo_runup_gsynth"))

SUBS <- function(d) list("Full sample"  = rep(TRUE, nrow(d)),
                         "Cash merger"  = d$pct_cash == 100,
                         "Stock merger" = d$pct_stk == 100)
months_all <- sort(unique(substr(
  unlist(lapply(paths, function(d) d$ann_tdate)), 1, 7)))
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
sel_rows <- function(b) as.vector(vapply(starts[b, ], function(st)
  ((st - 1L + 0:(BLOCK_MONTHS - 1L)) %% nm) + 1L, integer(BLOCK_MONTHS)))
rep_mean <- function(g, sel) colSums(g$S[sel, , drop = FALSE]) /
  pmax(colSums(g$N[sel, , drop = FALSE]), 1)

out <- list()
for (s in c("Full sample", "Cash merger", "Stock merger")) {
  ds <- lapply(paths, function(d) {
    d <- d[is.finite(car)][, month := substr(ann_tdate, 1, 7)]
    d[SUBS(d)[[s]]]
  })
  horizons <- sort(Reduce(intersect, lapply(ds, function(x)
    unique(x$event_date))))
  gs <- lapply(ds, grid_SN, horizons = horizons)
  reps <- lapply(gs, function(g) {
    m <- matrix(NA_real_, B_BOOT, length(horizons))
    for (b in seq_len(B_BOOT)) m[b, ] <- rep_mean(g, sel_rows(b))
    m
  })
  pt <- lapply(gs, function(g) colSums(g$S) / colSums(g$N))
  series <- list(
    "Treated"        = list(pt = pt[[1]], rp = reps[[1]]),
    "Random placebo" = list(pt = pt[[2]], rp = reps[[2]]),
    "Runup placebo"  = list(pt = pt[[3]], rp = reps[[3]]),
    "Reversion (runup - random)" =
      list(pt = pt[[3]] - pt[[2]], rp = reps[[3]] - reps[[2]]),
    "Corrected (treated - runup)" =
      list(pt = pt[[1]] - pt[[3]], rp = reps[[1]] - reps[[3]]))
  out[[s]] <- rbindlist(lapply(names(series), function(nmx)
    data.table(series = nmx, subsample = s, event_date = horizons,
               mean_car = series[[nmx]]$pt,
               boot_se = apply(series[[nmx]]$rp, 2, sd),
               lo = apply(series[[nmx]]$rp, 2, quantile, 0.025),
               hi = apply(series[[nmx]]$rp, 2, quantile, 0.975))))
}
pp <- rbindlist(out)
write.csv(pp, out_path("ma_runup_check.csv"), row.names = FALSE)

hz <- c(1, 21, 63, 126, 250)
cat("\nGsynth runup check (additive CATT %, 95% block-bootstrap CI):\n")
print(dcast(pp[event_date %in% hz,
               .(series, subsample, event_date,
                 cell = sprintf("%6.2f [%6.2f, %6.2f]", 100 * mean_car,
                                100 * lo, 100 * hi))],
            series + subsample ~ event_date, value.var = "cell"),
      row.names = FALSE)

## ---- figure -------------------------------------------------------------------
png(out_path("ma_runup_check.png"), width = 2600, height = 1000, res = 240)
par(mfrow = c(1, 3), mar = c(4, 4.2, 2.5, 2), mgp = c(2.4, 0.7, 0), las = 1)
show <- c("Treated" = "#2a78d6", "Random placebo" = "#6b6b6b",
          "Runup placebo" = "#e87316")
band <- c("Random placebo", "Runup placebo")
ylim <- 100 * range(pp[series %in% names(show), .(lo, hi, mean_car)],
                    finite = TRUE)
for (s in c("Full sample", "Cash merger", "Stock merger")) {
  plot(NA, xlim = c(-1, 255), ylim = ylim, xlab = "Event day",
       ylab = if (s == "Full sample")
         "Mean additive CATT from day -1 (%)" else "",
       main = s, cex.main = 0.95, font.main = 1, bty = "n", xaxt = "n")
  axis(1, at = c(0, 50, 100, 150, 200, 250))
  grid(nx = NA, ny = NULL, col = "#d8d8d8", lty = 1, lwd = 0.5)
  abline(h = 0, col = "#9a9a9a", lwd = 0.7)
  for (nmx in band) {
    x <- pp[series == nmx & subsample == s][order(event_date)]
    polygon(c(x$event_date, rev(x$event_date)), 100 * c(x$lo, rev(x$hi)),
            col = grDevices::adjustcolor(show[[nmx]], 0.15), border = NA)
    lines(x$event_date, 100 * x$mean_car, col = show[[nmx]], lwd = 1.8,
          lty = 2)
  }
  x <- pp[series == "Treated" & subsample == s][order(event_date)]
  lines(x$event_date, 100 * x$mean_car, col = show[["Treated"]], lwd = 2)
  if (s == "Full sample")
    legend("bottomleft",
           legend = c("treated (published gsynth fits)",
                      "random placebo (mechanical bias)",
                      "runup-matched placebo (+ selection reversion)",
                      "bands: 95% circular block bootstrap"),
           col = c(show, NA), lwd = c(2, 1.8, 1.8, NA),
           lty = c(1, 2, 2, NA), bty = "n", cex = 0.72)
}
dev.off()
cat("\nwrote", out_path("ma_runup_check.csv"), "and",
    out_path("ma_runup_check.png"), "\n")
