# Monthly-frequency (3-year-horizon) M&A results: pooled additive CATTs
# for three designs -- monthly-fit SC, monthly gsynth, and the hybrid
# (daily-fit SC weights evaluated on monthly returns) -- each with its
# random placebo (mechanical bias) and runup-matched placebo (selection
# reversion). CAR(h) = sum of att over event months 0..h (month -1 is
# the onset buffer). Circular block bootstrap over announcement months,
# shared draws; corrected effect = treated - runup placebo with
# replicate-level CIs.
#
# Writes output/ma_monthly_longrun.{csv,png}. Run from replication/:
#   Rscript ma/ma_monthly_longrun.R
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
deals[, ann_month := substr(ann_tdate, 1, 7)]

read_col <- function(dir, tag_col) {
  fs <- list.files(file.path("ma", "ma_refit_out", dir),
                   pattern = "^cohort_.*[.]csv$", full.names = TRUE)
  if (!length(fs)) stop("no fits in ", dir)
  d <- rbindlist(lapply(fs, fread), fill = TRUE)[status == "ok"]
  d <- d[event_month >= 0L]
  setorder(d, permno, ann_month, event_month)
  d[, car := cumsum(fifelse(is.finite(att), att, 0)),
    by = .(permno, ann_month)]
  d <- d[, .(tag = as.numeric(get(tag_col)),
             ann_month = as.character(ann_month), event_month, car)]
  merge(d, deals[, .(tag = permno, ann_month, pct_cash, pct_stk)],
        by = c("tag", "ann_month"), allow.cartesian = TRUE)
}
cols <- list(
  "Monthly SC|Treated"        = read_col("monthly_sc", "permno"),
  "Monthly SC|Random placebo" = read_col("monthly_placebo_sc",
                                         "matched_permno"),
  "Monthly SC|Runup placebo"  = read_col("monthly_placebo_runup_sc",
                                         "matched_permno"),
  "Gsynth|Treated"            = read_col("monthly_gsynth", "permno"),
  "Gsynth|Random placebo"     = read_col("monthly_placebo_gsynth",
                                         "matched_permno"),
  "Gsynth|Runup placebo"      = read_col("monthly_placebo_runup_gsynth",
                                         "matched_permno"),
  "Hybrid SC|Treated"         = read_col("hybrid_sc", "permno"),
  "Hybrid SC|Random placebo"  = read_col("hybrid_placebo_sc",
                                         "matched_permno"),
  "Hybrid SC|Runup placebo"   = read_col("hybrid_placebo_runup_sc",
                                         "matched_permno"))

SUBS <- function(d) list("Full sample"  = rep(TRUE, nrow(d)),
                         "Cash merger"  = d$pct_cash == 100,
                         "Stock merger" = d$pct_stk == 100)
months_all <- sort(unique(unlist(lapply(cols, function(d) d$ann_month))))
nm <- length(months_all)
nb <- ceiling(nm / BLOCK_MONTHS)
cat("announcement months:", nm, "\n")
starts <- matrix(sample.int(nm, B_BOOT * nb, replace = TRUE), B_BOOT, nb)

grid_SN <- function(ds, horizons) {
  agg <- ds[, .(s = sum(car), n = .N), by = .(ann_month, event_month)]
  S <- matrix(0, nm, length(horizons), dimnames = list(months_all, horizons))
  N <- S
  S[cbind(match(agg$ann_month, months_all),
          match(agg$event_month, horizons))] <- agg$s
  N[cbind(match(agg$ann_month, months_all),
          match(agg$event_month, horizons))] <- agg$n
  list(S = S, N = N)
}
sel_rows <- function(b) as.vector(vapply(starts[b, ], function(st)
  ((st - 1L + 0:(BLOCK_MONTHS - 1L)) %% nm) + 1L, integer(BLOCK_MONTHS)))
rep_mean <- function(g, sel) colSums(g$S[sel, , drop = FALSE]) /
  pmax(colSums(g$N[sel, , drop = FALSE]), 1)

out <- list()
for (s in c("Full sample", "Cash merger", "Stock merger")) {
  horizons <- 0:36
  gs <- lapply(cols, function(d) {
    d <- d[is.finite(car)]
    grid_SN(d[SUBS(d)[[s]]], horizons)
  })
  reps <- lapply(gs, function(g) {
    m <- matrix(NA_real_, B_BOOT, length(horizons))
    for (b in seq_len(B_BOOT)) m[b, ] <- rep_mean(g, sel_rows(b))
    m
  })
  pt <- lapply(gs, function(g) colSums(g$S) / colSums(g$N))
  ser <- list()
  for (nmx in names(cols)) ser[[nmx]] <- list(pt = pt[[nmx]],
                                              rp = reps[[nmx]])
  for (dsg in c("Monthly SC", "Gsynth", "Hybrid SC")) {
    tr <- paste0(dsg, "|Treated")
    ru <- paste0(dsg, "|Runup placebo")
    ser[[paste0(dsg, "|Corrected")]] <-
      list(pt = pt[[tr]] - pt[[ru]], rp = reps[[tr]] - reps[[ru]])
  }
  out[[s]] <- rbindlist(lapply(names(ser), function(nmx) {
    parts <- strsplit(nmx, "|", fixed = TRUE)[[1]]
    data.table(design = parts[1], series = parts[2], subsample = s,
               event_month = horizons, mean_car = ser[[nmx]]$pt,
               boot_se = apply(ser[[nmx]]$rp, 2, sd),
               lo = apply(ser[[nmx]]$rp, 2, quantile, 0.025),
               hi = apply(ser[[nmx]]$rp, 2, quantile, 0.975),
               n = as.integer(colSums(gs[[if (parts[2] == "Corrected")
                 paste0(parts[1], "|Treated") else nmx]]$N)))
  }))
}
pp <- rbindlist(out)
write.csv(pp, out_path("ma_monthly_longrun.csv"), row.names = FALSE)

hz <- c(3, 12, 24, 36)
cat("\nFull sample, additive CATT % [95% CI], by design and horizon (months):\n")
print(dcast(pp[subsample == "Full sample" & event_month %in% hz,
               .(design, series, event_month,
                 cell = sprintf("%6.2f [%6.2f,%6.2f]", 100 * mean_car,
                                100 * lo, 100 * hi))],
            design + series ~ event_month, value.var = "cell"),
      row.names = FALSE)
cat("\nCorrected effects at +36m by subsample:\n")
print(dcast(pp[series == "Corrected" & event_month == 36,
               .(design, subsample,
                 cell = sprintf("%6.2f [%6.2f,%6.2f]", 100 * mean_car,
                                100 * lo, 100 * hi))],
            design ~ subsample, value.var = "cell"), row.names = FALSE)

## ---- figure: random-placebo drift (left) + corrected effects (right) --------
dcols <- c("Monthly SC" = "#008300", "Gsynth" = "#2a78d6",
           "Hybrid SC" = "#7a49a5")
png(out_path("ma_monthly_longrun.png"), width = 2600, height = 1100,
    res = 240)
par(mfrow = c(1, 2), mar = c(4, 4.2, 2.8, 2), mgp = c(2.4, 0.7, 0),
    las = 1)
pa <- pp[subsample == "Full sample"]
for (panel in c("Random placebo", "Corrected")) {
  px <- pa[series == panel]
  ylim <- 100 * range(px$lo, px$hi, finite = TRUE)
  plot(NA, xlim = c(0, 38), ylim = ylim, xlab = "Event month",
       ylab = "Mean additive CATT from month 0 (%)",
       main = if (panel == "Random placebo")
         "Random-placebo drift (mechanical bias)" else
           "Corrected effect (treated - runup placebo)",
       cex.main = 0.95, font.main = 1, bty = "n")
  grid(nx = NA, ny = NULL, col = "#d8d8d8", lty = 1, lwd = 0.5)
  abline(h = 0, col = "#9a9a9a", lwd = 0.7)
  for (dsg in names(dcols)) {
    x <- px[design == dsg][order(event_month)]
    polygon(c(x$event_month, rev(x$event_month)),
            100 * c(x$lo, rev(x$hi)),
            col = grDevices::adjustcolor(dcols[[dsg]], 0.13), border = NA)
    lines(x$event_month, 100 * x$mean_car, col = dcols[[dsg]], lwd = 2)
  }
  if (panel == "Random placebo")
    legend("bottomleft", legend = names(dcols), col = dcols, lwd = 2,
           bty = "n", cex = 0.8)
}
dev.off()
cat("\nwrote", out_path("ma_monthly_longrun.csv"), "and",
    out_path("ma_monthly_longrun.png"), "\n")
