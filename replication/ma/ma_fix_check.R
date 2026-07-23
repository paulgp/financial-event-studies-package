# Do the estimand-preserving corrections kill the SC placebo drift?
# Four placebo columns, identical seeded placebo units, additive CATT,
# same block bootstrap:
#   SC baseline   — placebo_sc          (fixed-weight simplex, no fix)
#   SC + intercept— placebo_demean_sc   (est-window unit demeaning)
#   SC + ABK      — placebo_abk_sc      (prior-gross-return weighted
#                                        counterfactual; donor leg purged,
#                                        treated leg untouched)
#   Gsynth        — placebo_gsynth      (unit FE + factor averaging)
# Predictions (MEMO_longrun_bias.md): intercept collapses the drift to
# ~0; ABK flips it positive by the unit-side inflation it cannot touch.
# Writes output/ma_fix_check.{csv,png}.
#
# Run from replication/: Rscript ma/ma_fix_check.R
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

read_pl <- function(m) {
  fs <- list.files(file.path("ma", "ma_refit_out", m),
                   pattern = "^cohort_[0-9]+[.]csv$", full.names = TRUE)
  if (!length(fs)) stop("no fits in ma_refit_out/", m)
  d <- rbindlist(lapply(fs, fread), fill = TRUE)[status == "ok"]
  setorder(d, permno, date_index, event_date)
  d[, car_arith := cumsum(fifelse(is.finite(att), att, 0)),
    by = .(permno, date_index)]
  d[, base := car_arith[event_date == -2L], by = .(permno, date_index)]
  d <- d[event_date >= -1L,
         .(tag = as.numeric(matched_permno), date_index,
           ann_tdate = as.character(ann_tdate),
           event_date, car = car_arith - base)]
  merge(d, deals[, .(tag = permno, date_index, pct_cash, pct_stk)],
        by = c("tag", "date_index"), allow.cartesian = TRUE)
}
paths <- list("SC baseline"    = read_pl("placebo_sc"),
              "SC + intercept" = read_pl("placebo_demean_sc"),
              "SC + ABK"       = read_pl("placebo_abk_sc"),
              "Gsynth"         = read_pl("placebo_gsynth"))

SUBS <- function(d) list("Full sample"  = rep(TRUE, nrow(d)),
                         "Cash merger"  = d$pct_cash == 100,
                         "Stock merger" = d$pct_stk == 100)
months_all <- sort(unique(substr(
  unlist(lapply(paths, function(d) d$ann_tdate)), 1, 7)))
nm <- length(months_all)
nb <- ceiling(nm / BLOCK_MONTHS)
starts <- matrix(sample.int(nm, B_BOOT * nb, replace = TRUE), B_BOOT, nb)

pool_boot <- function(d, m) {
  d <- d[is.finite(car)]
  d[, month := substr(ann_tdate, 1, 7)]
  subs <- SUBS(d)
  out <- list()
  for (s in names(subs)) {
    ds <- d[subs[[s]]]
    horizons <- sort(unique(ds$event_date))
    agg <- ds[, .(s_car = sum(car), n = .N), by = .(month, event_date)]
    S <- matrix(0, nm, length(horizons),
                dimnames = list(months_all, horizons))
    N <- S
    S[cbind(match(agg$month, months_all),
            match(agg$event_date, horizons))] <- agg$s_car
    N[cbind(match(agg$month, months_all),
            match(agg$event_date, horizons))] <- agg$n
    reps <- matrix(NA_real_, B_BOOT, length(horizons))
    for (b in seq_len(B_BOOT)) {
      sel <- as.vector(vapply(starts[b, ], function(st)
        ((st - 1L + 0:(BLOCK_MONTHS - 1L)) %% nm) + 1L,
        integer(BLOCK_MONTHS)))
      reps[b, ] <- colSums(S[sel, , drop = FALSE]) /
        pmax(colSums(N[sel, , drop = FALSE]), 1)
    }
    out[[s]] <- data.table(
      design = m, subsample = s, event_date = horizons,
      mean_car = colSums(S) / colSums(N),
      boot_se = apply(reps, 2, sd),
      lo = apply(reps, 2, quantile, 0.025),
      hi = apply(reps, 2, quantile, 0.975),
      n = as.integer(colSums(N)))
  }
  rbindlist(out)
}
pp <- rbindlist(Map(pool_boot, paths, names(paths)))
write.csv(pp, out_path("ma_fix_check.csv"), row.names = FALSE)

hz <- c(1, 21, 63, 126, 250)
cat("\nPlacebo drift by design (additive CATT %, 95% CI):\n")
print(dcast(pp[subsample == "Full sample" & event_date %in% hz,
               .(design, event_date,
                 cell = sprintf("%6.2f [%6.2f, %6.2f]", 100 * mean_car,
                                100 * lo, 100 * hi))],
            design ~ event_date, value.var = "cell"), row.names = FALSE)
# per-day bias by era, from the raw fits
cat("\nmean daily att (bp/day), post 0..250, by era:\n")
for (m in names(paths)) {
  fs <- file.path("ma", "ma_refit_out",
                  c("SC baseline" = "placebo_sc",
                    "SC + intercept" = "placebo_demean_sc",
                    "SC + ABK" = "placebo_abk_sc",
                    "Gsynth" = "placebo_gsynth")[[m]])
  d <- rbindlist(lapply(list.files(fs, pattern = "^cohort_",
                                   full.names = TRUE), fread),
                 fill = TRUE)[status == "ok"]
  d[, year := year(ann_tdate)]
  cat(sprintf("  %-15s %6.1f overall | %6.1f pre-2001 | %6.1f post\n", m,
              1e4 * d[event_date >= 0 & is.finite(att), mean(att)],
              1e4 * d[event_date >= 0 & is.finite(att) & year <= 2000,
                      mean(att)],
              1e4 * d[event_date >= 0 & is.finite(att) & year >= 2001,
                      mean(att)]))
}

## ---- figure: full-sample placebo paths under each design --------------------
cols <- c("SC baseline" = "#008300", "SC + intercept" = "#7a49a5",
          "SC + ABK" = "#b8860b", "Gsynth" = "#2a78d6")
png(out_path("ma_fix_check.png"), width = 1700, height = 1100, res = 240)
par(mar = c(4, 4.2, 2.5, 9), mgp = c(2.4, 0.7, 0), las = 1)
pa <- pp[subsample == "Full sample"]
ylim <- 100 * range(pa$lo, pa$hi, finite = TRUE)
plot(NA, xlim = c(-1, 255), ylim = ylim, xlab = "Event day",
     ylab = "Placebo mean additive CATT from day -1 (%)",
     main = "Placebo drift under each counterfactual design (full sample)",
     cex.main = 0.95, font.main = 1, bty = "n", xaxt = "n")
axis(1, at = c(0, 50, 100, 150, 200, 250))
grid(nx = NA, ny = NULL, col = "#d8d8d8", lty = 1, lwd = 0.5)
abline(h = 0, col = "#9a9a9a", lwd = 0.7)
for (m in names(cols)) {
  x <- pa[design == m][order(event_date)]
  polygon(c(x$event_date, rev(x$event_date)), 100 * c(x$lo, rev(x$hi)),
          col = grDevices::adjustcolor(cols[[m]], 0.13), border = NA)
  lines(x$event_date, 100 * x$mean_car, col = cols[[m]], lwd = 2)
  yend <- 100 * x[event_date == 250, mean_car]
  text(255, yend, sprintf("%s %.1f%%", m, yend), pos = 4, cex = 0.72,
       col = cols[[m]], xpd = NA)
}
legend("bottomleft", legend = names(cols), col = cols, lwd = 2, bty = "n",
       cex = 0.75)
dev.off()
cat("\nwrote", out_path("ma_fix_check.csv"), "and",
    out_path("ma_fix_check.png"), "\n")
