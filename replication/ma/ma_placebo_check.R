# M&A placebo (randomization) check — gold-standard validation of the
# long-run inference. ma_refit_full.R with MA_REFIT_PLACEBO=1 fits one
# date-matched placebo donor per deal (real acquirers excluded from every
# pool), so the pooled placebo path carries the true calendar dependence
# and, under no effect, estimates zero plus any design bias (e.g. from the
# complete-531-day-coverage requirement). Two diagnostics:
#   1. bias: does the pooled placebo path hug zero?
#   2. inference: does the treated path lie outside the placebo path's 95%
#      block-bootstrap band (same 18-month announcement blocks as
#      ma_refit_longrun.R)?
#
# Writes output/ma_placebo_check.csv and output/ma_placebo_check.png.
#
# Run from replication/ (after the placebo fits): Rscript ma/ma_placebo_check.R [method]
source("config.R")
suppressMessages({library(haven); library(data.table)})
ma_work <- function(...) dind("M&A", "data", "work", ...)
ma_out  <- function(...) dind("M&A", "output", ...)

method <- commandArgs(trailingOnly = TRUE)
method <- if (length(method)) method[1] else "sc"
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

read_paths <- function(dir, tag_col) {
  fs <- list.files(file.path("ma", "ma_refit_out", dir),
                   pattern = "^cohort_[0-9]+[.]csv$", full.names = TRUE)
  if (!length(fs)) stop("no fits in ma_refit_out/", dir)
  d <- rbindlist(lapply(fs, fread), fill = TRUE)[status == "ok"]
  d[, base := car_log[event_date == -2L], by = .(permno, date_index)]
  # fread parses ann_tdate as IDate; keep it character so the month grid
  # (built via unlist, which strips the class) stays in "YYYY-MM-DD" form
  d <- d[event_date >= -1L,
         .(tag = as.numeric(get(tag_col)), date_index,
           ann_tdate = as.character(ann_tdate),
           event_date, car = car_log - base)]
  merge(d, deals[, .(tag = permno, date_index, pct_cash, pct_stk)],
        by = c("tag", "date_index"), allow.cartesian = TRUE)
}
paths <- list(Treated = read_paths(method, "permno"),
              Placebo = read_paths(paste0("placebo_", method),
                                   "matched_permno"))

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
      role = m, subsample = s, event_date = horizons,
      mean_car = colSums(S) / colSums(N),
      boot_se = apply(reps, 2, sd),
      lo = apply(reps, 2, quantile, 0.025),
      hi = apply(reps, 2, quantile, 0.975),
      n = as.integer(colSums(N)))
  }
  rbindlist(out)
}
pp <- rbindlist(Map(pool_boot, paths, names(paths)))
write.csv(pp, out_path("ma_placebo_check.csv"), row.names = FALSE)

hz <- c(1, 21, 63, 126, 250)
wide <- merge(pp[role == "Treated",
                 .(subsample, event_date, treated = mean_car)],
              pp[role == "Placebo",
                 .(subsample, event_date, placebo = mean_car,
                   plo = lo, phi = hi, pse = boot_se)],
              by = c("subsample", "event_date"))
wide[, outside := treated < plo | treated > phi]
cat(sprintf("\nPlacebo check, method = %s (%%, at selected horizons):\n",
            method))
print(dcast(wide[event_date %in% hz,
                 .(subsample, event_date,
                   cell = sprintf("t %6.2f | p %5.2f (%4.2f) %s",
                                  100 * treated, 100 * placebo, 100 * pse,
                                  fifelse(outside, "OUT", "in")))],
            subsample ~ event_date, value.var = "cell"), row.names = FALSE)
cat(sprintf("\nplacebo bias, |mean| at +250: %s\n",
            paste(sprintf("%s %.2f%%",
                          wide[event_date == 250, subsample],
                          100 * abs(wide[event_date == 250, placebo])),
                  collapse = " | ")))
cat(sprintf("treated outside placebo 95%% band at horizons >= +21: %d/%d\n",
            wide[event_date >= 21, sum(outside)],
            wide[event_date >= 21, .N]))

## ---- figure -------------------------------------------------------------------
png(out_path("ma_placebo_check.png"), width = 2600, height = 1000, res = 240)
par(mfrow = c(1, 3), mar = c(4, 4.2, 2.5, 2), mgp = c(2.4, 0.7, 0), las = 1)
ylim <- 100 * range(pp$lo, pp$hi, pp$mean_car, finite = TRUE)
for (s in c("Full sample", "Cash merger", "Stock merger")) {
  plot(NA, xlim = c(-1, 255), ylim = ylim, xlab = "Event day",
       ylab = if (s == "Full sample")
         "Mean log CAR from day -1 (%)" else "",
       main = s, cex.main = 0.95, font.main = 1, bty = "n", xaxt = "n")
  axis(1, at = c(0, 50, 100, 150, 200, 250))
  grid(nx = NA, ny = NULL, col = "#d8d8d8", lty = 1, lwd = 0.5)
  abline(h = 0, col = "#9a9a9a", lwd = 0.7)
  pl <- pp[role == "Placebo" & subsample == s][order(event_date)]
  polygon(c(pl$event_date, rev(pl$event_date)), 100 * c(pl$lo, rev(pl$hi)),
          col = "#d9d9d9", border = NA)
  lines(pl$event_date, 100 * pl$mean_car, col = "#6b6b6b", lwd = 2)
  tr <- pp[role == "Treated" & subsample == s][order(event_date)]
  lines(tr$event_date, 100 * tr$mean_car, col = "#2a78d6", lwd = 2)
  if (s == "Full sample")
    legend("bottomleft",
           legend = c(sprintf("treated (%s refit)", method),
                      "placebo mean", "placebo 95% block-bootstrap band"),
           col = c("#2a78d6", "#6b6b6b", "#d9d9d9"),
           lwd = c(2, 2, 8), bty = "n", cex = 0.75)
}
dev.off()
cat("\nwrote", out_path("ma_placebo_check.csv"), "and",
    out_path("ma_placebo_check.png"), "\n")
