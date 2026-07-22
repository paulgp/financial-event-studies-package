# M&A long-run effects — pooled [-1, h] log CAR paths out to +250 trading
# days after announcement, by deal type, for three counterfactuals: the
# published per-deal gsynth paths (sc_ma_1_7052.dta, no refit) and the
# feventr sc/apm refits (ma_refit_full.R path output). Deal sample and
# subsamples exactly as Table 6 (published 14,847 deals; 100% cash /
# 100% stock).
#
# Inference: deal CARs are cross-sectionally dependent — deals cluster on
# announcement dates and 250-day windows overlap across deals announced
# within a year of each other — so sd/sqrt(N) understates. The default
# bands here are a circular block bootstrap over ANNOUNCEMENT TIME:
# announcement months are resampled in blocks of 18 consecutive months
# (longer than the -30..+250 window), each block carrying all its deals,
# so within-date clustering and window overlap travel intact. Naive
# cross-deal SEs are reported alongside for the design-effect comparison.
#
# Writes output/ma_longrun_paths.csv (method x subsample x horizon: mean,
# naive se, block-bootstrap se and 95% band, n) and the figure
# output/ma_longrun.png, plus a horizons table on stdout.
#
# Run from replication/: Rscript ma/ma_refit_longrun.R
source("config.R")
suppressMessages({library(haven); library(data.table)})
ma_work <- function(...) dind("M&A", "data", "work", ...)
ma_out  <- function(...) dind("M&A", "output", ...)

B_BOOT <- 1000L
BLOCK_MONTHS <- 18L
set.seed(20260721)

## ---- deal sample (as ma_refit_compare.R / table6.R) -------------------------
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
cat("published deal sample:", nrow(deals), "\n")

## ---- per-deal rebased CAR paths ---------------------------------------------
read_refit_paths <- function(m) {
  fs <- list.files(file.path("ma", "ma_refit_out", m),
                   pattern = "^cohort_[0-9]+[.]csv$", full.names = TRUE)
  if (!length(fs)) return(NULL)
  d <- rbindlist(lapply(fs, fread), fill = TRUE)[status == "ok"]
  d[, base := car_log[event_date == -2L], by = .(permno, date_index)]
  # fread parses ann_tdate as IDate; keep it character so the month grid
  # (built via unlist, which strips the class) stays in "YYYY-MM-DD" form
  d <- d[event_date >= -1L, .(permno = as.numeric(permno), date_index,
                              ann_tdate = as.character(ann_tdate),
                              event_date, car = car_log - base)]
  # one acquirer path can serve several same-day deals (multi-target
  # announcements): expand 1:m to deal level, as table6.R Part A does
  merge(d, deals[, .(permno, date_index, pct_cash, pct_stk)],
        by = c("permno", "date_index"), allow.cartesian = TRUE)
}

read_gsynth_paths <- function() {
  sc <- as.data.table(read_dta(
    ma_work("sdc_ma_sl_m_gsynth", "sc_ma_1_7052.dta"),
    col_select = c("permno", "event_date", "daret_treated", "daret_sc",
                   "ann_tdate")))
  sc <- sc[event_date >= -1L & event_date <= 250L]
  sc[, ann_tdate := as.character(as.Date(ann_tdate))]
  setorder(sc, permno, ann_tdate, event_date)
  # zero-fill missing returns (delistings late in the window), mirroring
  # the panel construction and table6.R Part B; returns <= -100% make
  # log1p non-finite -- Stata's sum() skips those days as missing, so
  # treat them as zero contributions the same way
  l1 <- function(x) {
    z <- log1p(fcoalesce(x, 0))
    fifelse(is.finite(z), z, 0)
  }
  sc[, car := cumsum(l1(daret_treated)) - cumsum(l1(daret_sc)),
     by = .(permno, ann_tdate)]
  sc <- sc[, n_days := .N, by = .(permno, ann_tdate)][n_days == 252L]
  merge(sc[, .(permno, ann_tdate, event_date, car)],
        deals[, .(permno, ann_tdate, date_index, pct_cash, pct_stk)],
        by = c("permno", "ann_tdate"), allow.cartesian = TRUE)
}

paths <- list(Gsynth = read_gsynth_paths(),
              SC = read_refit_paths("sc"),
              APM = read_refit_paths("apm"))
paths <- Filter(Negate(is.null), paths)

SUBS <- function(d) list("Full sample"  = rep(TRUE, nrow(d)),
                         "Cash merger"  = d$pct_cash == 100,
                         "Stock merger" = d$pct_stk == 100)

## ---- pooled means + circular block bootstrap over announcement months -------
months_all <- sort(unique(substr(
  unlist(lapply(paths, function(d) d$ann_tdate)), 1, 7)))
nm <- length(months_all)
nb <- ceiling(nm / BLOCK_MONTHS)
cat(sprintf("bootstrap: %d announcement months, %d blocks of %d, B = %d\n",
            nm, nb, BLOCK_MONTHS, B_BOOT))
# one shared set of block draws so bands are comparable across methods
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
    # full month x horizon grids (empty months as zeros) so blocks respect
    # real calendar spacing
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
      method = m, subsample = s, event_date = horizons,
      mean_car = colSums(S) / colSums(N),
      se = ds[, .(v = sd(car) / sqrt(.N)), by = event_date][
        order(event_date), v],
      boot_se = apply(reps, 2, sd),
      lo = apply(reps, 2, quantile, 0.025),
      hi = apply(reps, 2, quantile, 0.975),
      n = as.integer(colSums(N)))
  }
  rbindlist(out)
}
pp <- rbindlist(Map(pool_boot, paths, names(paths)))
write.csv(pp, out_path("ma_longrun_paths.csv"), row.names = FALSE)

hz <- c(1, 21, 63, 126, 250)
cat("\nMean [-1,+h] log CAR (%), naive cross-deal SE / block-bootstrap SE:\n")
tabh <- dcast(pp[event_date %in% hz,
                 .(method, subsample, event_date,
                   cell = sprintf("%6.2f (%.2f / %.2f)", 100 * mean_car,
                                  100 * se, 100 * boot_se))],
              method + subsample ~ event_date, value.var = "cell")
print(tabh, row.names = FALSE)
infl <- pp[event_date >= 21, .(design_effect = median(boot_se / se)),
           by = .(method, subsample)]
cat("\nmedian SE inflation (block bootstrap / naive), horizons >= +21:\n")
print(dcast(infl, method ~ subsample, value.var = "design_effect"),
      digits = 2, row.names = FALSE)

## ---- figure: small multiples, block-bootstrap 95% bands ---------------------
cols <- c(Gsynth = "#2a78d6", SC = "#008300", APM = "#e87ba4")
png(out_path("ma_longrun.png"), width = 2600, height = 1000, res = 240)
par(mfrow = c(1, 3), mar = c(4, 4.2, 2.5, 5.5), mgp = c(2.4, 0.7, 0),
    las = 1)
ylim <- 100 * range(pp$lo, pp$hi, finite = TRUE)
for (s in c("Full sample", "Cash merger", "Stock merger")) {
  plot(NA, xlim = c(-1, 285), ylim = ylim, xlab = "Event day",
       ylab = if (s == "Full sample")
         "Mean log CAR from day -1 (%)" else "",
       main = sprintf("%s (n=%s)", s,
                      format(pp[subsample == s & event_date == 250 &
                                  method == names(paths)[1], n],
                             big.mark = ",")),
       cex.main = 0.95, font.main = 1, bty = "n", xaxt = "n")
  axis(1, at = c(0, 50, 100, 150, 200, 250))
  grid(nx = NA, ny = NULL, col = "#d8d8d8", lty = 1, lwd = 0.5)
  abline(h = 0, col = "#9a9a9a", lwd = 0.7)
  for (m in names(paths)) {
    x <- pp[method == m & subsample == s][order(event_date)]
    polygon(c(x$event_date, rev(x$event_date)),
            100 * c(x$lo, rev(x$hi)),
            col = grDevices::adjustcolor(cols[[m]], 0.15), border = NA)
    lines(x$event_date, 100 * x$mean_car, col = cols[[m]], lwd = 2)
    yend <- 100 * x[event_date == 250, mean_car]
    points(250, yend, pch = 16, col = cols[[m]], cex = 0.8)
    text(255, yend, sprintf("%s %.1f%%", m, yend), pos = 4, cex = 0.72,
         col = "#2b2b2b", xpd = NA)
  }
  if (s == "Full sample")
    legend("bottomleft",
           legend = c(names(cols[names(paths)]),
                      "bands: 95% circular block bootstrap",
                      "(18-month announcement blocks)"),
           col = c(cols[names(paths)], NA, NA),
           lwd = c(rep(2, length(paths)), NA, NA), bty = "n", cex = 0.75)
}
dev.off()
cat("\nwrote", out_path("ma_longrun_paths.csv"), "and",
    out_path("ma_longrun.png"), "\n")
