# Gsynth version of the M&A placebo check. Placebo paths come from
# ma_refit_full.R gsynth with MA_REFIT_PLACEBO=1 (feventr gsynth engine,
# one date-matched non-acquirer per deal; the per-cohort seed matches the
# SC placebo run, so the placebo UNITS are identical across methods).
# Treated paths are the paper's saved per-deal gsynth fits
# (sc_ma_1_7052.dta), so the treated line matches Figure 8 / the repo's
# long-run figures; note the treated fits are the paper's gsynth spec
# while the placebos use feventr's wrapper of the same CRAN gsynth.
#
# Additive CATT headline, log in the CSV. Writes
# output/ma_placebo_gsynth.{csv,png} and, if the SC placebo CSV exists,
# a combined two-band figure output/ma_longrun_placebo_bands.png where
# each method is judged against its own placebo band.
#
# Run from replication/: Rscript ma/ma_placebo_gsynth.R
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

read_placebo <- function(dir) {
  fs <- list.files(file.path("ma", "ma_refit_out", dir),
                   pattern = "^cohort_[0-9]+[.]csv$", full.names = TRUE)
  if (!length(fs)) stop("no fits in ma_refit_out/", dir)
  d <- rbindlist(lapply(fs, fread), fill = TRUE)[status == "ok"]
  setorder(d, permno, date_index, event_date)
  d[, car_arith := cumsum(fifelse(is.finite(att), att, 0)),
    by = .(permno, date_index)]
  d[, `:=`(base_l = car_log[event_date == -2L],
           base_a = car_arith[event_date == -2L]),
    by = .(permno, date_index)]
  d <- d[event_date >= -1L,
         .(tag = as.numeric(matched_permno), date_index,
           ann_tdate = as.character(ann_tdate),
           event_date, log = car_log - base_l, arith = car_arith - base_a)]
  d <- melt(d, measure.vars = c("arith", "log"), variable.name = "metric",
            value.name = "car")
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
  l1 <- function(x) {
    z <- log1p(fcoalesce(x, 0))
    fifelse(is.finite(z), z, 0)
  }
  sc[, `:=`(log = cumsum(l1(daret_treated)) - cumsum(l1(daret_sc)),
            arith = cumsum(fcoalesce(daret_treated, 0)) -
              cumsum(fcoalesce(daret_sc, 0))),
     by = .(permno, ann_tdate)]
  sc <- sc[, n_days := .N, by = .(permno, ann_tdate)][n_days == 252L]
  sc <- melt(sc[, .(permno, ann_tdate, event_date, arith, log)],
             measure.vars = c("arith", "log"), variable.name = "metric",
             value.name = "car")
  merge(sc, deals[, .(permno, ann_tdate, date_index, pct_cash, pct_stk)],
        by = c("permno", "ann_tdate"), allow.cartesian = TRUE)
}

paths <- list(Treated = read_gsynth_treated(),
              Placebo = read_placebo("placebo_gsynth"))

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
pp <- rbindlist(lapply(names(paths), function(m)
  rbindlist(lapply(c("arith", "log"), function(g)
    pool_boot(paths[[m]][metric == g], m)[, metric := g]))))
write.csv(pp, out_path("ma_placebo_gsynth.csv"), row.names = FALSE)

pa <- pp[metric == "arith"]
hz <- c(1, 21, 63, 126, 250)
wide <- merge(pa[role == "Treated",
                 .(subsample, event_date, treated = mean_car)],
              pa[role == "Placebo",
                 .(subsample, event_date, placebo = mean_car,
                   plo = lo, phi = hi, pse = boot_se)],
              by = c("subsample", "event_date"))
wide[, outside := treated < plo | treated > phi]
cat("\nGsynth placebo check (additive CATT %, selected horizons):\n")
print(dcast(wide[event_date %in% hz,
                 .(subsample, event_date,
                   cell = sprintf("t %6.2f | p %5.2f (%4.2f) %s",
                                  100 * treated, 100 * placebo, 100 * pse,
                                  fifelse(outside, "OUT", "in")))],
            subsample ~ event_date, value.var = "cell"), row.names = FALSE)
cat(sprintf("\nplacebo drift at +250: %s\n",
            paste(sprintf("%s %.1f%%",
                          wide[event_date == 250, subsample],
                          100 * wide[event_date == 250, placebo]),
                  collapse = " | ")))
cat(sprintf("treated outside placebo 95%% band at horizons >= +21: %d/%d\n",
            wide[event_date >= 21, sum(outside)],
            wide[event_date >= 21, .N]))

## ---- figure: gsynth treated vs gsynth placebo band --------------------------
png(out_path("ma_placebo_gsynth.png"), width = 2600, height = 1000,
    res = 240)
par(mfrow = c(1, 3), mar = c(4, 4.2, 2.5, 2), mgp = c(2.4, 0.7, 0),
    las = 1)
ylim <- 100 * range(pa$lo, pa$hi, pa$mean_car, finite = TRUE)
for (s in c("Full sample", "Cash merger", "Stock merger")) {
  plot(NA, xlim = c(-1, 255), ylim = ylim, xlab = "Event day",
       ylab = if (s == "Full sample")
         "Mean additive CATT from day -1 (%)" else "",
       main = s, cex.main = 0.95, font.main = 1, bty = "n", xaxt = "n")
  axis(1, at = c(0, 50, 100, 150, 200, 250))
  grid(nx = NA, ny = NULL, col = "#d8d8d8", lty = 1, lwd = 0.5)
  abline(h = 0, col = "#9a9a9a", lwd = 0.7)
  pl <- pa[role == "Placebo" & subsample == s][order(event_date)]
  polygon(c(pl$event_date, rev(pl$event_date)), 100 * c(pl$lo, rev(pl$hi)),
          col = "#d9d9d9", border = NA)
  lines(pl$event_date, 100 * pl$mean_car, col = "#6b6b6b", lwd = 2,
        lty = 2)
  tr <- pa[role == "Treated" & subsample == s][order(event_date)]
  lines(tr$event_date, 100 * tr$mean_car, col = "#2a78d6", lwd = 2)
  if (s == "Full sample")
    legend("bottomleft",
           legend = c("treated (published gsynth fits)",
                      "placebo mean (feventr gsynth)",
                      "placebo 95% block-bootstrap band"),
           col = c("#2a78d6", "#6b6b6b", "#d9d9d9"),
           lwd = c(2, 2, 8), lty = c(1, 2, 1), bty = "n", cex = 0.75)
}
dev.off()
cat("\nwrote", out_path("ma_placebo_gsynth.csv"), "and",
    out_path("ma_placebo_gsynth.png"), "\n")
