# Post-decimalization variant of the long-run + placebo-bounds figure:
# identical machinery to ma_refit_longrun.R / ma_placebo_check.R /
# ma_longrun_placebo_fig.R but restricted to deals announced 2001+,
# where the bid-ask-bounce bias in the synthetic side is small
# (placebo mean att -4bp/day vs -19bp/day pre-decimalization).
# Additive CATT only. Writes output/ma_longrun_placebo_2001.png and
# output/ma_longrun_placebo_2001.csv.
#
# Run from replication/: Rscript ma/ma_longrun_placebo_2001.R
source("config.R")
suppressMessages({library(haven); library(data.table)})
ma_work <- function(...) dind("M&A", "data", "work", ...)
ma_out  <- function(...) dind("M&A", "output", ...)

B_BOOT <- 1000L
BLOCK_MONTHS <- 18L
ERA_MIN <- 2001L
set.seed(20260721)

## ---- deal sample (as ma_refit_longrun.R) ------------------------------------
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
deals <- deals[year(as.Date(ann_tdate)) >= ERA_MIN]
cat("deal sample,", ERA_MIN, "+:", nrow(deals), "\n")

## ---- per-deal additive CATT paths -------------------------------------------
read_refit <- function(m, tag_col = "permno") {
  fs <- list.files(file.path("ma", "ma_refit_out", m),
                   pattern = "^cohort_[0-9]+[.]csv$", full.names = TRUE)
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

read_gsynth <- function() {
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

paths <- list(Gsynth = read_gsynth(), SC = read_refit("sc"),
              APM = read_refit("apm"),
              Placebo = read_refit("placebo_sc", "matched_permno"))

SUBS <- function(d) list("Full sample"  = rep(TRUE, nrow(d)),
                         "Cash merger"  = d$pct_cash == 100,
                         "Stock merger" = d$pct_stk == 100)
months_all <- sort(unique(substr(
  unlist(lapply(paths, function(d) d$ann_tdate)), 1, 7)))
nm <- length(months_all)
nb <- ceiling(nm / BLOCK_MONTHS)
cat(sprintf("bootstrap: %d months, %d blocks of %d, B = %d\n",
            nm, nb, BLOCK_MONTHS, B_BOOT))
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
      method = m, subsample = s, event_date = horizons,
      mean_car = colSums(S) / colSums(N),
      boot_se = apply(reps, 2, sd),
      lo = apply(reps, 2, quantile, 0.025),
      hi = apply(reps, 2, quantile, 0.975),
      n = as.integer(colSums(N)))
  }
  rbindlist(out)
}
pp <- rbindlist(Map(pool_boot, paths, names(paths)))
write.csv(pp, out_path("ma_longrun_placebo_2001.csv"), row.names = FALSE)

hz <- c(1, 21, 63, 126, 250)
cat("\n2001+ additive CATT (%), mean (block-bootstrap SE):\n")
print(dcast(pp[event_date %in% hz,
               .(method, subsample, event_date,
                 cell = sprintf("%6.2f (%.2f)", 100 * mean_car,
                                100 * boot_se))],
            method + subsample ~ event_date, value.var = "cell"),
      row.names = FALSE)

## ---- figure -------------------------------------------------------------------
cols <- c(Gsynth = "#2a78d6", SC = "#008300", APM = "#e87ba4")
png(out_path("ma_longrun_placebo_2001.png"), width = 2600, height = 1000,
    res = 240)
par(mfrow = c(1, 3), mar = c(4, 4.2, 2.5, 5.5), mgp = c(2.4, 0.7, 0),
    las = 1)
pb_all <- pp[method == "Placebo"]
tr_all <- pp[method != "Placebo"]
ylim <- 100 * range(pb_all$lo, pb_all$hi, tr_all$mean_car, finite = TRUE)
for (s in c("Full sample", "Cash merger", "Stock merger")) {
  plot(NA, xlim = c(-1, 285), ylim = ylim, xlab = "Event day",
       ylab = if (s == "Full sample")
         "Mean additive CATT from day -1 (%)" else "",
       main = sprintf("%s, %d+ (n=%s)", s, ERA_MIN,
                      format(tr_all[subsample == s & event_date == 250 &
                                      method == "Gsynth", n],
                             big.mark = ",")),
       cex.main = 0.95, font.main = 1, bty = "n", xaxt = "n")
  axis(1, at = c(0, 50, 100, 150, 200, 250))
  grid(nx = NA, ny = NULL, col = "#d8d8d8", lty = 1, lwd = 0.5)
  abline(h = 0, col = "#9a9a9a", lwd = 0.7)
  pb <- pb_all[subsample == s][order(event_date)]
  polygon(c(pb$event_date, rev(pb$event_date)),
          100 * c(pb$lo, rev(pb$hi)), col = "#dedede", border = NA)
  lines(pb$event_date, 100 * pb$mean_car, col = "#6b6b6b", lwd = 2,
        lty = 2)
  ends <- data.table(
    lab = c("Gsynth", "SC", "APM", "placebo"),
    col = c("#2b2b2b", "#2b2b2b", "#2b2b2b", "#6b6b6b"),
    yv = c(vapply(c("Gsynth", "SC", "APM"), function(m)
      100 * tr_all[method == m & subsample == s & event_date == 250,
                   mean_car], numeric(1)),
      100 * pb[event_date == 250, mean_car]))
  for (m in c("Gsynth", "SC", "APM")) {
    x <- tr_all[method == m & subsample == s][order(event_date)]
    lines(x$event_date, 100 * x$mean_car, col = cols[[m]], lwd = 2)
    points(250, ends[lab == m, yv], pch = 16, col = cols[[m]], cex = 0.8)
  }
  ends[, y := yv]
  setorder(ends, y)
  gap <- 0.032 * diff(ylim)
  for (i in seq_len(nrow(ends))[-1])
    ends$y[i] <- max(ends$y[i], ends$y[i - 1] + gap)
  text(255, ends$y, sprintf("%s %.1f%%", ends$lab, ends$yv), pos = 4,
       cex = 0.72, col = ends$col, xpd = NA)
  if (s == "Full sample")
    legend("bottomleft",
           legend = c("Gsynth", "SC", "APM",
                      "placebo mean (SC pipeline)",
                      "placebo 95% block-bootstrap band",
                      "(18-month announcement blocks)"),
           col = c(cols, "#6b6b6b", "#dedede", NA),
           lwd = c(2, 2, 2, 2, 8, NA),
           lty = c(1, 1, 1, 2, 1, NA), bty = "n", cex = 0.72)
}
dev.off()
cat("\nwrote", out_path("ma_longrun_placebo_2001.png"), "and",
    out_path("ma_longrun_placebo_2001.csv"), "\n")
