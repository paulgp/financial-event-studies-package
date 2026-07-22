# Long-run M&A effects against the placebo null — one figure combining
# the ma_refit_longrun.R treated paths (additive CATT, all three
# counterfactuals) with the ma_placebo_check.R placebo band (date-matched
# non-acquirers run through the identical SC pipeline). The grey band is
# the 95% circular block bootstrap of the PLACEBO path, so a treated line
# outside it rejects the no-effect null under announcement-time
# resampling. The placebo is fit with the SC pipeline: it is the exact
# null benchmark for the SC line and an approximate one for gsynth/apm
# (whose designs load less/more on the same microstructure bias).
#
# Reads output/ma_longrun_paths.csv and output/ma_placebo_check.csv;
# writes output/ma_longrun_placebo.png. Run from replication/ after both
# upstream scripts.
source("config.R")
suppressMessages(library(data.table))

lr <- fread(out_path("ma_longrun_paths.csv"))[metric == "arith"]
pl <- fread(out_path("ma_placebo_check.csv"))[metric == "arith" &
                                                role == "Placebo"]

cols <- c(Gsynth = "#2a78d6", SC = "#008300", APM = "#e87ba4")
png(out_path("ma_longrun_placebo.png"), width = 2600, height = 1000,
    res = 240)
par(mfrow = c(1, 3), mar = c(4, 4.2, 2.5, 5.5), mgp = c(2.4, 0.7, 0),
    las = 1)
ylim <- 100 * range(pl$lo, pl$hi, lr$mean_car, finite = TRUE)
for (s in c("Full sample", "Cash merger", "Stock merger")) {
  plot(NA, xlim = c(-1, 285), ylim = ylim, xlab = "Event day",
       ylab = if (s == "Full sample")
         "Mean additive CATT from day -1 (%)" else "",
       main = sprintf("%s (n=%s)", s,
                      format(lr[subsample == s & event_date == 250 &
                                  method == "Gsynth", n],
                             big.mark = ",")),
       cex.main = 0.95, font.main = 1, bty = "n", xaxt = "n")
  axis(1, at = c(0, 50, 100, 150, 200, 250))
  grid(nx = NA, ny = NULL, col = "#d8d8d8", lty = 1, lwd = 0.5)
  abline(h = 0, col = "#9a9a9a", lwd = 0.7)
  pb <- pl[subsample == s][order(event_date)]
  polygon(c(pb$event_date, rev(pb$event_date)),
          100 * c(pb$lo, rev(pb$hi)), col = "#dedede", border = NA)
  lines(pb$event_date, 100 * pb$mean_car, col = "#6b6b6b", lwd = 2,
        lty = 2)
  ends <- data.table(
    lab = c("Gsynth", "SC", "APM", "placebo"),
    col = c("#2b2b2b", "#2b2b2b", "#2b2b2b", "#6b6b6b"),
    y = c(vapply(c("Gsynth", "SC", "APM"), function(m)
      100 * lr[method == m & subsample == s & event_date == 250,
               mean_car], numeric(1)),
      100 * pb[event_date == 250, mean_car]))
  for (m in c("Gsynth", "SC", "APM")) {
    x <- lr[method == m & subsample == s][order(event_date)]
    lines(x$event_date, 100 * x$mean_car, col = cols[[m]], lwd = 2)
    points(250, ends[lab == m, y], pch = 16, col = cols[[m]], cex = 0.8)
  }
  # spread end labels to a minimum vertical gap (text only; markers stay
  # at the true values)
  setorder(ends, y)
  gap <- 0.032 * diff(ylim)
  for (i in seq_len(nrow(ends))[-1])
    ends$y[i] <- max(ends$y[i], ends$y[i - 1] + gap)
  ends[, txt := sprintf("%s %.1f%%", lab,
                        c(vapply(c("Gsynth", "SC", "APM"), function(m)
                          100 * lr[method == m & subsample == s &
                                     event_date == 250, mean_car],
                          numeric(1)),
                          100 * pb[event_date == 250, mean_car])[
                            match(lab, c("Gsynth", "SC", "APM",
                                         "placebo"))])]
  text(255, ends$y, ends$txt, pos = 4, cex = 0.72, col = ends$col,
       xpd = NA)
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
cat("wrote", out_path("ma_longrun_placebo.png"), "\n")
