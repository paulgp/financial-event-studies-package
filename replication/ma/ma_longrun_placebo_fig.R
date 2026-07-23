# Long-run M&A effects against method-matched placebo nulls — the
# ma_refit_longrun.R treated paths (additive CATT) overlaid on TWO
# placebo bands, each from date-matched non-acquirers run through the
# corresponding pipeline (identical placebo units across methods, by
# seeding): the SC placebo band (ma_placebo_check.R) and the gsynth
# placebo band (ma_placebo_gsynth.R). Each treated line should be judged
# against its own band; APM has no placebo run. The contrast IS the
# result: the SC placebo drifts to -30% by +250d (noise-selection bias)
# while the gsynth placebo stays within a point of zero, so gsynth's
# long-run drift is signal while SC's is mostly design bias.
#
# Reads output/ma_longrun_paths.csv, output/ma_placebo_check.csv and
# output/ma_placebo_gsynth.csv; writes output/ma_longrun_placebo.png.
# Run from replication/ after the upstream scripts.
source("config.R")
suppressMessages(library(data.table))

lr <- fread(out_path("ma_longrun_paths.csv"))[metric == "arith"]
plc <- rbind(
  fread(out_path("ma_placebo_check.csv"))[metric == "arith" &
    role == "Placebo"][, method := "SC"],
  fread(out_path("ma_placebo_gsynth.csv"))[metric == "arith" &
    role == "Placebo"][, method := "Gsynth"],
  fill = TRUE)

cols <- c(Gsynth = "#2a78d6", SC = "#008300", APM = "#e87ba4")
pcols <- c(Gsynth = "#7aa7dd", SC = "#66a866")

png(out_path("ma_longrun_placebo.png"), width = 2600, height = 1000,
    res = 240)
par(mfrow = c(1, 3), mar = c(4, 4.2, 2.5, 5.5), mgp = c(2.4, 0.7, 0),
    las = 1)
ylim <- 100 * range(plc$lo, plc$hi, lr$mean_car, finite = TRUE)
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
  for (m in c("SC", "Gsynth")) {
    pb <- plc[method == m & subsample == s][order(event_date)]
    polygon(c(pb$event_date, rev(pb$event_date)),
            100 * c(pb$lo, rev(pb$hi)),
            col = grDevices::adjustcolor(cols[[m]], 0.13), border = NA)
    lines(pb$event_date, 100 * pb$mean_car, col = pcols[[m]], lwd = 1.8,
          lty = 2)
  }
  ends <- data.table(
    lab = c("Gsynth", "SC", "APM", "SC placebo", "Gsynth placebo"),
    col = c("#2b2b2b", "#2b2b2b", "#2b2b2b", pcols[["SC"]],
            pcols[["Gsynth"]]),
    yv = c(vapply(c("Gsynth", "SC", "APM"), function(m)
      100 * lr[method == m & subsample == s & event_date == 250,
               mean_car], numeric(1)),
      100 * plc[method == "SC" & subsample == s & event_date == 250,
                mean_car],
      100 * plc[method == "Gsynth" & subsample == s & event_date == 250,
                mean_car]))
  for (m in c("Gsynth", "SC", "APM")) {
    x <- lr[method == m & subsample == s][order(event_date)]
    lines(x$event_date, 100 * x$mean_car, col = cols[[m]], lwd = 2)
    points(250, ends[lab == m, yv], pch = 16, col = cols[[m]], cex = 0.8)
  }
  # spread end labels to a minimum vertical gap (text only; markers stay
  # at the true values)
  ends[, y := yv]
  setorder(ends, y)
  gap <- 0.032 * diff(ylim)
  for (i in seq_len(nrow(ends))[-1])
    ends$y[i] <- max(ends$y[i], ends$y[i - 1] + gap)
  text(255, ends$y, sprintf("%s %.1f%%", ends$lab, ends$yv), pos = 4,
       cex = 0.68, col = ends$col, xpd = NA)
  if (s == "Full sample")
    legend("bottomleft",
           legend = c("Gsynth", "SC", "APM (no placebo run)",
                      "placebo means (dashed) with 95%",
                      "block-bootstrap bands, method-matched",
                      "(same placebo units in both)"),
           col = c(cols, NA, NA, NA),
           lwd = c(2, 2, 2, NA, NA, NA), bty = "n", cex = 0.72)
}
dev.off()
cat("wrote", out_path("ma_longrun_placebo.png"), "\n")
