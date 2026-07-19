# Compare the full-run gsynth/cfm/apm index-inclusion fits: pooled CATT
# paths on the cohorts all three methods fit (equal-weight across cohorts,
# cumulative from day -100, cross-cohort SEs), per-cohort CAR(+20)
# correlations, and the selected factor counts.
#
# Reads index_inclusion/{gsynth,cfm,apm}_out/ (build with
# table5_gsynth_full.R and cfm_apm_full.R). Writes
# output/ii_cfm_apm_catt_paths.csv, output/ii_cfm_apm_summary.csv, and the
# figure output/ii_cfm_apm_catt.png.
#
# Run from replication/: Rscript index_inclusion/cfm_apm_compare.R
suppressMessages(library(data.table))
source("config.R")

read_method <- function(dir) {
  fs <- list.files(file.path("index_inclusion", dir),
                   pattern = "^cohort_[0-9]+[.]csv$", full.names = TRUE)
  if (!length(fs)) stop("no cohort fits in index_inclusion/", dir)
  rbindlist(lapply(fs, fread), fill = TRUE)[!is.na(att)]
}
paths <- list(Gsynth = read_method("gsynth_out"), CFM = read_method("cfm_out"),
              APM = read_method("apm_out"))
ids <- Reduce(intersect, lapply(paths, function(d) unique(d$index_anndate)))
cat("cohorts:", paste(names(paths), vapply(paths, function(d)
  length(unique(d$index_anndate)), 0L), collapse = " | "),
  "| common", length(ids), "\n")

catt <- list(); car <- list()
for (m in names(paths)) {
  d <- paths[[m]][index_anndate %in% ids]
  setorder(d, index_anndate, event_date)
  d[, car := cumsum(att), by = index_anndate]
  catt[[m]] <- d[, .(catt = mean(car), se = sd(car) / sqrt(.N)),
                 by = event_date]
  car[[m]] <- d[event_date == 20, .(index_anndate, car20 = car)]
}

hz <- c(-50, -1, 0, 5, 20)
summ <- rbindlist(lapply(names(catt), function(m)
  cbind(method = m, catt[[m]][event_date %in% hz])))
summ[, `:=`(catt = 100 * catt, se = 100 * se)]
cat("\nPooled CATT (%, cumulative from -100; cross-cohort SEs):\n")
print(dcast(summ, method ~ event_date, value.var = c("catt", "se")),
      digits = 3)

cw <- Reduce(function(a, b) merge(a, b, by = "index_anndate"),
             lapply(names(car), function(m)
               setnames(copy(car[[m]]), "car20", m)))
P <- cor(as.matrix(cw[, -1]))
S <- cor(as.matrix(cw[, -1]), method = "spearman")
cat("\nPer-cohort CAR(+20) correlations (pearson below / spearman above):\n")
M <- P; M[upper.tri(M)] <- S[upper.tri(S)]
print(round(M, 3))
cat("\nSelected factor counts (identical selection inputs for cfm/apm):\n")
for (m in c("CFM", "APM")) {
  cat(m, ": ")
  print(table(paths[[m]][index_anndate %in% ids, r[1], by = index_anndate]$V1))
}

write.csv(rbindlist(lapply(names(catt), function(m)
  cbind(method = m, catt[[m]]))), out_path("ii_cfm_apm_catt_paths.csv"),
  row.names = FALSE)
write.csv(summ, out_path("ii_cfm_apm_summary.csv"), row.names = FALSE)

# figure: pooled CATT paths (categorical palette validated per dataviz skill)
cols <- c(Gsynth = "#2a78d6", CFM = "#008300", APM = "#e87ba4")
png(out_path("ii_cfm_apm_catt.png"), width = 2000, height = 1250, res = 240)
par(mar = c(4.2, 4.5, 2.5, 7), mgp = c(2.6, 0.7, 0), las = 1)
ylim <- 100 * range(sapply(catt, function(x) range(x$catt)))
plot(NA, xlim = c(-100, 30), ylim = ylim, xlab = "Event day",
     ylab = "Pooled CATT (%, cumulative from day -100)",
     main = sprintf("S&P 500 index inclusion: pooled CATT, %d common cohorts",
                    length(ids)),
     cex.main = 0.95, font.main = 1, bty = "n", xaxt = "n")
axis(1, at = seq(-100, 20, 20))
grid(nx = NA, ny = NULL, col = "#d8d8d8", lty = 1, lwd = 0.5)
abline(h = 0, col = "#9a9a9a", lwd = 0.7)
abline(v = 0, col = "#9a9a9a", lty = 2, lwd = 0.7)
text(0, ylim[2], "announcement", pos = 4, cex = 0.62, col = "#6b6b6b")
for (m in names(catt)) {
  x <- catt[[m]][order(event_date)]
  lines(x$event_date, 100 * x$catt, col = cols[[m]], lwd = 2)
  points(20, 100 * x[event_date == 20, catt], pch = 16, col = cols[[m]],
         cex = 0.8)
  text(21.5, 100 * x[event_date == 20, catt],
       sprintf("%s  %.2f%%", m, 100 * x[event_date == 20, catt]),
       pos = 4, cex = 0.72, col = "#2b2b2b", xpd = NA)
}
legend("topleft", legend = names(cols), col = cols, lwd = 2, bty = "n",
       cex = 0.75)
dev.off()
cat("\nwrote", out_path("ii_cfm_apm_catt_paths.csv"), "+ summary + figure\n")
