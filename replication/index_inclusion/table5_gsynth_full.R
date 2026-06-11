# Table 5 Gsynth column — FULL RUN over all 635 announcement-date cohorts.
#
# Stage 1 extracts per-cohort slices of the 5 panel .dta files (44GB total)
# into the cache outside the repo (one streaming pass per panel; skipped when
# already cached). Stage 2 fits feventr gsynth per cohort (force="unit",
# r=c(1,10), window=c(-100,20), est_window=c(-280,-101), se="none") with
# parallel::mclapply(mc.cores=8), checkpointing one CSV per cohort to
# index_inclusion/gsynth_out/ (already-done cohorts are skipped, so the run
# is restartable; the pilot's cohorts 1-15 count). Stage 3 aggregates the
# day +1 effects into the Table 5 Gsynth rows and rewrites
# output/table5.csv + output/table5_comparison.csv.
#
# Published-vintage note: the original run covered 613/635 cohorts (21 gsynth
# failures + cohort 635 skipped by a loop bug). The comparison Gsynth column
# is therefore aggregated over cohorts present in the saved sc_ii_siblis.dta
# (published vintage); the all-successful-cohorts variant is reported
# alongside in output/table5_gsynth_run_summary.csv.
#
# Run from replication/:
#   nohup Rscript index_inclusion/table5_gsynth_full.R > output/table5_gsynth_full.log 2>&1 &

suppressMessages({library(haven); library(data.table); library(parallel)})
source("config.R")
source("index_inclusion/betas_common.R")
source("index_inclusion/gsynth_cohort.R")

outdir <- "index_inclusion/gsynth_out"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
cache <- feventr_cache("cohorts")
python <- Sys.getenv("FEVENTR_PYTHON", unname(Sys.which("python3")))

# ---- stage 1: per-cohort extraction (one pass per panel, skipped if cached) --
panels <- list(c("panel_ii_1_272.dta", 1, 272),
               c("panel_ii_273_300.dta", 273, 300),
               c("panel_ii_300_400.dta", 301, 400),
               c("panel_ii_400_500.dta", 401, 500),
               c("panel_ii_500_635.dta", 501, 635))
for (p in panels) {
  rc <- system2(python,
                c("index_inclusion/extract_cohorts.py",
                  shQuote(ii_work(file.path("siblis_anndates", p[1]))), p[2], p[3],
                  shQuote(cache)))
  if (rc != 0) stop("extraction failed for ", p[1])
}

# ---- stage 2: per-cohort gsynth with checkpointing ---------------------------
cohorts <- 1:635
st <- mclapply(cohorts, run_cohort, outdir = outdir, cache = cache,
               mc.cores = 8, mc.preschedule = FALSE)
status <- vapply(st, function(x) if (is.character(x)) x else "error", "")
cat("cohort status counts:\n"); print(table(sub(":.*", "", status)))
failed <- cohorts[!status %in% c("done")]
if (length(failed)) cat("not-done cohorts:", paste(failed, collapse = " "), "\n")

# ---- stage 3: aggregate day +1 effects into Table 5 --------------------------
res <- rbindlist(lapply(cohorts, function(i) {
  f <- file.path(outdir, sprintf("cohort_%d.csv", i))
  if (file.exists(f)) fread(f)[event_date == 1] else NULL
}))
ev <- as.data.table(read_dta(ii_work("include_event_date_siblis.dta")))
ann <- unique(ev[, .(anndate)])[order(anndate)][, index_anndate := .I]
res <- merge(res, ann, by = "index_anndate")
res[, group := decade_group(year(anndate))]

sc <- as.data.table(read_dta(ii_work("sc_ii_siblis.dta")))
pub_cohorts <- unique(sc$index_anndate)

agg <- function(d) d[!is.na(group),
                     .(Gsynth = 100 * weighted.mean(att, n_treat),
                       n_events = sum(n_treat), n_cohorts = .N), by = group]
g_pub <- agg(res[index_anndate %in% pub_cohorts])  # published vintage
g_all <- agg(res)                                  # every successful cohort
summ <- merge(g_pub, g_all, by = "group", suffixes = c("_pubvintage", "_all"))
dec_lab <- c("1980-1989", "1990-1999", "2000-2009", "2010-2020")
summ[, row := dec_lab[group]]
print(summ[order(group)], digits = 4)
write.csv(summ[order(group)], out_path("table5_gsynth_run_summary.csv"),
          row.names = FALSE)

# rewrite table5.csv / table5_comparison.csv with the re-estimated Gsynth rows
tab <- as.data.table(read.csv(out_path("table5.csv")))
g <- g_pub[, .(row = dec_lab[group], estimate_new = Gsynth)]
tab <- merge(tab, g, by = "row", all.x = TRUE)
tab[col == "Gsynth",
    `:=`(estimate = estimate_new,
         provisional = "feventr full run (published-vintage cohorts)")]
tab[, estimate_new := NULL]
col_order <- c("Diff-in-Means", "Market", "CAPM", "FF3F", "Gsynth")
tab <- tab[order(row, match(col, col_order))]
write.csv(tab, out_path("table5.csv"), row.names = FALSE)

tg <- as.data.table(read_target(5))
cmp <- merge(tg[, .(row, col, estimate_pub = estimate)],
             tab[, .(row, col, estimate_rep = estimate, provisional)],
             by = c("row", "col"))
cmp[, diff := estimate_rep - estimate_pub]
cmp[, ok := abs(diff) <= 0.1 + 0.005]
cmp <- cmp[order(row, match(col, col_order))]
print(cmp, digits = 4)
cat(sprintf("\nTable 5 (final): %d/%d cells within tolerance\n",
            sum(cmp$ok), nrow(cmp)))
write.csv(cmp, out_path("table5_comparison.csv"), row.names = FALSE)
