# Table 1 full run: 4 panels (none/assignment/timing/both) x 50 sims,
# seeds 1234..1283 (one per sim, identical to the published runs).
# Parallelized over sims with parallel::mclapply(mc.cores = 6); each sim is
# checkpointed to simulations/sim_out/bias_<selection>_<sim>.csv and skipped
# if already present (the pilot's panel A/D sims 1-3 are reused).
# Run from replication/: Rscript simulations/table1_full.R
# Final outputs: output/table1.csv and output/table1_comparison.csv.

source("config.R")
source("simulations/table1_common.R")

n_sims <- 50L
out_dir <- "simulations/sim_out"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

jobs <- expand.grid(sim = seq_len(n_sims), selection = unname(SELECTIONS),
                    stringsAsFactors = FALSE)
todo <- jobs[!file.exists(mapply(sim_out_file, out_dir,
                                 jobs$selection, jobs$sim)), ]
cat(sprintf("%d/%d sim cells to run\n", nrow(todo), nrow(jobs)))

if (nrow(todo) > 0) {
  res <- parallel::mclapply(seq_len(nrow(todo)), function(k) {
    sel <- todo$selection[k]; s <- todo$sim[k]
    t0 <- Sys.time()
    ok <- tryCatch({
      run_one_sim_checkpointed(out_dir, sel, s, gsynth_parallel = FALSE)
      TRUE
    }, error = function(e) {
      message(sprintf("[%s sim %d] ERROR: %s", sel, s, conditionMessage(e)))
      FALSE
    })
    cat(sprintf("[%s sim %02d] %s in %.1fs\n", sel, s,
                if (ok) "done" else "FAILED",
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    ok
  }, mc.cores = 6, mc.preschedule = FALSE)
  if (!all(unlist(res))) stop("some sims failed; re-run to retry")
}

# --- aggregate to the published table layout ---------------------------------
tab <- do.call(rbind, lapply(names(SELECTIONS), function(panel) {
  sel <- SELECTIONS[[panel]]
  dfs <- lapply(seq_len(n_sims), function(s) {
    fn <- sim_out_file(out_dir, sel, s)
    if (!file.exists(fn)) stop("missing checkpoint: ", fn)
    read.csv(fn)
  })
  agg <- aggregate_panel(dfs)
  data.frame(table = 1L, panel = panel, agg)
}))
write.csv(tab, out_path("table1.csv"), row.names = FALSE)

# --- cell-by-cell comparison vs targets ---------------------------------------
targets <- read_target(1)
cmp <- merge(targets[, c("table", "panel", "row", "col", "estimate", "units")],
             tab, by = c("table", "panel", "row", "col"),
             suffixes = c("_published", "_reproduced"), all.x = TRUE)
cmp$diff <- cmp$estimate_reproduced - cmp$estimate_published
# Published cells are rounded to 2 dp. Tolerance: +-0.1pp for bias cells
# (plus rounding); coverage cells are sim proportions (granularity 0.02 on
# the event day) with unseeded gsynth-bootstrap noise -> 0.025 + rounding.
cmp$tol <- ifelse(grepl("Coverage", cmp$col), 0.025 + 0.005, 0.1 + 0.005)
cmp$ok <- abs(cmp$diff) <= cmp$tol
cmp <- cmp[order(match(cmp$panel, names(SELECTIONS)),
                 match(cmp$row, MODEL_LABELS), cmp$col), ]
write.csv(cmp[, c("table", "panel", "row", "col", "units",
                  "estimate_published", "estimate_reproduced",
                  "diff", "tol", "ok")],
          out_path("table1_comparison.csv"), row.names = FALSE)

cat(sprintf("\n%d/%d cells within tolerance\n", sum(cmp$ok, na.rm = TRUE), nrow(cmp)))
print(cmp[!cmp$ok | is.na(cmp$ok),
          c("panel", "row", "col", "estimate_published", "estimate_reproduced", "diff")],
      digits = 3, row.names = FALSE)
cat("\nDone: output/table1.csv, output/table1_comparison.csv\n")
