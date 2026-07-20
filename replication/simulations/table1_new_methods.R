# Table 1 companion — the two latent-factor additions (cfm, apm) through the
# published simulation designs (table1_common.R: two-factor selection DGP,
# 50 sims per panel, seeds 1234..1283), fit via feventr::event_study() on
# the same recentered panels: window c(0, 10), est_window c(-239, -1).
#
# cfm: analytic SEs; per-day p-values are normal plug-in. Note the estimand
# difference: cfm targets the SYSTEMATIC effect (intercept + loading break
# over the whole 11-day window), so a one-day +3% jump is smeared to
# ~0.27%/day across the window by construction — expect large event-day
# "bias" and nonzero post-day "bias" with a roughly correct 11-day total.
# apm: imputes day-specific counterfactual means (like gsynth); weighted
# bootstrap with 200 draws to match the published gsynth nboots.
#
# Checkpointed per (panel, sim) into sim_out/newm_*.csv; safe to rerun.
# Writes output/table1_new_methods.csv and prints alongside the replicated
# Table 1 rows.
#
# Run from replication/: Rscript simulations/table1_new_methods.R
source("config.R")
source("simulations/table1_common.R")
suppressPackageStartupMessages(library(parallel))

NSIM <- 50L
NEW_LABELS <- c(cfm = "Causal Factor Model", apm = "Aggregated Proj. Matrix")

run_one_sim_new <- function(selection, sim_idx) {
  sim <- simulate_events(selection = selection, seed = seed_for_sim(sim_idx))
  long <- sim_long_panel(sim)
  event_t <- N_PRE + 1L
  true_param <- c(TAU, rep(0, N_POST))
  tr_ids <- as.character(sort(unique(long$id[long$treated])))
  base <- list(data = long, unit = "id", time = "t", ret = "ret",
               treated = tr_ids, event_time = event_t,
               window = c(0, N_POST), est_window = c(-N_PRE, -1),
               returns = "simple")
  f_c <- do.call(event_study, c(base, list(method = "cfm")))
  p_c <- 2 * pnorm(-abs(f_c$att / f_c$se$att))
  f_a <- do.call(event_study, c(base, list(method = "apm", reps = 200,
                                           seed = seed_for_sim(sim_idx))))
  p_a <- 2 * pnorm(-abs(f_a$att / f_a$se$att))
  data.frame(true_param = true_param,
             bias_cfm = unname(f_c$att) - true_param,
             coverage_cfm = unname(p_c) < 0.05,
             bias_apm = unname(f_a$att) - true_param,
             coverage_apm = unname(p_a) < 0.05)
}

new_out_file <- function(selection, sim_idx)
  file.path("simulations/sim_out", sprintf("newm_%s_%02d.csv", selection, sim_idx))

tasks <- expand.grid(sim_idx = seq_len(NSIM), selection = unname(SELECTIONS),
                     stringsAsFactors = FALSE)
st <- mclapply(seq_len(nrow(tasks)), function(i) {
  data.table::setDTthreads(1L)
  sel <- tasks$selection[i]; k <- tasks$sim_idx[i]
  fn <- new_out_file(sel, k)
  if (file.exists(fn)) return("skip")
  res <- tryCatch(run_one_sim_new(sel, k), error = function(e) e)
  if (inherits(res, "error")) return(paste0("failed: ", conditionMessage(res)))
  tmp <- paste0(fn, ".tmp")
  write.csv(res, tmp, row.names = FALSE)
  file.rename(tmp, fn)
  "done"
}, mc.cores = 6L)
print(table(unlist(st)))
stopifnot(!any(grepl("failed", unlist(st))))

# aggregate with the published formulas (see aggregate_panel)
aggregate_new <- function(sim_dfs) {
  rows <- list()
  for (m in c("cfm", "apm")) {
    b <- sapply(sim_dfs, function(d) d[[paste0("bias_", m)]])
    cv <- sapply(sim_dfs, function(d) as.numeric(d[[paste0("coverage_", m)]]))
    rows[[m]] <- c(
      `All Periods: E(Bias)`       = mean(colSums(b)) / 11 * 100,
      `All Periods: MAD`           = mean(abs(colSums(b))) / 11 * 100,
      `All Periods: RMSE`          = sqrt(mean(colSums(b^2))) * 100,
      `Treated Periods: E(Bias)`   = mean(b[1, ]) * 100,
      `Treated Periods: MAD`       = mean(abs(b[1, ])) * 100,
      `Treated Periods: Coverage`  = mean(cv[1, ]),
      `Untreated Periods: E(Bias)` = mean(colSums(b[-1, , drop = FALSE])) / 10 * 100,
      `Untreated Periods: MAD`     = mean(abs(colSums(b[-1, , drop = FALSE]))) / 10 * 100,
      `Untreated Periods: Coverage` = mean(colSums(cv[-1, , drop = FALSE])) / 10
    )
  }
  do.call(rbind, lapply(c("cfm", "apm"), function(m)
    data.frame(row = NEW_LABELS[[m]], col = names(rows[[m]]),
               estimate = unname(rows[[m]]))))
}

out <- list()
for (p in seq_along(SELECTIONS)) {
  sel <- unname(SELECTIONS)[p]
  dfs <- lapply(seq_len(NSIM), function(k) read.csv(new_out_file(sel, k)))
  out[[p]] <- cbind(panel = names(SELECTIONS)[p], aggregate_new(dfs))
}
new <- do.call(rbind, out)
write.csv(new, out_path("table1_new_methods.csv"), row.names = FALSE)

old <- read.csv(out_path("table1.csv"))
comb <- rbind(old[, c("panel", "row", "col", "estimate")],
              new[, c("panel", "row", "col", "estimate")])
for (p in unique(comb$panel)) {
  cat("\n==", p, "==\n")
  w <- reshape(comb[comb$panel == p, c("row", "col", "estimate")],
               idvar = "row", timevar = "col", direction = "wide")
  names(w) <- sub("estimate[.]", "", names(w))
  print(w, digits = 2, row.names = FALSE)
}
