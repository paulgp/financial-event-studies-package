# Table 1 pilot: 3 sims of Panels A (selection = "none") and D ("both"),
# end-to-end including gsynth. Validates
#   (1) per-sim bias/coverage vectors against the published per-sim
#       bias_estimates_*.csv saved in the PEAD_DinD Dropbox (deterministic
#       columns should match bit-exactly; gsynth point estimates near-exactly;
#       gsynth coverage has unseeded parametric-bootstrap noise),
#   (2) feventr::event_study(method = "mean") ATT == fixest coefficients,
#   (3) the 3-sim aggregation against published Panel A/D targets
#       (directional sanity only at 3 sims).
# Run from replication/: Rscript simulations/table1_pilot.R

source("config.R")
source("simulations/table1_common.R")

out_dir <- "simulations/sim_out"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pilot_sels <- c("none", "both")
pilot_sims <- 1:3

pub_file <- function(selection, sim_idx) {
  flag <- c(none = "FALSE_FALSE", assignment = "TRUE_FALSE",
            timing = "FALSE_TRUE", both = "TRUE_TRUE")[[selection]]
  dind("output/simulations/selection_2factors",
       sprintf("bias_estimates_500_10_240_one-shot-treatment_%s_%d.csv",
               flag, sim_idx))
}

for (sel in pilot_sels) {
  for (s in pilot_sims) {
    t0 <- Sys.time()
    run_one_sim_checkpointed(out_dir, sel, s)
    cat(sprintf("[%s sim %d] done in %.1fs\n", sel, s,
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))

    mine <- read.csv(sim_out_file(out_dir, sel, s))
    pubf <- pub_file(sel, s)
    if (file.exists(pubf)) {
      pub <- read.csv(pubf)
      for (cn in c("bias_simple", "bias_market", "bias_ff")) {
        d <- max(abs(mine[[cn]] - pub[[cn]]))
        cat(sprintf("  %-16s max |diff| vs published = %.3g %s\n",
                    cn, d, ifelse(d < 1e-10, "(bit-exact)", "")))
      }
      d <- max(abs(mine$bias_gsynth - pub$bias_gsynth))
      cat(sprintf("  %-16s max |diff| vs published = %.3g\n", "bias_gsynth", d))
      for (cn in c("coverage_simple", "coverage_market", "coverage_ff",
                   "coverage_gsynth")) {
        nd <- sum(mine[[cn]] != pub[[cn]])
        cat(sprintf("  %-16s rejections flipped vs published = %d/11\n", cn, nd))
      }
    } else cat("  (published per-sim file not found: ", pubf, ")\n")
  }
}

# --- feventr method="mean" cross-check (sim 1, panel A) -----------------------
sim <- simulate_events(selection = "none", seed = seed_for_sim(1))
long <- sim_long_panel(sim)
event_t <- N_PRE + 1L
fit <- feventr::event_study(
  long, unit = "id", time = "t", ret = "ret",
  treated = as.character(unique(long$id[long$treated])),
  event_time = event_t, method = "mean", returns = "simple",
  window = c(0, N_POST), est_window = c(-N_PRE, -1), se = "none"
)
post <- long[long$t >= event_t, ]
ct <- coeftable(feols(ret ~ -1 + i(t) + i(t, treated), data = post, vcov = "HC1"))
fe <- ct[grep(":treated", rownames(ct)), "Estimate"]
cat(sprintf("\nfeventr mean ATT vs feols coefficients: max |diff| = %.3g\n",
            max(abs(unname(fit$att) - unname(fe)))))

# --- 3-sim aggregation vs published panels (directional) ----------------------
targets <- read_target(1)
for (sel in pilot_sels) {
  panel <- names(SELECTIONS)[SELECTIONS == sel]
  dfs <- lapply(pilot_sims, function(s) read.csv(sim_out_file(out_dir, sel, s)))
  agg <- aggregate_panel(dfs)
  tg <- targets[targets$panel == panel, c("row", "col", "estimate")]
  cmp <- merge(agg, tg, by = c("row", "col"), suffixes = c("_pilot3", "_published"))
  cat("\n==", panel, "(3 sims vs published 50-sim values) ==\n")
  print(cmp[order(cmp$row, cmp$col), ], digits = 3, row.names = FALSE)
}
cat("\nPilot complete.\n")
