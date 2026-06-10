# Shared machinery for Table 1 (simulation bias/coverage) replication.
#
# DGP: feventr::simulate_events(), which mirrors the published
# code/simulations_selection_TL.R sim_data() bit-exactly (same RNG
# consumption order) for seeds 1234..1283.
#
# Estimators mirror the published estimate_bias()/estimate_excess_ret_models()
# exactly:
#   - simple means: fixest::feols(ret ~ -1 + i(t) + i(t, treated), vcov="HC1")
#     on the post panel (event day + 10). feventr method="mean" produces the
#     identical ATT coefficients (cross-checked in the pilot); fixest is used
#     here because the published coverage column comes from its HC1 p-values.
#   - CAPM ("market"): per-firm OLS of exret = ret - rf on MKT over the 239
#     pre-event days, abnormal = ret - rf - fitted; per-period mean over
#     treated via feols(~ -1 + i(t), HC1).
#   - correct factor ("ff"): per-firm OLS of RAW ret on MKT + SMB (note the
#     published rf asymmetry), abnormal = ret - rf - fitted.
#   - gsynth: gsynth::gsynth(ret ~ D, force = "unit", r = c(1, 40), CV = TRUE,
#     se = TRUE, inference = "parametric", nboots = 200) — the exact published
#     call (CV/nboots are gsynth defaults there). We call gsynth directly (not
#     through feventr's engine) because the published numbers come from this
#     call verbatim, including its internal CV; parallel = FALSE because the
#     full run parallelizes over sims with mclapply.

suppressPackageStartupMessages({
  library(feventr)
  library(fixest)
  library(gsynth)
})

N_PRE  <- 239L  # estimation days before the event day
N_POST <- 10L   # post-event days after the event day (11-day reported window)
TAU    <- 0.03

SELECTIONS <- c(
  "Panel A: Random Assignment + Random Timing"      = "none",
  "Panel B: Assignment Selection + Random Timing"   = "assignment",
  "Panel C: Random Assignment + Timing Selection"   = "timing",
  "Panel D: Assignment Selection + Timing Selection" = "both"
)

MODEL_LABELS <- c(simple = "Simple Means", market = "CAPM",
                  ff = "Correct Factor Structure", gsynth = "Gsynth (PCA)")

seed_for_sim <- function(sim_idx) 1234L + sim_idx - 1L

# Build the recentered long panel for one simulation: t = 1..250 with the
# event day at t = 240 (original: filter to [event-239, event+10], recenter).
sim_long_panel <- function(sim) {
  ev <- sim$event_time
  d <- sim$data[sim$data$t >= ev - N_PRE & sim$data$t <= ev + N_POST, ]
  d$t <- d$t - (ev - N_PRE) + 1L          # 1..250, event at 240
  d$treated <- sim$betas$treated[d$id]
  f <- sim$factors[sim$factors$t >= ev - N_PRE & sim$factors$t <= ev + N_POST, ]
  f$t <- f$t - (ev - N_PRE) + 1L
  merge(d, f, by = "t")
}

# Per-firm OLS coefficients of y ~ X (with intercept) over the pre period.
firm_ols <- function(pre, yname, xnames) {
  ids <- sort(unique(pre$id))
  out <- matrix(NA_real_, length(ids), length(xnames) + 1L,
                dimnames = list(ids, c("alpha", xnames)))
  for (j in seq_along(ids)) {
    di <- pre[pre$id == ids[j], ]
    X <- cbind(1, as.matrix(di[, xnames, drop = FALSE]))
    out[j, ] <- qr.coef(qr(X), di[[yname]])
  }
  out
}

# One simulation -> the published bias_estimates data.frame (11 rows).
run_one_sim <- function(selection, sim_idx, gsynth_parallel = FALSE) {
  sim <- simulate_events(selection = selection, seed = seed_for_sim(sim_idx))
  long <- sim_long_panel(sim)
  event_t <- N_PRE + 1L                                # 240
  true_param <- c(TAU, rep(0, N_POST))

  # --- simple means (post panel, all firms) ---------------------------------
  post <- long[long$t >= event_t, ]
  fit_simple <- feols(ret ~ -1 + i(t) + i(t, treated), data = post, vcov = "HC1")
  ct <- coeftable(fit_simple)
  tr_rows <- grep(":treated", rownames(ct))
  est_simple <- ct[tr_rows, "Estimate"]
  p_simple   <- ct[tr_rows, "Pr(>|t|)"]

  # --- CAPM and correct-factor abnormal returns -----------------------------
  pre <- long[long$t < event_t, ]
  pre$exret <- pre$ret - pre$rf
  b_capm <- firm_ols(pre, "exret", "mktrf")            # exret ~ MKT
  b_ff   <- firm_ols(pre, "ret",  c("mktrf", "smb"))   # RAW ret ~ MKT + SMB
  post_tr <- post[post$treated, ]
  idx <- match(post_tr$id, as.integer(rownames(b_capm)))
  post_tr$ar_capm <- post_tr$ret - post_tr$rf -
    (b_capm[idx, "alpha"] + b_capm[idx, "mktrf"] * post_tr$mktrf)
  post_tr$ar_ff <- post_tr$ret - post_tr$rf -
    (b_ff[idx, "alpha"] + b_ff[idx, "mktrf"] * post_tr$mktrf +
       b_ff[idx, "smb"] * post_tr$smb)
  ct_m <- coeftable(feols(ar_capm ~ -1 + i(t), data = post_tr, vcov = "HC1"))
  ct_f <- coeftable(feols(ar_ff   ~ -1 + i(t), data = post_tr, vcov = "HC1"))

  # --- gsynth (published call; CV=TRUE/nboots=200 are its defaults) ---------
  gdat <- long
  gdat$treated <- as.numeric(gdat$treated)
  gdat$time <- gdat$t - event_t
  gdat$D <- gdat$treated * (gdat$t >= event_t)
  gfit <- suppressWarnings(gsynth(
    ret ~ D, data = gdat, index = c("id", "time"),
    force = "unit", r = c(1, 40), CV = TRUE,
    se = TRUE, inference = "parametric", nboots = 200,
    parallel = gsynth_parallel
  ))
  att_g <- gfit$att[as.character(0:N_POST)]
  # est.att rownames are shifted +1 vs att names (gsynth labels the first
  # treated period "1"); rows 1..11 = event day + 10 post, identical to the
  # published positional extraction est.att[240:250, ].
  p_g   <- gfit$est.att[as.character(seq_len(N_POST + 1L)), "p.value"]

  data.frame(
    true_param      = true_param,
    bias_simple     = unname(est_simple) - true_param,
    coverage_simple = unname(p_simple) < 0.05,
    bias_market     = unname(ct_m[, "Estimate"]) - true_param,
    bias_ff         = unname(ct_f[, "Estimate"]) - true_param,
    coverage_market = unname(ct_m[, "Pr(>|t|)"]) < 0.05,
    coverage_ff     = unname(ct_f[, "Pr(>|t|)"]) < 0.05,
    bias_gsynth     = unname(att_g) - true_param,
    coverage_gsynth = unname(p_g) < 0.05
  )
}

sim_out_file <- function(out_dir, selection, sim_idx)
  file.path(out_dir, sprintf("bias_%s_%02d.csv", selection, sim_idx))

run_one_sim_checkpointed <- function(out_dir, selection, sim_idx, ...) {
  fn <- sim_out_file(out_dir, selection, sim_idx)
  if (file.exists(fn)) return(invisible(fn))
  res <- run_one_sim(selection, sim_idx, ...)
  tmp <- paste0(fn, ".tmp")
  write.csv(res, tmp, row.names = FALSE)
  file.rename(tmp, fn)
  invisible(fn)
}

# Aggregation (mirrors make_tables_selection_TL.R lines 73-153):
#   All:       E(Bias) = mean_sims(sum_t bias)/11*100; MAD = mean(|sum_t|)/11*100;
#              RMSE = sqrt(mean_sims(sum_t bias^2))*100   (no /11)
#   Treated:   t = 1 (event day): mean*100, MAD*100, coverage = mean(0/1)
#   Untreated: sums over t = 2..11, /10; coverage = mean(sum 0/1)/10
aggregate_panel <- function(sim_dfs) {
  models <- c("simple", "market", "ff", "gsynth")
  rows <- list()
  for (m in models) {
    b <- sapply(sim_dfs, function(d) d[[paste0("bias_", m)]])      # 11 x nsim
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
  do.call(rbind, lapply(models, function(m)
    data.frame(row = MODEL_LABELS[[m]], col = names(rows[[m]]),
               estimate = unname(rows[[m]]))))
}
