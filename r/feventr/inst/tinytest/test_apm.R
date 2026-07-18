# Aggregated projection matrix (Lei & Ross) via the GitHub apm package

if (!requireNamespace("apm", quietly = TRUE))
  exit_file("apm not available")

# strong-factor DGP with a single treated unit: 60 donors, 2 factors,
# constant tau on the event window
set.seed(3)
N0 <- 60; T0 <- 40; T1 <- 6; Tn <- T0 + T1; tau <- 0.02
Fm <- matrix(rnorm(Tn * 2, 0, 0.01), Tn, 2)
Lam <- matrix(rnorm((N0 + 1) * 2, 1, 0.5), N0 + 1, 2)
ret <- Lam %*% t(Fm) + matrix(rnorm((N0 + 1) * Tn, 0, 0.003), N0 + 1, Tn)
ret[N0 + 1, (T0 + 1):Tn] <- ret[N0 + 1, (T0 + 1):Tn] + tau
long <- data.frame(id = rep(seq_len(N0 + 1), times = Tn),
                   t = rep(seq_len(Tn), each = N0 + 1),
                   ret = as.vector(ret))

fit <- feventr::event_study(long, "id", "t", "ret",
                            treated = as.character(N0 + 1), event_time = T0 + 1,
                            method = "apm", window = c(0, 5),
                            est_window = c(-40, -1), returns = "simple",
                            r = 2, reps = 100, seed = 7)

# recovers the injected constant effect up to idiosyncratic noise
# (deterministic given the seeds)
expect_true(max(abs(fit$att - tau)) < 0.015)
expect_equal(fit$diagnostics$info$r, 2L)
expect_equal(fit$diagnostics$info$n_cohorts, 2L)
expect_true(is.finite(fit$diagnostics$info$pre_rmse))

# native inference: multinomial weighted bootstrap with stored draws
expect_equal(fit$se$method, "bootstrap")
expect_equal(fit$se$reps, 100L)
expect_equal(dim(fit$se$draws), c(100L, 6L))
expect_equal(unname(fit$se$att), unname(apply(fit$se$draws, 2, sd)),
             tolerance = 1e-12)
expect_equal(fit$att_avg_se, sd(rowMeans(fit$se$draws)), tolerance = 1e-12)
expect_equal(dim(vcov(fit)), c(6L, 6L))

# same seed reproduces the bootstrap exactly
fit2 <- feventr::event_study(long, "id", "t", "ret",
                             treated = as.character(N0 + 1), event_time = T0 + 1,
                             method = "apm", window = c(0, 5),
                             est_window = c(-40, -1), returns = "simple",
                             r = 2, reps = 100, seed = 7)
expect_equal(fit$se$att, fit2$se$att)

# eigenvalue-ratio selection with the default r range finds the true rank
fit_ah <- feventr::event_study(long, "id", "t", "ret",
                               treated = as.character(N0 + 1), event_time = T0 + 1,
                               method = "apm", window = c(0, 5),
                               est_window = c(-40, -1), returns = "simple",
                               se = "none")
expect_equal(fit_ah$diagnostics$info$r, 2L)
expect_equivalent(unname(fit_ah$att), unname(fit$att), tolerance = 1e-10)

# placebo inference reruns the engine on donors
fit_pl <- feventr::event_study(long, "id", "t", "ret",
                               treated = as.character(N0 + 1), event_time = T0 + 1,
                               method = "apm", window = c(0, 5),
                               est_window = c(-40, -1), returns = "simple",
                               r = 2, se = "placebo", reps = 9, seed = 1)
expect_equal(fit_pl$se$method, "placebo")
expect_equal(length(fit_pl$se$att), 6L)

# guards
expect_error(feventr::event_study(long, "id", "t", "ret",
                                  treated = as.character(N0 + 1),
                                  event_time = T0 + 1, method = "apm",
                                  window = c(0, 5), est_window = c(-40, -1),
                                  returns = "simple", se = "tstat"),
             pattern = "tstat")
expect_error(feventr::event_study(long, "id", "t", "ret",
                                  treated = as.character(N0 + 1),
                                  event_time = T0 + 1, method = "apm",
                                  window = c(0, 5), est_window = c(-40, -1),
                                  returns = "simple", se = "conformal"),
             pattern = "conformal")
expect_error(feventr::event_study(long, "id", "t", "ret",
                                  treated = as.character(N0 + 1),
                                  event_time = T0 + 1, method = "apm",
                                  window = c(0, 5), est_window = c(-40, -1),
                                  returns = "simple", r = 45),
             pattern = "between 1 and")

# batch mode with the minimal events table
b <- feventr::event_study_batch(long, "id", "t", "ret",
                                events = data.frame(unit = as.character(N0 + 1),
                                                    event_time = T0 + 1),
                                method = "apm", window = c(0, 5),
                                est_window = c(-40, -1), returns = "simple",
                                r = 2)
expect_equivalent(unname(b$att), unname(fit$att), tolerance = 1e-10)
