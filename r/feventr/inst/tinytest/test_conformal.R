# Conformal inference: brute-force references (method "mean" has an analytic
# refit-under-the-null), determinism, CI properties, and the warm-start
# invariances for the placebo loop and solver.

set.seed(11)
N0 <- 40; N1 <- 4; T0 <- 60; T1 <- 5
fac <- rnorm(T0 + T1, 0, 0.008)
Y <- matrix(rnorm((N0 + N1) * (T0 + T1), 0, 0.01), N0 + N1) +
  matrix(rep(fac, each = N0 + N1), N0 + N1)
tau <- 0.015
Y[(N0 + 1):(N0 + N1), (T0 + 1):(T0 + T1)] <-
  Y[(N0 + 1):(N0 + N1), (T0 + 1):(T0 + T1)] + tau
rownames(Y) <- c(paste0("d", 1:N0), paste0("x", 1:N1))
colnames(Y) <- 1:(T0 + T1)
long <- data.frame(id = rep(rownames(Y), times = ncol(Y)),
                   t = rep(1:ncol(Y), each = nrow(Y)),
                   ret = as.vector(Y))

fit_mean <- feventr::event_study(long, "id", "t", "ret",
                                 treated = paste0("x", 1:N1),
                                 event_time = T0 + 1, method = "mean",
                                 window = c(0, T1 - 1),
                                 est_window = c(-T0, -1),
                                 returns = "simple", se = "conformal")

# --- brute-force pointwise reference -------------------------------------
# for "mean", the refit-under-h0 counterfactual is the donor mean, so the
# residuals are (treated mean - donor mean), with h0 subtracted at the post
# slot; p(h0) = rank of the post |residual| among all T0+1 |residuals|
d <- colMeans(Y[(N0 + 1):(N0 + N1), ]) - colMeans(Y[1:N0, ])
p_hand <- function(j, h0) {
  e <- c(d[1:T0], d[T0 + j] - h0)
  mean(abs(e) >= abs(e[T0 + 1]) - 1e-12)
}
for (j in 1:T1)
  expect_equal(unname(fit_mean$se$p[j]), p_hand(j, 0), info = paste("p0, j =", j))

# CI endpoints sit where p_hand crosses alpha: just inside accepted, just
# outside rejected (p is a step function; probe a hair beyond the bisection
# tolerance on each side)
alpha <- 0.05
for (j in 1:T1) {
  lo <- fit_mean$se$ci[j, 1]; hi <- fit_mean$se$ci[j, 2]
  eps <- 0.05 * stats::sd(d[1:T0])
  expect_true(p_hand(j, lo + eps) > alpha, info = paste("inside lo, j =", j))
  expect_true(p_hand(j, lo - eps) <= alpha, info = paste("outside lo, j =", j))
  expect_true(p_hand(j, hi - eps) > alpha, info = paste("inside hi, j =", j))
  expect_true(p_hand(j, hi + eps) <= alpha, info = paste("outside hi, j =", j))
  expect_true(lo < fit_mean$att[j] && fit_mean$att[j] < hi)
}

# --- brute-force joint (moving-block) reference ---------------------------
p_joint_hand <- function(h0) {
  e <- abs(c(d[1:T0], d[(T0 + 1):(T0 + T1)] - h0))
  Tt <- T0 + T1
  blk <- (T0 + 1):Tt
  s <- vapply(1:Tt, function(sh) mean(e[((blk - 1 + sh) %% Tt) + 1]), 0)
  mean(s >= s[Tt] - 1e-12)
}
expect_equal(fit_mean$se$avg_p, p_joint_hand(0))
expect_true(fit_mean$se$avg_ci[1] < mean(fit_mean$att) &&
              mean(fit_mean$att) < fit_mean$se$avg_ci[2])

# deterministic: identical refit, no seed involved
fit2 <- feventr::event_study(long, "id", "t", "ret",
                             treated = paste0("x", 1:N1),
                             event_time = T0 + 1, method = "mean",
                             window = c(0, T1 - 1), est_window = c(-T0, -1),
                             returns = "simple", se = "conformal")
expect_identical(fit_mean$se, fit2$se)

# --- conformal with the synthetic-control engine --------------------------
fit_sc <- feventr::event_study(long, "id", "t", "ret",
                               treated = paste0("x", 1:N1),
                               event_time = T0 + 1, method = "sc",
                               window = c(0, T1 - 1), est_window = c(-T0, -1),
                               returns = "simple", se = "conformal")
expect_true(all(is.finite(fit_sc$se$ci)))
expect_true(all(fit_sc$se$ci[, 1] < fit_sc$att & fit_sc$att < fit_sc$se$ci[, 2]))
# true constant effect is covered by the joint CI
expect_true(fit_sc$se$avg_ci[1] < tau && tau < fit_sc$se$avg_ci[2])
# p-values live on the exact-permutation grid
expect_true(all(fit_sc$se$p >= 1 / (T0 + 1) - 1e-12 & fit_sc$se$p <= 1))
# methods understand conformal fits
expect_equal(unname(confint(fit_sc)), unname(fit_sc$se$ci))
expect_error(confint(fit_sc, level = 0.9), pattern = "level")
expect_stdout(print(fit_sc), pattern = "conformal CI")
expect_true(all(c("lower", "upper", "p") %in% names(summary(fit_sc))))
# guarded combinations
expect_error(feventr::event_study(long, "id", "t", "ret",
                                  treated = paste0("x", 1:N1),
                                  event_time = T0 + 1, method = "gsynth",
                                  window = c(0, T1 - 1),
                                  est_window = c(-T0, -1),
                                  returns = "simple", se = "conformal"),
             pattern = "conformal")

# --- warm starts ----------------------------------------------------------
# solver: warm start changes the path, not the solution
A <- t(Y[1:N0, 1:T0]); b <- colMeans(Y[(N0 + 1):(N0 + N1), 1:T0])
cold <- feventr::solve_simplex_ls(A, b)
warm <- feventr::solve_simplex_ls(A, b, w0 = cold$w)
expect_equal(warm$objective, cold$objective, tolerance = 1e-8)
# a junk warm start is projected/ignored, never breaks
junk <- feventr::solve_simplex_ls(A, b, w0 = rep(-1, N0))
expect_equal(junk$objective, cold$objective, tolerance = 1e-8)

# placebo: assignments are pre-drawn, so parallel draws are bit-identical
# to serial ones under the same seed (refits consume no RNG)
refit_sc <- function(Y, N0, T0, w0 = NULL) feventr:::eng_sc(Y, N0, T0, w0 = w0)
p_ser <- feventr:::inf_placebo(Y, N0, T0, n_treated = N1, refit = refit_sc,
                               reps = 20, seed = 7)
p_par <- feventr:::inf_placebo(Y, N0, T0, n_treated = N1, refit = refit_sc,
                               reps = 20, seed = 7, cores = 2L)
expect_identical(p_par$draws, p_ser$draws)
expect_identical(p_par$att, p_ser$att)
