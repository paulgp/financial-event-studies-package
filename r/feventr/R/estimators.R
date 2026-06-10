# Estimator engines: pure functions on matrices --------------------------------
#
# Shared contract: every engine takes Y ((N0 + N1) x T, controls in rows
# 1..N0, estimation periods in columns 1..T0) and returns
#   list(y0hat   = counterfactual path for the treated-unit mean (length T),
#        tau     = post-period effects, colMeans(treated)[post] - y0hat[post],
#        weights = list(omega, lambda, beta),
#        info    = engine-specific diagnostics)
# Engines never touch data.frames, dates, or fit objects: a Python port maps
# 1:1. The counterfactual *path* (not just tau) is the contract because
# compound and log CARs need paths.

eng_mean <- function(Y, N0, T0) {
  post <- seq.int(T0 + 1L, ncol(Y))
  y0hat <- colMeans(Y[seq_len(N0), , drop = FALSE])
  trt <- colMeans(Y[-seq_len(N0), , drop = FALSE])
  list(y0hat = y0hat, tau = trt[post] - y0hat[post],
       weights = list(omega = NULL, lambda = NULL, beta = NULL), info = list())
}

eng_did <- function(Y, N0, T0) {
  post <- seq.int(T0 + 1L, ncol(Y))
  ctrl <- colMeans(Y[seq_len(N0), , drop = FALSE])
  trt <- colMeans(Y[-seq_len(N0), , drop = FALSE])
  shift <- mean(trt[seq_len(T0)]) - mean(ctrl[seq_len(T0)])
  y0hat <- ctrl + shift
  list(y0hat = y0hat, tau = trt[post] - y0hat[post],
       weights = list(omega = NULL, lambda = NULL, beta = NULL), info = list())
}

# Factor-model abnormal returns. F is a T x K factor matrix aligned to the
# columns of Y. With `beta = NULL`, per-treated-unit OLS of the estimation
# columns on (1, F); with fixed `beta` (length K) the loadings are imposed and
# alpha = 0, which is the market-adjusted convention (beta = 1 on the market
# return). The caller is responsible for the excess-return convention: CAPM
# means passing excess returns / Mkt-RF as the data prescribe.
eng_factor <- function(Y, N0, T0, F, beta = NULL) {
  F <- as.matrix(F)
  post <- seq.int(T0 + 1L, ncol(Y))
  trt <- Y[-seq_len(N0), , drop = FALSE]
  if (!is.null(beta)) {
    fit_units <- matrix(rep(as.vector(F %*% beta), each = nrow(trt)),
                        nrow = nrow(trt))
    coefs <- matrix(c(0, beta), nrow(trt), ncol(F) + 1L, byrow = TRUE)
  } else {
    X <- cbind(1, F[seq_len(T0), , drop = FALSE])
    B <- qr.solve(X, t(trt[, seq_len(T0), drop = FALSE]))  # (K+1) x n1
    fit_units <- t(cbind(1, F) %*% B)                      # n1 x T
    coefs <- t(B)
  }
  dimnames(coefs) <- list(rownames(trt), c("alpha", colnames(F)))
  y0hat <- colMeans(fit_units)
  list(y0hat = y0hat, tau = colMeans(trt)[post] - y0hat[post],
       weights = list(omega = NULL, lambda = NULL, beta = coefs),
       info = list(y0hat_units = fit_units))
}

# Cumulative effects over the event window from per-period mean return paths.
#   sum:      arithmetic sum of effects (cumsum of treated - synthetic)
#   compound: prod(1 + r_treated) - prod(1 + r_synthetic), cumulative
#   log:      cumsum of log(1 + r_treated) - log(1 + r_synthetic)  [Table 6]
car_from_paths <- function(treated, synthetic, cumulate) {
  switch(cumulate,
         sum      = cumsum(treated - synthetic),
         compound = cumprod(1 + treated) - cumprod(1 + synthetic),
         log      = cumsum(log1p(treated) - log1p(synthetic)),
         stop("unknown `cumulate`: ", cumulate))
}
