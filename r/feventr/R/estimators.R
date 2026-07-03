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

# Synthetic control. Simplex weights on the donors, matched on the
# estimation-window paths: per-period returns (`Ymatch = NULL`, the paper's
# convention) or a caller-supplied t0 x N0 matching matrix (e.g. cumulated
# pre-event paths). Weights are then applied to per-period returns.
eng_sc <- function(Y, N0, T0, Ymatch = NULL, V = NULL, solver = "hybrid",
                   support_size = NULL, max_iter = 2000, w0 = NULL) {
  post <- seq.int(T0 + 1L, ncol(Y))
  A <- if (is.null(Ymatch)) t(Y[seq_len(N0), seq_len(T0), drop = FALSE]) else Ymatch
  trt <- colMeans(Y[-seq_len(N0), , drop = FALSE])
  b <- if (is.null(Ymatch)) trt[seq_len(T0)] else attr(Ymatch, "b")
  sol <- solve_simplex_ls(A, b, V = V, method = solver,
                          support_size = support_size, max_iter = max_iter,
                          w0 = w0)
  omega <- stats::setNames(sol$w, rownames(Y)[seq_len(N0)])
  y0hat <- as.vector(crossprod(Y[seq_len(N0), , drop = FALSE], sol$w))
  list(y0hat = y0hat, tau = trt[post] - y0hat[post],
       weights = list(omega = omega, lambda = NULL, beta = NULL),
       info = list(solver = sol$method, iterations = sol$iterations,
                   objective = sol$objective, n_support = length(sol$support),
                   pre_rmse = sqrt(mean((trt[seq_len(T0)] - y0hat[seq_len(T0)])^2))))
}

# Ridge-augmented synthetic control (augsynth-style). Closed-form
# augmentation on top of the SC weights:
#   w_aug = w_sc + X0 (X0'X0 + lambda I)^{-1} (x1 - X0' w_sc),
# a cheap t0 x t0 solve (Ben-Michael, Feller & Rothstein 2021). With
# `lambda = NULL`, lambda is chosen by leave-one-out CV over estimation
# periods (SC weights held fixed), on a log-spaced grid scaled by the top
# singular value of X0, as in augsynth's ridge_lambda.R.
eng_ridge_sc <- function(Y, N0, T0, Ymatch = NULL, V = NULL, lambda = NULL,
                         solver = "hybrid", support_size = NULL,
                         max_iter = 2000, n_lambda = 20L, w0 = NULL) {
  base <- eng_sc(Y, N0, T0, Ymatch = Ymatch, V = V, solver = solver,
                 support_size = support_size, max_iter = max_iter, w0 = w0)
  w_sc <- unname(base$weights$omega)
  # Augment against the same matched representation the base weights minimize:
  # per-period returns (Ymatch = NULL, the paper's convention) or the caller's
  # matching matrix (e.g. cumulated pre-event paths for match_on = 'cumret').
  # Building X0/x1 from raw returns regardless of Ymatch would correct a
  # different imbalance than the one w_sc was chosen to balance.
  if (is.null(Ymatch)) {
    X0 <- t(Y[seq_len(N0), seq_len(T0), drop = FALSE])   # t0 x n0
    x1 <- colMeans(Y[-seq_len(N0), seq_len(T0), drop = FALSE])
  } else {
    X0 <- Ymatch                                          # t0 x n0 (matched)
    x1 <- attr(Ymatch, "b")
  }
  ridge_w <- function(l, X, y, w) {
    resid <- y - as.vector(X %*% w)
    w + as.vector(crossprod(X, solve(tcrossprod(X) + diag(l, nrow(X)), resid)))
  }
  if (is.null(lambda)) {
    smax <- svd(X0, nu = 0, nv = 0)$d[1]
    grid <- exp(seq(log(smax^2), log(smax^2 * 1e-8), length.out = n_lambda))
    cv_err <- vapply(grid, function(l) {
      sum(vapply(seq_len(T0), function(j) {
        wj <- ridge_w(l, X0[-j, , drop = FALSE], x1[-j], w_sc)
        (x1[j] - sum(X0[j, ] * wj))^2
      }, 0))
    }, 0)
    lambda <- grid[which.min(cv_err)]
  }
  w <- ridge_w(lambda, X0, x1, w_sc)
  trt <- colMeans(Y[-seq_len(N0), , drop = FALSE])
  y0hat <- as.vector(crossprod(Y[seq_len(N0), , drop = FALSE], w))
  post <- seq.int(T0 + 1L, ncol(Y))
  base$info$lambda <- lambda
  base$info$pre_rmse <- sqrt(mean((trt[seq_len(T0)] - y0hat[seq_len(T0)])^2))
  list(y0hat = y0hat, tau = trt[post] - y0hat[post],
       weights = list(omega = stats::setNames(w, rownames(Y)[seq_len(N0)]),
                      omega_sc = base$weights$omega, lambda = NULL, beta = NULL),
       info = base$info)
}

# Synthetic difference-in-differences. Port of synthdid_estimate()'s
# no-covariate path (Arkhangelsky et al. 2021; synthdid package, dual
# BSD-3/GPL>=2): collapsed (N0+1) x (T0+1) form, Frank-Wolfe for time
# weights lambda and unit weights omega with the package's zeta
# regularization and sparsify-then-refit, intercepts via demeaning.
eng_sdid <- function(Y, N0, T0, max_iter = 10000L, sparsify = TRUE) {
  Tn <- ncol(Y)
  N1 <- nrow(Y) - N0
  T1 <- Tn - T0
  pre <- seq_len(T0)
  post <- seq.int(T0 + 1L, Tn)
  Y0 <- Y[seq_len(N0), , drop = FALSE]
  trt <- colMeans(Y[-seq_len(N0), , drop = FALSE])

  noise <- stats::sd(as.vector(apply(Y0[, pre, drop = FALSE], 1, diff)))
  zeta_omega <- ((N1 * T1)^(1 / 4)) * noise
  zeta_lambda <- 1e-6 * noise
  min_dec <- (1e-5 * noise)^2

  fw_reg <- function(A, b, zeta, max_iter, w0 = NULL) {
    # demeaned (intercept) FW with synthdid's objective-decrease stopping:
    # eta = nrow(A) zeta^2; tracked value zeta^2 ||w||^2 + ||err||^2 / nrow(A)
    A <- sweep(A, 2, colMeans(A))
    b <- b - mean(b)
    eta <- nrow(A) * zeta^2
    n <- ncol(A)
    w <- if (is.null(w0)) rep(1 / n, n) else w0
    Aw <- as.vector(A %*% w)
    val_old <- Inf
    for (it in seq_len(max_iter)) {
      st <- fw_step_state(A, b, w, Aw, eta)
      w <- st$w
      Aw <- st$Aw
      val <- zeta^2 * sum(w^2) + sum((Aw - b)^2) / nrow(A)
      if (val_old - val <= min_dec) break
      val_old <- val
    }
    w
  }
  spars <- function(v) { v[v <= max(v) / 4] <- 0; v / sum(v) }

  # time weights lambda: regress post-mean on pre columns across control units
  A_l <- Y0[, pre, drop = FALSE]
  b_l <- rowMeans(Y0[, post, drop = FALSE])
  lam <- fw_reg(A_l, b_l, zeta_lambda, if (sparsify) 100L else max_iter)
  if (sparsify) lam <- fw_reg(A_l, b_l, zeta_lambda, max_iter, w0 = spars(lam))

  # unit weights omega: regress treated pre means on donor pre paths
  A_o <- t(Y0[, pre, drop = FALSE])
  b_o <- trt[pre]
  om <- fw_reg(A_o, b_o, zeta_omega, if (sparsify) 100L else max_iter)
  if (sparsify) om <- fw_reg(A_o, b_o, zeta_omega, max_iter, w0 = spars(om))

  ctrl <- as.vector(crossprod(Y0, om))
  level <- sum(lam * (trt[pre] - ctrl[pre]))
  y0hat <- ctrl + level
  list(y0hat = y0hat, tau = trt[post] - y0hat[post],
       weights = list(omega = stats::setNames(om, rownames(Y)[seq_len(N0)]),
                      lambda = stats::setNames(lam, colnames(Y)[pre]),
                      beta = NULL),
       info = list(zeta_omega = zeta_omega, zeta_lambda = zeta_lambda,
                   noise = noise,
                   pre_rmse = sqrt(mean((trt[pre] - y0hat[pre])^2))))
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
