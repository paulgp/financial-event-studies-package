suppressMessages(library(osqp))
set.seed(42)

# augsynth's synth_qp, verbatim (V = identity), from R/fit_synth.R
synth_qp <- function(X1, X0) {
  V <- diag(ncol(X0))
  Pmat <- X0 %*% V %*% t(X0)
  qvec <- -t(X1) %*% V %*% t(X0)
  n0 <- nrow(X0)
  A <- rbind(rep(1, n0), diag(n0))
  l <- c(1, numeric(n0))
  u <- c(1, rep(1, n0))
  settings <- osqp::osqpSettings(verbose = FALSE, eps_rel = 1e-8, eps_abs = 1e-8)
  sol <- osqp::solve_osqp(P = Pmat, q = qvec, A = A, l = l, u = u, pars = settings)
  sol$x
}

# Frank-Wolfe with exact line search (synthdid solver.R pattern)
# minimize ||A w - b||^2 over the simplex; A is t0 x n0 (columns = donors)
fw_synth <- function(A, b, maxiter = 2000, tol = 1e-12) {
  n0 <- ncol(A)
  w <- rep(1 / n0, n0)
  Aw <- A %*% w
  for (it in 1:maxiter) {
    err <- Aw - b
    grad <- crossprod(A, err)          # O(n0 * t0)
    i <- which.min(grad)
    dA <- A[, i] - Aw
    denom <- sum(dA^2)
    if (denom < 1e-30) break
    step <- max(0, min(1, -sum(err * dA) / denom))
    if (step < tol) break
    w <- (1 - step) * w
    w[i] <- w[i] + step
    Aw <- (1 - step) * Aw + step * A[, i]
  }
  list(w = as.numeric(w), iters = it)
}

obj <- function(A, b, w) sum((A %*% w - b)^2)

run <- function(n0, t0, run_osqp = TRUE) {
  # low-rank factor structure: donors and target share 3 factors + noise
  k <- 3
  F_ <- matrix(rnorm(t0 * k), t0, k)
  L <- matrix(abs(rnorm(n0 * k)), n0, k)
  A <- F_ %*% t(L) + matrix(rnorm(t0 * n0, sd = 0.5), t0, n0)  # t0 x n0
  w_true <- rep(0, n0); w_true[1:5] <- 0.2
  b <- A %*% w_true + rnorm(t0, sd = 0.1)

  t_fw <- system.time(fw <- fw_synth(A, b))[["elapsed"]]
  o_fw <- obj(A, b, fw$w)

  if (run_osqp) {
    t_qp <- system.time(w_qp <- synth_qp(matrix(b, ncol = 1), t(A)))[["elapsed"]]
    o_qp <- obj(A, b, w_qp)
    cat(sprintf("n0=%5d t0=%d | osqp: %7.2fs obj=%.5f | FW: %6.2fs obj=%.5f (%d iters) | speedup %.0fx\n",
                n0, t0, t_qp, o_qp, t_fw, o_fw, fw$iters, t_qp / t_fw))
  } else {
    cat(sprintf("n0=%5d t0=%d | osqp: skipped | FW: %6.2fs obj=%.5f (%d iters)\n",
                n0, t0, t_fw, o_fw, fw$iters))
  }
}

run(500, 100)
run(2000, 100)
run(5000, 100)
