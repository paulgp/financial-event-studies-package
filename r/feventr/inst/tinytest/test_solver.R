# Solver correctness: solve_simplex_ls vs CVXR reference, across methods,
# V-matrix forms, and ridge penalties; synthetic engines vs independent
# implementations (synthdid, augsynth) on small problems.

set.seed(8)
mk_problem <- function(n0, t0, k = 3) {
  F_ <- matrix(rnorm(t0 * k), t0, k)
  L <- matrix(abs(rnorm(n0 * k)), n0, k)
  A <- F_ %*% t(L) + matrix(rnorm(t0 * n0, sd = 0.5), t0, n0)
  w_true <- rep(0, n0); w_true[1:5] <- 0.2
  list(A = A, b = as.vector(A %*% w_true + rnorm(t0, sd = 0.1)))
}
objective <- function(A, b, w, eta = 0) sum((A %*% w - b)^2) + eta * sum(w^2)

p <- mk_problem(150, 40)

if (requireNamespace("CVXR", quietly = TRUE)) {
  w_ref <- feventr:::ref_simplex_ls(p$A, p$b)
  o_ref <- objective(p$A, p$b, w_ref)
  # 1e-4 relative is the project gate for hybrid/qp; pure FW is allowed its
  # documented sublinear tail (it stalls slightly above the exact optimum)
  for (m in c("hybrid", "fw", "qp")) {
    s <- feventr::solve_simplex_ls(p$A, p$b, method = m, max_iter = 5000)
    tol_m <- if (m == "fw") 5e-2 else 1e-4
    expect_true(abs(s$objective - o_ref) / o_ref < tol_m,
                info = paste("method", m))
    expect_equal(sum(s$w), 1, tolerance = 1e-6)
    expect_true(all(s$w >= -1e-10))
  }
  # exact methods recover the (near-)exact weights, not just the objective
  s_qp <- feventr::solve_simplex_ls(p$A, p$b, method = "qp")
  expect_true(max(abs(s_qp$w - w_ref)) < 1e-3)

  # V as a diagonal vector and as a matrix
  v <- runif(40, 0.5, 2)
  w_ref_v <- feventr:::ref_simplex_ls(p$A, p$b, V = v)
  s_v <- feventr::solve_simplex_ls(p$A, p$b, V = v, method = "hybrid")
  expect_true(abs(objective(sqrt(v) * p$A, sqrt(v) * p$b, s_v$w) -
                    objective(sqrt(v) * p$A, sqrt(v) * p$b, w_ref_v)) /
                objective(sqrt(v) * p$A, sqrt(v) * p$b, w_ref_v) < 1e-4)
  s_vm <- feventr::solve_simplex_ls(p$A, p$b, V = diag(v), method = "qp")
  expect_true(max(abs(s_vm$w - s_v$w)) < 1e-3)

  # ridge penalty (eta > 0)
  w_ref_e <- feventr:::ref_simplex_ls(p$A, p$b, eta = 5)
  s_e <- feventr::solve_simplex_ls(p$A, p$b, eta = 5, method = "qp")
  expect_true(max(abs(s_e$w - w_ref_e)) < 1e-3)
} else {
  exit_file("CVXR not available")
}

# hybrid support restriction binds and stays near-exact
p2 <- mk_problem(400, 30)
s_full <- feventr::solve_simplex_ls(p2$A, p2$b, method = "qp")
s_hyb <- feventr::solve_simplex_ls(p2$A, p2$b, method = "hybrid",
                                   support_size = 60)
expect_true(length(s_hyb$support) <= 60)
expect_true((s_hyb$objective - s_full$objective) /
              s_full$objective < 1e-4)

# --- eng_sc / eng_sdid vs independent packages -------------------------------
set.seed(15)
sim <- feventr::simulate_events(n_units = 80, n_pre = 60, n_candidate = 1,
                                n_post = 8, tau = 0.02, seed = 31)
args <- list(data = sim$data, unit = "id", time = "t", ret = "ret",
             treated = sim$events$unit, event_time = sim$event_time,
             window = c(0, 8), est_window = c(-60, -1), returns = "simple",
             se = "none")

fit_sc <- do.call(feventr::event_study, c(args, list(method = "sc")))
expect_equal(sum(fit_sc$weights$omega), 1, tolerance = 1e-6)
expect_true(fit_sc$diagnostics$info$pre_rmse < 0.01)
# SC weights solve the same problem as the reference solver
if (requireNamespace("CVXR", quietly = TRUE)) {
  pn <- fit_sc$panel
  A <- t(pn$Y[seq_len(pn$N0), seq_len(pn$T0)])
  b <- colMeans(pn$Y[-seq_len(pn$N0), seq_len(pn$T0), drop = FALSE])
  w_ref <- feventr:::ref_simplex_ls(A, b)
  expect_true((objective(A, b, fit_sc$weights$omega) - objective(A, b, w_ref)) /
                objective(A, b, w_ref) < 1e-4)
}

# augsynth agreement: same SCM objective, so ATT paths should align closely
if (requireNamespace("augsynth", quietly = TRUE)) {
  dd <- sim$data
  dd$trt <- as.integer(dd$id %in% as.integer(sim$events$unit) &
                         dd$t >= sim$event_time)
  dd <- dd[dd$t >= sim$event_time - 60 & dd$t <= sim$event_time + 8, ]
  aug <- suppressWarnings(suppressMessages(
    augsynth::augsynth(ret ~ trt, unit = id, time = t, data = dd,
                       progfunc = "None", scm = TRUE)))
  att_aug <- predict(aug, att = TRUE)[as.character(sim$event_time:(sim$event_time + 8))]
  expect_equal(unname(fit_sc$att), unname(att_aug), tolerance = 1e-3)
}

# synthdid agreement: our eng_sdid ports its algorithm
if (requireNamespace("synthdid", quietly = TRUE)) {
  fit_sdid <- do.call(feventr::event_study, c(args, list(method = "sdid")))
  pn <- fit_sdid$panel
  est <- synthdid::synthdid_estimate(pn$Y, pn$N0, pn$T0)
  expect_equal(fit_sdid$att_avg, unname(c(est)), tolerance = 1e-6)
  w <- attr(est, "weights")
  expect_equal(unname(fit_sdid$weights$omega), unname(w$omega), tolerance = 1e-4)
  expect_equal(unname(fit_sdid$weights$lambda), unname(w$lambda), tolerance = 1e-4)
}

# ridge augmentation: fixed lambda matches the closed form; CV runs
fit_r <- do.call(feventr::event_study, c(args, list(method = "ridge", lambda = 1)))
pn <- fit_r$panel
X0 <- t(pn$Y[seq_len(pn$N0), seq_len(pn$T0)])
x1 <- colMeans(pn$Y[-seq_len(pn$N0), seq_len(pn$T0), drop = FALSE])
w_sc <- unname(fit_r$weights$omega_sc)
w_man <- w_sc + as.vector(crossprod(X0, solve(tcrossprod(X0) + diag(1, pn$T0),
                                              x1 - as.vector(X0 %*% w_sc))))
expect_equal(unname(fit_r$weights$omega), w_man, tolerance = 1e-8)
fit_rcv <- do.call(feventr::event_study, c(args, list(method = "ridge")))
expect_true(is.finite(fit_rcv$diagnostics$info$lambda))
# augmentation can only improve pre-period fit
expect_true(fit_rcv$diagnostics$info$pre_rmse <= fit_sc$diagnostics$info$pre_rmse + 1e-12)

# the closed-form leave-one-out CV selects the same lambda as the naive
# n_lambda x T0 double loop it replaced (issue 16)
pnr <- fit_rcv$panel
X0r <- t(pnr$Y[seq_len(pnr$N0), seq_len(pnr$T0), drop = FALSE])
x1r <- colMeans(pnr$Y[-seq_len(pnr$N0), seq_len(pnr$T0), drop = FALSE])
wsc_r <- unname(fit_rcv$weights$omega_sc)
ridge_w <- function(l, X, y, w) {
  rr <- y - as.vector(X %*% w)
  w + as.vector(crossprod(X, solve(tcrossprod(X) + diag(l, nrow(X)), rr)))
}
smax_r <- svd(X0r, nu = 0, nv = 0)$d[1]
grid_r <- exp(seq(log(smax_r^2), log(smax_r^2 * 1e-8), length.out = 20))
cv_naive <- vapply(grid_r, function(l) sum(vapply(seq_len(pnr$T0), function(j) {
  wj <- ridge_w(l, X0r[-j, , drop = FALSE], x1r[-j], wsc_r)
  (x1r[j] - sum(X0r[j, ] * wj))^2
}, 0)), 0)
expect_equal(fit_rcv$diagnostics$info$lambda, grid_r[which.min(cv_naive)],
             tolerance = 1e-10)

# ridge conformal/placebo (which reuse the frozen lambda instead of re-CVing
# on every refit) still produce finite inference
fit_rc_conf <- do.call(feventr::event_study,
                       c(args[setdiff(names(args), "se")],
                         list(method = "ridge", se = "conformal")))
expect_true(all(is.finite(fit_rc_conf$se$ci)))

# match_on = 'cumret': the ridge augmentation must correct the *cumret*
# imbalance the base weights minimize, not the raw-return imbalance (issue 3).
# The augmented weights equal the closed form built from the cumret matching
# matrix, not from per-period returns.
fit_rc <- do.call(feventr::event_study,
                  c(args, list(method = "ridge", lambda = 1, match_on = "cumret")))
pnc <- fit_rc$panel
pre_c <- seq_len(pnc$T0)
Xc <- apply(pnc$Y[seq_len(pnc$N0), pre_c, drop = FALSE], 1, cumsum)   # t0 x n0
xc1 <- cumsum(colMeans(pnc$Y[-seq_len(pnc$N0), pre_c, drop = FALSE]))
wsc_c <- unname(fit_rc$weights$omega_sc)
w_man_c <- wsc_c + as.vector(crossprod(Xc, solve(tcrossprod(Xc) + diag(1, pnc$T0),
                                                 xc1 - as.vector(Xc %*% wsc_c))))
expect_equal(unname(fit_rc$weights$omega), w_man_c, tolerance = 1e-8)
# and it differs from the (buggy) raw-return augmentation
X0c <- t(pnc$Y[seq_len(pnc$N0), pre_c, drop = FALSE])
x1c <- colMeans(pnc$Y[-seq_len(pnc$N0), pre_c, drop = FALSE])
w_raw <- wsc_c + as.vector(crossprod(X0c, solve(tcrossprod(X0c) + diag(1, pnc$T0),
                                                x1c - as.vector(X0c %*% wsc_c))))
expect_true(max(abs(w_man_c - w_raw)) > 1e-6)

# placebo inference populates draws and SEs
fit_pl <- do.call(feventr::event_study,
                  c(args[setdiff(names(args), "se")],
                    list(method = "sc", se = "placebo", reps = 20, seed = 1)))
expect_equal(fit_pl$se$method, "placebo")
expect_equal(dim(fit_pl$se$draws), c(20L, 9L))
expect_true(all(is.finite(fit_pl$se$att)))
expect_true(is.finite(fit_pl$att_avg_se))

# gsynth wrapper maps back to the engine contract
if (requireNamespace("gsynth", quietly = TRUE)) {
  fit_g <- do.call(feventr::event_study,
                   c(args, list(method = "gsynth", r = c(0, 3), force = "none")))
  expect_equal(length(fit_g$att), 9L)
  expect_true(fit_g$diagnostics$info$pre_rmse < 0.02)
  expect_true(abs(fit_g$att[["0"]] - sim$tau) < 0.02)
}

# --- event_betas diagnostics --------------------------------------------------
bet <- feventr::event_betas(fit_sc, sim$factors[, c("t", "mktrf", "smb")])
expect_inherits(bet, "fes_betas")
expect_true(all(c("treated", "control", "control_weighted_sc") %in% bet$group))
# treated mean loadings recover the DGP betas (N(1, 0.3^2))
tm <- bet[bet$group == "treated" & bet$stat == "mean", ]
expect_true(abs(tm$mktrf - mean(sim$betas$b_mkt[sim$betas$treated])) < 0.15)
# weighted control loadings equal the manual weighted mean
pn <- fit_sc$panel
pre <- seq_len(pn$T0)
X <- cbind(1, as.matrix(sim$factors[match(pn$time_values, sim$factors$t),
                                    c("mktrf", "smb")])[pre, ])
Bman <- t(qr.solve(X, t(pn$Y[seq_len(pn$N0), pre])))
wm <- bet[bet$group == "control_weighted_sc" & bet$stat == "mean", "mktrf"]
expect_equal(wm, unname(weighted.mean(Bman[, 2], fit_sc$weights$omega)),
             tolerance = 1e-10)
