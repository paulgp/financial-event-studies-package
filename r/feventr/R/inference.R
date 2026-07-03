# Inference ---------------------------------------------------------------
#
# t-statistic inference, matching the conventions of the replication scripts:
# - mean/did: two-sample t between treated and control returns; the average
#   ATT over the event window uses the pooled-variance firm-day t-test
#   (Stata `ttest ret, by(group)` on stacked post-window observations).
# - market/factor: one-sample t of treated units' abnormal returns against 0
#   (the published factor-model columns never touch the control sample).
# Per-period SEs use unequal-variance (Welch) forms for CI paths.

inf_tstat <- function(Y, N0, T0, eng, method) {
  post <- seq.int(T0 + 1L, ncol(Y))
  if (method %in% c("market", "factor")) {
    ar_full <- Y[-seq_len(N0), , drop = FALSE] - eng$info$y0hat_units
    ar <- ar_full[, post, drop = FALSE]
    if (nrow(ar) == 1L) {
      # single treated unit: the cross-sectional sd is NA (one observation per
      # period), so use the classic single-firm event-study SE — the time
      # series sd of the estimation-window abnormal returns (MacKinlay 1997),
      # with df reduced by the number of estimated factor coefficients
      sigma <- stats::sd(ar_full[, seq_len(T0)])
      att_se <- rep(sigma, length(post))
      avg_se <- sigma / sqrt(length(post))
      df <- T0 - ncol(eng$weights$beta)
    } else {
      att_se <- apply(ar, 2L, stats::sd) / sqrt(nrow(ar))
      pooled <- as.vector(ar)
      avg_se <- stats::sd(pooled) / sqrt(length(pooled))
      df <- length(pooled) - 1L
    }
  } else {
    trt <- Y[-seq_len(N0), post, drop = FALSE]
    ctl <- Y[seq_len(N0), post, drop = FALSE]
    v1 <- apply(trt, 2L, stats::var); v0 <- apply(ctl, 2L, stats::var)
    att_se <- sqrt(v1 / nrow(trt) + v0 / nrow(ctl))
    x <- as.vector(trt); y <- as.vector(ctl)
    sp2 <- ((length(x) - 1L) * stats::var(x) + (length(y) - 1L) * stats::var(y)) /
      (length(x) + length(y) - 2L)
    avg_se <- sqrt(sp2 * (1 / length(x) + 1 / length(y)))
    df <- length(x) + length(y) - 2L
  }
  list(att = att_se, avg = avg_se, df = df, method = "tstat",
       reps = NULL, draws = NULL)
}

# Placebo inference (Arkhangelsky et al. 2021 style, as in Stata sdid
# vce(placebo)): repeatedly reassign treatment to a random subset of the
# donors (same number as actually treated), refit the estimator on donors
# only, and use the dispersion of the placebo estimates as the SE.
# All placebo assignments are drawn up front (the refits consume no RNG), so
# the draws are identical for any `cores` and match the previous serial
# implementation seed-for-seed; refits then run in parallel.
# Deliberately NOT warm-started from the full-fit weights: unlike conformal
# refits (same target, one coordinate shifted by h0), each placebo draw
# targets a different pseudo-treated group's mean, so the full-fit solution
# carries no information about the draw's optimum (measured: no speedup,
# sometimes slower via extra KKT-screen rounds).
inf_placebo <- function(Y, N0, T0, n_treated, refit, reps = 100, seed = NULL,
                        cores = 1L) {
  if (N0 <= n_treated)
    stop("placebo inference needs more donors than treated units")
  if (!is.null(seed)) set.seed(seed)
  fakes <- lapply(seq_len(reps), function(k) sample.int(N0, n_treated))
  fit_one <- function(fake) {
    Yp <- rbind(Y[setdiff(seq_len(N0), fake), , drop = FALSE],
                Y[fake, , drop = FALSE])
    refit(Yp, N0 - n_treated, T0)$tau
  }
  par_ok <- cores > 1L && .Platform$OS.type != "windows"
  if (cores > 1L && .Platform$OS.type == "windows")
    warning("parallel placebo refits are not supported on Windows; ",
            "running serially (mc.cores forced to 1)")
  taus <- if (par_ok) parallel::mclapply(fakes, fit_one, mc.cores = cores)
          else lapply(fakes, fit_one)
  # mclapply reports a forked error as a try-error slot (only a warning) and a
  # killed fork as NULL; unchecked, rbind coerces to character (cryptic
  # rowMeans crash) or silently drops rows (SEs over fewer than `reps` draws)
  bad <- vapply(taus, function(t) is.null(t) || inherits(t, "try-error"), TRUE)
  if (any(bad)) {
    first <- taus[bad][[1]]
    detail <- if (inherits(first, "try-error"))
      conditionMessage(attr(first, "condition")) else "fork killed (NULL result)"
    stop(sum(bad), " of ", reps, " placebo refit(s) failed; first: ", detail)
  }
  draws <- do.call(rbind, taus)
  avg <- rowMeans(draws)
  list(att = apply(draws, 2L, stats::sd),
       avg = stats::sd(avg),
       df = reps - 1L, method = "placebo", reps = reps, draws = draws)
}

# Conformal inference (Chernozhukov, Wuthrich & Zhu 2021, JASA): impose the
# null h0, refit with the null-adjusted post period(s) *included in the
# fit*, and test whether the post residual(s) look exchangeable with the
# fit residuals. Two exact shortcuts replace Monte Carlo permutation
# sampling (the fastaugsynth tricks):
#   - pointwise (one post period): a permutation only chooses which single
#     residual lands in the post slot, so the exact permutation
#     distribution of the |statistic| is just |residuals| — enumerate it.
#   - joint constant-effect null over the whole post window: moving-block
#     (cyclic-shift) permutations; all T blocks enumerated.
# CI endpoints come from bracket-expansion + bisection on the p-value (the
# p-value is a step function; granularity 1/(T0+1)). Everything is
# deterministic — no seed, no reps. Each refit reuses the previous
# solution as a Frank-Wolfe warm start: across h0 trials only the linear
# term of the objective moves, so the active donor set barely changes.
inf_conformal <- function(Y, N0, T0, refit, att, level = 0.95) {
  Tn <- ncol(Y)
  T1 <- Tn - T0
  alpha <- 1 - level
  pre <- seq_len(T0)
  wenv <- new.env(parent = emptyenv())
  wenv$w <- NULL

  # residuals over (pre periods + `cols`) from a refit under null h0; the
  # included columns all count as estimation periods, with a duplicated
  # dummy post column to satisfy the engine contract (its tau is unused)
  resids <- function(cols, h0) {
    k <- length(cols)
    Yn <- Y[, c(pre, cols), drop = FALSE]
    Yn[-seq_len(N0), T0 + seq_len(k)] <- Yn[-seq_len(N0), T0 + seq_len(k)] - h0
    Ya <- cbind(Yn, Yn[, T0 + k])
    f <- refit(Ya, N0, T0 + k, w0 = wenv$w)
    om <- f$weights$omega_sc
    if (is.null(om)) om <- f$weights$omega
    if (!is.null(om) && all(om >= -1e-12)) wenv$w <- unname(om)
    e <- colMeans(Ya[-seq_len(N0), , drop = FALSE]) - f$y0hat
    e[seq_len(T0 + k)]
  }

  p_point <- function(j, h0) {
    e <- resids(T0 + j, h0)
    mean(abs(e) >= abs(e[T0 + 1L]) - 1e-12)
  }
  p_joint <- function(h0) {
    e <- abs(resids(seq.int(T0 + 1L, Tn), h0))
    Tt <- length(e)
    block <- seq.int(Tt - T1 + 1L, Tt)
    s <- vapply(seq_len(Tt),
                function(sh) mean(e[((block - 1L + sh) %% Tt) + 1L]), 0)
    mean(s >= s[Tt] - 1e-12)   # shift Tt is the identity block
  }

  # CI = {h0 : p(h0) > alpha}. Bracket-expansion assumes the center is
  # accepted, but the joint constant-effect null h0 = mean(att) can itself be
  # rejected when per-period effects are heterogeneous — then the confidence
  # set is empty and bisection would otherwise fabricate a tight interval
  # around a rejected value. Test the center first and return NA if rejected.
  ci_bounds <- function(pfun, center, scale) {
    if (pfun(center) <= alpha + 1e-12) return(c(NA_real_, NA_real_))
    one <- function(dir) {
      acc <- center
      rej <- NULL
      k <- scale
      for (i in seq_len(60L)) {
        h <- center + dir * k
        if (pfun(h) <= alpha + 1e-12) { rej <- h; break }
        acc <- h
        k <- 2 * k
      }
      if (is.null(rej)) return(dir * Inf)
      while (abs(rej - acc) > scale * 0.01) {
        mid <- (acc + rej) / 2
        if (pfun(mid) <= alpha + 1e-12) rej <- mid else acc <- mid
      }
      (acc + rej) / 2
    }
    c(one(-1), one(1))
  }

  ci <- matrix(NA_real_, T1, 2L)
  p0 <- numeric(T1)
  for (j in seq_len(T1)) {
    scale <- max(stats::sd(resids(T0 + j, att[j])[pre]), 1e-10)
    p0[j] <- p_point(j, 0)
    ci[j, ] <- ci_bounds(function(h) p_point(j, h), att[j], scale)
  }
  avg <- mean(att)
  scale <- max(stats::sd(resids(seq.int(T0 + 1L, Tn), avg)[pre]), 1e-10)
  list(att = NULL, avg = NULL, ci = ci,
       avg_ci = ci_bounds(p_joint, avg, scale),
       p = p0, avg_p = p_joint(0),
       level = level, method = "conformal", reps = NULL, draws = NULL)
}
