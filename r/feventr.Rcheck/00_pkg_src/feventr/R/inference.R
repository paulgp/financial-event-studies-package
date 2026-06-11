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
    ar <- (Y[-seq_len(N0), , drop = FALSE] - eng$info$y0hat_units)[, post, drop = FALSE]
    att_se <- apply(ar, 2L, stats::sd) / sqrt(nrow(ar))
    pooled <- as.vector(ar)
    avg_se <- stats::sd(pooled) / sqrt(length(pooled))
    df <- length(pooled) - 1L
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
inf_placebo <- function(Y, N0, T0, n_treated, refit, reps = 100, seed = NULL) {
  if (N0 <= n_treated)
    stop("placebo inference needs more donors than treated units")
  if (!is.null(seed)) set.seed(seed)
  Tpost <- ncol(Y) - T0
  draws <- matrix(NA_real_, reps, Tpost)
  for (k in seq_len(reps)) {
    fake <- sample.int(N0, n_treated)
    Yp <- rbind(Y[setdiff(seq_len(N0), fake), , drop = FALSE],
                Y[fake, , drop = FALSE])
    draws[k, ] <- refit(Yp, N0 - n_treated, T0)$tau
  }
  avg <- rowMeans(draws)
  list(att = apply(draws, 2L, stats::sd),
       avg = stats::sd(avg),
       df = reps - 1L, method = "placebo", reps = reps, draws = draws)
}
