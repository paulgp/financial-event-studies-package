# Generalized synthetic control: wrapper around CRAN gsynth ----------------
#
# Builds a pseudo long panel from the engine matrices, calls gsynth::gsynth
# (Xu 2017: interactive fixed effects, factor number by cross-validation,
# parametric bootstrap), and maps the result back to the engine contract.
# Treatment turns on at the first event-window column, matching the
# replication scripts' `treated = treat * (event_date >= onset)` convention.

eng_gsynth <- function(Y, N0, T0, r = c(0, 5), force = "unit",
                       se = FALSE, nboots = 1000, parallel = FALSE) {
  if (!requireNamespace("gsynth", quietly = TRUE))
    stop("method 'gsynth' requires the gsynth package")
  Tn <- ncol(Y)
  N <- nrow(Y)
  post <- seq.int(T0 + 1L, Tn)
  long <- data.frame(
    id = rep(seq_len(N), times = Tn),
    time = rep(seq_len(Tn), each = N),
    ret = as.vector(Y),
    D = as.numeric(rep(seq_len(N) > N0, times = Tn) &
                     rep(seq_len(Tn) > T0, each = N))
  )
  est <- gsynth::gsynth(ret ~ D, data = long, index = c("id", "time"),
                        force = force, CV = TRUE, r = r, se = se,
                        inference = "parametric", nboots = nboots,
                        parallel = parallel)
  # Y.ct: T x Ntr counterfactuals for the treated units
  yct <- if (is.matrix(est$Y.ct)) rowMeans(est$Y.ct) else as.vector(est$Y.ct)
  trt <- colMeans(Y[-seq_len(N0), , drop = FALSE])
  out_se <- NULL
  if (se) {
    est_att <- est$est.att[post, , drop = FALSE]
    out_se <- list(att = est_att[, "S.E."],
                   avg = unname(est$est.avg[1, "S.E."]),
                   method = "bootstrap", reps = nboots, draws = NULL)
  }
  list(y0hat = yct, tau = trt[post] - yct[post],
       weights = list(omega = NULL, lambda = NULL, beta = NULL),
       info = list(r = est$r.cv, force = force, MSPE = est$MSPE,
                   pre_rmse = sqrt(mean((trt[seq_len(T0)] - yct[seq_len(T0)])^2)),
                   se = out_se))
}
