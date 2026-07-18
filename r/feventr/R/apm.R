# Aggregated projection matrix: wrapper around GitHub apm ------------------
#
# Lei & Ross (arXiv:2312.07520): spectral estimation of counterfactual
# outcome means under a low-rank factor model with general missingness,
# fixed-T / large-N asymptotics. Event-study mapping: donors form the
# fully-observed cohort, treated units observe only the estimation columns,
# and the treated cohort's missing event-window Y(0) means are imputed via
# the aggregated projection matrix. tau = observed treated mean - imputed
# counterfactual mean. SEs are the package's multinomial weighted bootstrap
# over units, conditional on the realized treated path (which carries no
# cross-sectional sampling variation when a single unit is treated).

eng_apm <- function(Y, N0, T0, r = NULL, se = FALSE, nboots = 200L,
                    seed = NULL) {
  if (!requireNamespace("apm", quietly = TRUE))
    stop("method 'apm' requires the apm package: ",
         "remotes::install_github('brad-ross/apm', subdir = 'r')")
  Tn <- ncol(Y)
  pre <- seq_len(T0)
  post <- seq.int(T0 + 1L, Tn)
  if (N0 < 2L) stop("method 'apm' needs at least 2 donors")

  if (length(r) == 1L) {
    r_use <- as.integer(r)
    if (r_use < 1L || r_use > min(T0 - 1L, N0 - 1L))
      stop("method 'apm' needs `r` between 1 and ", min(T0 - 1L, N0 - 1L),
           " for these windows (got ", r_use, ")")
  } else {
    Xc <- Y[seq_len(N0), , drop = FALSE]
    Xc <- Xc - rowMeans(Xc)
    kmax <- min(if (length(r)) max(r) else 8L, min(T0, Tn - T0) - 2L, N0 - 1L)
    if (kmax < 1L)
      stop("method 'apm': no admissible factor count for these windows")
    mu <- pmax(eigen(crossprod(Xc), symmetric = TRUE,
                     only.values = TRUE)$values, 0)
    r_use <- er_factor_count(mu, kmax)
  }

  # apm 0.1.0 calls setorderv/setnames/data.table() unqualified without
  # importing them, so it only works with data.table on the search path;
  # attach it for the duration of the call and restore the search path
  if (!"package:data.table" %in% search()) {
    attachNamespace("data.table")
    on.exit(detach("package:data.table"), add = TRUE)
  }

  units <- rownames(Y)
  if (is.null(units)) units <- paste0("u", seq_len(nrow(Y)))
  oid <- sprintf("t%04d", seq_len(Tn))
  # donors contribute every column, treated units only the estimation
  # columns: post-event treated observations are Y(1), not the Y(0) the
  # factor model is fit on, so they must stay out of the panel
  panel_df <- data.frame(
    unit_id = c(rep(units[seq_len(N0)], times = Tn),
                rep(units[-seq_len(N0)], times = T0)),
    outcome_id = c(rep(oid, each = N0),
                   rep(oid[pre], each = nrow(Y) - N0)),
    y = c(as.vector(Y[seq_len(N0), , drop = FALSE]),
          as.vector(Y[-seq_len(N0), pre, drop = FALSE])))
  panel <- apm::UnbalancedPanel$new(panel_df, "unit_id", "outcome_id", "y",
                                    model_rank = r_use, min_cohort_size = 1,
                                    sort_cohorts_lexicographically = TRUE)
  uc <- as.data.frame(panel$get_unit_cohorts())
  trc <- unique(uc$cohort_id[uc$unit_id %in% units[-seq_len(N0)]])
  if (length(trc) != 1L)
    stop("method 'apm': treated units did not form a single cohort")

  wb <- NULL
  if (se) {
    if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
    wb <- apm::get_weighted_bootstrap_draws(panel$get_num_units(),
                                            as.integer(nboots),
                                            type = "multinomial",
                                            seed = as.integer(seed))
  }
  comps <- apm::est_target_param_components(
    panel,
    est_specs = list(apm = list(factor_model_estimator = "principal_components",
                                include_outcome_fes = FALSE, r = r_use)),
    bootstrap = wb, num_threads = 1L)
  ome <- comps$outcome_means$apm
  y0hat <- as.vector(ome$mean_outcomes()[trc, ])
  trt <- colMeans(Y[-seq_len(N0), , drop = FALSE])
  tau <- trt[post] - y0hat[post]

  out_se <- NULL
  if (se) {
    draws <- t(vapply(seq_len(ome$num_bootstraps()),
                      function(b) trt[post] - ome$mean_outcomes(b)[trc, post],
                      numeric(length(post))))
    out_se <- list(att = apply(draws, 2, stats::sd),
                   avg = stats::sd(rowMeans(draws)),
                   method = "bootstrap", reps = as.integer(nboots),
                   draws = draws)
  }
  list(y0hat = y0hat, tau = tau,
       weights = list(omega = NULL, lambda = NULL, beta = NULL),
       info = list(r = r_use, n_cohorts = panel$get_num_cohorts(),
                   pre_rmse = sqrt(mean((trt[pre] - y0hat[pre])^2)),
                   se = out_se))
}
