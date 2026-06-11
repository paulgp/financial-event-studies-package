# event_study(): single-event estimation entry point ---------------------------

#' Estimate the effect of a single event on financial returns
#'
#' Fits one of the event-study estimators on a long unit x time return panel
#' and returns a common fit object: the ATT path over the event window,
#' cumulative effects, standard errors, and diagnostics.
#'
#' @param data Long panel: one row per unit-time, with return observations.
#' @param unit,time,ret Column names (strings) in `data`.
#' @param treated Character vector of treated unit ids.
#' @param event_time The value of `time` that is event day 0. Event time is
#'   positional within the panel's sorted unique times (trading periods).
#'   For partial event days (e.g. Geithner's 3pm-close day 0), `data` must
#'   already contain the partial-day return as the day-0 observation.
#' @param method Estimator. `"mean"` (difference in means), `"did"`,
#'   `"market"` (market-adjusted: loading fixed at 1, no intercept, on the
#'   single factor column), `"factor"` (per-unit OLS loadings on the supplied
#'   factor columns; CAPM = one column, FF3F = three), `"sc"`, `"ridge"`,
#'   `"sdid"`, `"gsynth"`.
#' @param window Event window in trading-period offsets, e.g. `c(0, 10)`.
#' @param est_window Pre-event estimation/matching window; must end before
#'   `window` starts. A gap between the two is allowed and those periods are
#'   not loaded (e.g. an excluded placebo window).
#' @param align How window offsets map to panel times: `"position"` (default;
#'   offsets count positions among the panel's sorted unique times — trading
#'   periods) or `"value"` (numeric time values used directly as offsets;
#'   required when the time column is already an event-time index and some
#'   periods are absent from the data, where positional counting would
#'   silently shift the windows).
#' @param returns `"simple"` or `"log"`; required, no default. Records the
#'   return convention of `ret` and selects the default cumulation.
#' @param cumulate How cumulative effects accumulate over the event window:
#'   `"sum"`, `"compound"` (gross-return products differenced), or `"log"`
#'   (sum of log gross-return differences). `"auto"` = sum for log returns,
#'   compound for simple returns.
#' @param factors Data frame with the panel `time` column plus numeric factor
#'   columns, covering every period in both windows. Required for
#'   `"market"`/`"factor"`. Excess-return adjustment (subtracting rf) is the
#'   caller's responsibility.
#' @param donors Optional character vector restricting the donor pool;
#'   default all non-treated units. Units without complete history over the
#'   loaded windows are dropped (recorded in `diagnostics$dropped`).
#' @param match_on For synthetic methods: match on per-period returns
#'   (`"ret"`, the paper's convention) or cumulated pre-event paths
#'   (`"cumret"`).
#' @param V Optional time-weighting matrix (length-T0 vector or T0 x T0
#'   matrix) for the synthetic-control objective; default identity.
#' @param solver Simplex solver for `"sc"`/`"ridge"`: `"hybrid"`
#'   (Frank-Wolfe + support-restricted QP polish, default), `"fw"`, `"qp"`.
#' @param lambda Ridge penalty for `"ridge"`; `NULL` = cross-validated.
#' @param r Factor-number range for `"gsynth"` cross-validation.
#' @param force Fixed effects for `"gsynth"`: `"unit"` (default), `"none"`,
#'   or `"two-way"`.
#' @param se Inference: `"auto"` maps mean/did/market/factor to `"tstat"`,
#'   sc/ridge/sdid to `"placebo"`, gsynth to `"bootstrap"`. `"conformal"`
#'   (methods mean/did/sc/ridge/sdid, `V = NULL`) is the Chernozhukov,
#'   Wuthrich & Zhu (2021) refit-under-the-null procedure with exact
#'   (deterministic, no Monte Carlo) permutation distributions: per-period
#'   p-values and CIs via single-post enumeration, plus a moving-block joint
#'   test and CI for a constant effect over the window. Conformal fits store
#'   CIs and p-values rather than SEs. `"none"` skips.
#' @param reps Placebo repetitions / bootstrap draws (`NULL` = method
#'   default: 100 placebo, 1000 bootstrap).
#' @param level Confidence level for conformal CIs (fixed at fit time;
#'   ignored by other `se` types, which set the level in [confint()]).
#' @param cores Placebo refits run with `parallel::mclapply(mc.cores =
#'   cores)`; all placebo assignments are drawn up front, so results are
#'   identical for any `cores`.
#' @param keep_data Store the panel matrices in the fit (needed for plots,
#'   `event_betas()`, placebo inference re-runs).
#' @param seed Optional RNG seed for placebo/bootstrap inference (conformal
#'   inference is deterministic).
#' @return An object of class `fes_fit`.
#' @export
event_study <- function(data, unit, time, ret, treated, event_time,
                        method = c("mean", "did", "market", "factor",
                                   "sc", "ridge", "sdid", "gsynth"),
                        window = c(0, 10), est_window = c(-250, -11), returns,
                        cumulate = c("auto", "sum", "compound", "log"),
                        align = c("position", "value"),
                        factors = NULL, donors = NULL,
                        match_on = c("ret", "cumret"), V = NULL,
                        solver = c("hybrid", "fw", "qp"), lambda = NULL,
                        r = c(0, 5), force = c("unit", "none", "two-way"),
                        se = c("auto", "placebo", "bootstrap", "tstat",
                               "conformal", "none"),
                        reps = NULL, level = 0.95, cores = 1L,
                        keep_data = TRUE, seed = NULL) {
  method <- match.arg(method)
  cumulate <- match.arg(cumulate)
  match_on <- match.arg(match_on)
  solver <- match.arg(solver)
  se <- match.arg(se)
  returns <- match.arg(returns, c("simple", "log"))
  if (cumulate == "auto") cumulate <- if (returns == "log") "sum" else "compound"
  if (se == "auto")
    se <- switch(method,
                 mean = , did = , market = , factor = "tstat",
                 sc = , ridge = , sdid = "placebo",
                 gsynth = "bootstrap")

  align <- match.arg(align)
  p <- fes_panel(data, unit, time, ret, treated, event_time,
                 window = window, est_window = est_window, donors = donors,
                 align = align)

  F <- if (method %in% c("market", "factor")) align_factors(factors, time, p)
  force <- match.arg(force)
  # matching matrix for synthetic methods (recomputed per placebo draw)
  mk_match <- if (match_on == "cumret") {
    function(Y, N0, T0) {
      A <- apply(Y[seq_len(N0), seq_len(T0), drop = FALSE], 1, cumsum)
      structure(A, b = cumsum(colMeans(Y[-seq_len(N0), seq_len(T0), drop = FALSE])))
    }
  } else function(Y, N0, T0) NULL
  # every closure takes a `w0` warm start; engines without one ignore it
  refit <- switch(method,
    mean   = function(Y, N0, T0, w0 = NULL) eng_mean(Y, N0, T0),
    did    = function(Y, N0, T0, w0 = NULL) eng_did(Y, N0, T0),
    market = {
      if (ncol(F) != 1L) stop("method 'market' needs exactly one factor column")
      function(Y, N0, T0, w0 = NULL) eng_factor(Y, N0, T0, F, beta = 1)
    },
    factor = function(Y, N0, T0, w0 = NULL) eng_factor(Y, N0, T0, F),
    sc     = function(Y, N0, T0, w0 = NULL)
      eng_sc(Y, N0, T0, Ymatch = mk_match(Y, N0, T0), V = V, solver = solver,
             w0 = w0),
    ridge  = function(Y, N0, T0, w0 = NULL)
      eng_ridge_sc(Y, N0, T0, Ymatch = mk_match(Y, N0, T0), V = V,
                   lambda = lambda, solver = solver, w0 = w0),
    sdid   = function(Y, N0, T0, w0 = NULL) eng_sdid(Y, N0, T0),
    gsynth = function(Y, N0, T0, w0 = NULL)
      eng_gsynth(Y, N0, T0, r = r, force = force,
                 se = identical(se, "bootstrap"),
                 nboots = if (is.null(reps)) 1000L else reps))
  eng <- refit(p$Y, p$N0, p$T0)

  post <- seq.int(p$T0 + 1L, ncol(p$Y))
  ev_times <- p$times[post]
  att <- eng$tau
  treated_path <- colMeans(p$Y[-seq_len(p$N0), , drop = FALSE])
  car <- car_from_paths(treated_path[post], eng$y0hat[post], cumulate)
  names(att) <- names(car) <- ev_times

  se_out <- switch(se,
    none = NULL,
    tstat = inf_tstat(p$Y, p$N0, p$T0, eng, method),
    placebo = inf_placebo(p$Y, p$N0, p$T0, n_treated = length(p$treated),
                          refit = refit,
                          reps = if (is.null(reps)) 100L else reps,
                          seed = seed, cores = cores),
    conformal = {
      if (method %in% c("market", "factor", "gsynth"))
        stop("se = 'conformal' is available for methods mean/did/sc/ridge/sdid")
      if (!is.null(V))
        stop("`V` is not supported with conformal inference")
      inf_conformal(p$Y, p$N0, p$T0, refit = refit, att = att, level = level)
    },
    bootstrap = {
      if (method != "gsynth")
        stop("se = 'bootstrap' is only available for method 'gsynth'")
      eng$info$se
    })
  if (!is.null(se_out) && !is.null(se_out$att)) names(se_out$att) <- ev_times
  if (!is.null(se_out) && !is.null(se_out$ci)) {
    dimnames(se_out$ci) <- list(ev_times, c("lower", "upper"))
    names(se_out$p) <- ev_times
  }

  structure(list(
    att = att,
    att_avg = mean(att),
    att_avg_se = se_out$avg,
    car = car,
    se = se_out,
    paths = list(treated = stats::setNames(treated_path, p$times),
                 synthetic = stats::setNames(eng$y0hat, p$times)),
    weights = eng$weights,
    diagnostics = list(n_treated = length(p$treated), n_donors = p$N0,
                       dropped = p$dropped, info = eng$info[setdiff(names(eng$info), "y0hat_units")]),
    panel = if (keep_data) p else NULL,
    method = method,
    conventions = list(returns = returns, cumulate = cumulate,
                       match_on = match_on, window = window,
                       est_window = est_window),
    call = match.call()
  ), class = "fes_fit")
}
