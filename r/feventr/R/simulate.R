# Simulation DGP ----------------------------------------------------------
#
# Two-factor DGP with selection-on-assignment and selection-on-timing from
# Goldsmith-Pinkham & Lyu, Section 4.1 (Table 1). Mirrors the published
# sim_data() in PEAD_DinD/code/simulations_selection_TL.R exactly, including
# its RNG consumption order (residuals -> beta_mkt -> beta_smb -> factor-day
# sample -> [dead timing-logit draw] -> assignment draw), so that
# `simulate_events(seed = s)` reproduces the paper's simulated panels
# bit-for-bit for seeds 1234..1283.

#' Simulate an event-study return panel with a two-factor structure
#'
#' Returns r_it = beta_mkt_i * MktRF_t + beta_smb_i * SMB_t + rf_t + eps_it,
#' with a one-shot treatment effect `tau` added on the event day for treated
#' units. Factor days are sampled jointly (whole rows, preserving the
#' contemporaneous MktRF-SMB correlation) without replacement from `factors`.
#' The residual sd defaults to the paper's rule that idiosyncratic variance
#' equals systematic variance: sqrt(var(MktRF) + var(SMB) + 2 cov).
#'
#' Selection mechanisms (paper Section 4.1): `"assignment"` treats units by a
#' logit in beta_smb scaled so that on average `treat_share` are treated
#' (lower loadings more likely treated); `"timing"` places the event on the
#' day with the largest SMB realization among the `n_candidate` candidate
#' days (the paper's logit-timing draw is vestigial dead code in the
#' published DGP, but it is replayed here because it consumes RNG).
#'
#' @param n_units Number of firms.
#' @param n_pre Estimation periods before the earliest candidate event day.
#' @param n_candidate Candidate event days (event day 0 is the first when
#'   timing is random).
#' @param n_post Post-event periods simulated after the candidate window.
#' @param tau Treatment effect added on the event day.
#' @param treat_share Average treated share.
#' @param beta_mean,beta_sd Length-2 (mkt, smb) loading distribution params.
#' @param resid_sd Idiosyncratic sd; `NULL` = paper's 50%-variance rule.
#' @param selection One of `"none"`, `"assignment"`, `"timing"`, `"both"`.
#' @param factors Data frame with columns `mktrf`, `smb`, `rf` (decimal daily
#'   returns); default the bundled `ff_daily` (Ken French daily factors,
#'   1926-07-01..2022-11-30, the paper's sample).
#' @param seed Optional RNG seed.
#' @return A list of class `fes_sim`: `data` (long panel: id, t, ret),
#'   `events` (treated unit ids + event_time, ready for `event_study()`),
#'   `betas` (per-unit loadings and treatment flag), `factors` (per-period
#'   realized mktrf, smb, rf), `event_time`, `tau`.
#' @export
simulate_events <- function(n_units = 500, n_pre = 239, n_candidate = 250,
                            n_post = 10, tau = 0.03, treat_share = 0.1,
                            beta_mean = c(mkt = 1, smb = 1),
                            beta_sd = c(mkt = 0.3, smb = 0.3),
                            resid_sd = NULL,
                            selection = c("none", "assignment", "timing", "both"),
                            factors = NULL, seed = NULL) {
  selection <- match.arg(selection)
  sel_assign <- selection %in% c("assignment", "both")
  sel_timing <- selection %in% c("timing", "both")
  if (is.null(factors)) {
    # robust lazy-data access whether or not the package is attached
    e <- new.env()
    utils::data("ff_daily", package = "feventr", envir = e)
    factors <- e$ff_daily
  }
  stopifnot(all(c("mktrf", "smb", "rf") %in% names(factors)))
  if (is.null(resid_sd))
    resid_sd <- sqrt(stats::var(factors$mktrf) + stats::var(factors$smb) +
                       2 * stats::cov(factors$mktrf, factors$smb))
  if (!is.null(seed)) set.seed(seed)

  P <- n_pre + 1L                       # event day position when timing is random
  T <- P + n_candidate + n_post
  n <- n_units
  # RNG order below mirrors the published sim_data() exactly.
  e_it <- matrix(stats::rnorm(n * T, 0, resid_sd), n, T)
  beta <- cbind(stats::rnorm(n, beta_mean[1], beta_sd[1]),
                stats::rnorm(n, beta_mean[2], beta_sd[2]))
  ran_ind <- sample(nrow(factors), T, replace = FALSE)
  F <- as.matrix(factors[ran_ind, c("mktrf", "smb")])
  rf <- factors$rf[ran_ind]
  r_it <- beta %*% t(F) + e_it + rep(1, n) %*% t(rf)

  if (sel_timing) {
    sub <- rank(F[P:(P + n_candidate - 1L), "smb"])
    sub <- max(sub) - sub + 1
    scale_factor <- log(1 / n_candidate) / mean(sub)
    selected <- stats::rbinom(n_candidate, 1, stats::plogis(scale_factor * sub)) == 1
    # which.min keeps this scalar: rounded SMB values tie (only ~600 distinct
    # values), and which(sub == min(sub)) would return every tied day, making
    # event_time a length > 1 vector that corrupts the events table
    event_time <- P - 1L + unname(which.min(sub))  # the max-SMB candidate day
  } else {
    event_time <- P
  }

  if (sel_assign) {
    scale_factor <- log(treat_share) / mean(beta[, 2])
    treated <- stats::rbinom(n, 1, stats::plogis(scale_factor * beta[, 2])) == 1
  } else {
    treated <- stats::rbinom(n, 1, treat_share) == 1
  }
  r_it[treated, event_time] <- r_it[treated, event_time] + tau

  structure(list(
    data = data.frame(id = rep(seq_len(n), times = T),
                      t = rep(seq_len(T), each = n),
                      ret = as.vector(r_it)),
    events = data.frame(unit = as.character(which(treated)),
                        event_time = event_time),
    betas = data.frame(id = seq_len(n), b_mkt = beta[, 1], b_smb = beta[, 2],
                       treated = treated),
    factors = data.frame(t = seq_len(T), mktrf = F[, "mktrf"],
                         smb = F[, "smb"], rf = rf),
    event_time = event_time,
    tau = tau
  ), class = "fes_sim")
}
