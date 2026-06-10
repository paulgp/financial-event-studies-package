#' feventr: Causal Inference for Financial Event Studies
#'
#' Implements the estimators studied in Goldsmith-Pinkham & Lyu, "Causal
#' Inference in Financial Event Studies": difference baselines, factor-model
#' abnormal returns, fast synthetic control, ridge-augmented synthetic
#' control, synthetic difference-in-differences, and generalized synthetic
#' control, with a common panel input and fit object.
#'
#' Provenance: the Frank-Wolfe simplex solver ports the algorithm of the
#' synthdid package (Arkhangelsky, Athey, Hirshberg, Imbens & Wager; dual
#' BSD-3/GPL>=2); ridge augmentation follows augsynth (Ben-Michael, Feller &
#' Rothstein; comparison benchmarks pinned at commit 982f650b).
#'
#' @importFrom stats var cov sd rnorm rbinom plogis qnorm setNames coef vcov confint
#' @importFrom data.table as.data.table dcast setnames
#' @keywords internal
"_PACKAGE"

.datatable.aware <- TRUE

#' Daily Fama-French factors, July 1926 - November 2022
#'
#' Daily Mkt-RF, SMB, HML and RF in decimal units from the Ken French data
#' library (the vintage used by the paper's simulations). Used as the default
#' factor sample for [simulate_events()].
#'
#' @format Data frame with 25,378 rows: `date`, `mktrf`, `smb`, `hml`, `rf`.
#' @source Kenneth R. French data library,
#'   \url{https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html}
"ff_daily"
