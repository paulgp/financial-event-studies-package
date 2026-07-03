# Cross-sectional event study with pre-event benchmarking (Cohn, Johnson, Liu
# & Wardlaw 2026, JFE) --------------------------------------------------------
#
# A different axis from the time-series estimators in this package (mean/did/
# market/factor/sc/ridge/sdid/gsynth). Those estimate a counterfactual return
# path for treated unit(s) and report per-period abnormal returns. Here the
# estimand is a single cross-sectional slope: on the event day, regress every
# firm's short-horizon return on a firm characteristic, and ask whether that
# slope is abnormal relative to the same regression run on pre-event days.
#
# The threat this addresses is confounding contemporaneous news: firms with
# similar characteristics comove, so the regression's own (even industry-
# clustered) standard error tests "is there a relationship today" rather than
# "did the event cause one", and over-rejects on ordinary days. The counter-
# factual is the empirical distribution of the slope on pre-event days.
#
# Numerics (PCA covariance, Cholesky/Woodbury GLS, empirical-CDF + parametric
# statistics) port the authors' Stata `csestudy` / Mata reference and its R
# and Python translations. Kept self-contained (like calendar_time): the GLS
# covariance is a *rolling* pre-event window that shifts with each pseudo-event
# date, and validity is assessed per date, neither of which fits fes_panel's
# single balanced treated-vs-donor matrix.

# PCA covariance decomposition of a T_pre x N pre-event return matrix `R`
# (rows = pre-event days, columns = firms). Returns loadings V (N x npc),
# factor variances lam (npc), idiosyncratic variances d (N), so that
#   Omega = V diag(lam) V' + diag(d).
# Demean per firm, economy SVD, first npc right singular vectors, ddof = 1
# variances, floor d at 1e-15 (matching the reference).
.cse_pca <- function(R, npc) {
  A <- sweep(R, 2L, colMeans(R), "-")
  sv <- svd(A)                                   # sv$v is N x min(T_pre, N)
  V <- sv$v[, seq_len(npc), drop = FALSE]
  scores <- A %*% V
  lam <- apply(scores, 2L, stats::var)           # ddof = 1
  resid <- A - tcrossprod(scores, V)             # A - scores %*% t(V)
  d <- pmax(apply(resid, 2L, stats::var), 1e-15) # ddof = 1, per firm
  list(V = V, lam = lam, d = d)
}

# OLS slope(s) of y on X (intercept is the LAST column of X, added by caller).
.cse_ols <- function(y, X) unname(qr.solve(X, y))

# GLS via Cholesky of Omega = V diag(lam) V' + diag(d): whiten y, X by L^-1
# (Omega = L L') then OLS on the whitened data. Default; most precise.
.cse_gls_chol <- function(y, X, V, lam, d) {
  Omega <- tcrossprod(V %*% diag(lam, nrow = length(lam)), V)  # V diag(lam) V'
  diag(Omega) <- diag(Omega) + d
  Omega <- (Omega + t(Omega)) / 2                # force symmetry (numerical)
  U <- chol(Omega)                               # upper: Omega = t(U) %*% U
  yw <- backsolve(U, y, transpose = TRUE)        # L^-1 y, L = t(U)
  Xw <- backsolve(U, X, transpose = TRUE)
  unname(qr.solve(Xw, yw))
}

# GLS via the Woodbury identity: only an npc x npc inverse, no N x N factor.
#   Omega^-1 = D^-1 - D^-1 V M V' D^-1,  M = (Lambda^-1 + V' D^-1 V)^-1
# Faster for large cross-sections, slightly less precise.
.cse_gls_woodbury <- function(y, X, V, lam, d) {
  d_inv <- 1 / d
  VtDinv <- sweep(t(V), 2L, d_inv, "*")          # (npc x N): V' D^-1
  M <- solve(diag(1 / pmax(lam, 1e-15), nrow = length(lam)) + VtDinv %*% V)
  oinv <- function(B) (d_inv * B) - (d_inv * V) %*% (M %*% (VtDinv %*% B))
  unname(solve(crossprod(X, oinv(X)), crossprod(X, oinv(y))))
}

# Significance from the (L+1) x k coefficient matrix whose first row is the
# event-date estimate and remaining L rows the pre-event estimates.
#   - empirical CDF p-value: fraction of dates whose |coef - pre_mean| is at
#     least the event's, weak inequality over all L+1 dates (denominator L+1);
#   - parametric two-tailed p-value with an (L+1)/L variance adjustment on a
#     t(L-1) reference (matching the reference implementation).
.cse_significance <- function(all_betas, beta_event) {
  pre <- all_betas[-1L, , drop = FALSE]
  pre_mean <- colMeans(pre, na.rm = TRUE)
  pre_sd <- sqrt(apply(pre, 2L, stats::var, na.rm = TRUE))
  Lp1 <- nrow(all_betas)
  dev <- abs(sweep(all_betas, 2L, pre_mean, "-"))
  ge <- sweep(dev, 2L, abs(beta_event - pre_mean), FUN = ">=")
  p_cdf <- colSums(ge, na.rm = TRUE) / Lp1
  z <- abs(beta_event - pre_mean) / (pre_sd * sqrt(Lp1 / (Lp1 - 1)))
  p_parametric <- 2 * stats::pt(z, df = Lp1 - 2, lower.tail = FALSE)
  list(p_cdf = p_cdf, p_parametric = p_parametric,
       pre_mean = pre_mean, pre_sd = pre_sd)
}

#' Cross-sectional event study with pre-event benchmarking
#'
#' Estimates the event-day cross-sectional relationship between short-horizon
#' returns and firm characteristics, and benchmarks it against the distribution
#' of the same relationship on a window of pre-event ("pseudo-event") days.
#' This is the Cohn, Johnson, Liu & Wardlaw (2026) estimator, a different axis
#' from the treated-versus-counterfactual estimators in [event_study()]: the
#' unit of analysis is a single cross-sectional slope on the event day, not a
#' per-unit abnormal-return path.
#'
#' On each date the model regresses `ret` on the characteristics `chars` (plus
#' an intercept) across all firms present that day:
#' \deqn{r_{i,t} = a_t + b_t' x_{i,t} + e_{i,t}.}
#' Confounding contemporaneous news makes the errors `e` cross-sectionally
#' correlated among firms with similar `x`, so the regression's own standard
#' error tests the wrong null (\eqn{b_t = 0}, "any relationship today") instead
#' of the null of interest ("the *event* caused the relationship"), and
#' over-rejects on ordinary days. Instead of trusting that standard error, the
#' event-day slope \eqn{b_0} is compared to the empirical distribution of
#' \eqn{b_t} over the pre-event window: the counterfactual is the pre-event
#' slopes, not a fitted return path.
#'
#' With `method = "gls"`, each date's regression is generalized least squares
#' using a return covariance estimated from the leading principal components of
#' that date's own preceding-return window (\eqn{\Omega = V\,\mathrm{diag}
#' (\lambda)\,V' + \mathrm{diag}(d)}), which down-weights common-factor
#' comovement and sharply increases power over OLS. The GLS covariance window
#' *rolls* with each pseudo-event date, so the panel must extend
#' `-pre_window[1]` further back than `pre_window[1]` itself.
#'
#' @param data Long panel: one row per unit-time.
#' @param unit,time,ret Column names (strings) in `data`. `ret` is the
#'   short-horizon (e.g. one-day) event return that is regressed on `chars`;
#'   for `method = "gls"` it is also the series the covariance is built from,
#'   so it should be a per-period (not cumulated) return.
#' @param event_time The value of `time` that is the event day. Time is
#'   positional within the panel's sorted unique times (trading periods), as in
#'   [event_study()].
#' @param chars Character vector of characteristic column names (the regressors,
#'   observed on each date; lag them yourself if desired). `NULL` fits an
#'   intercept-only model (mean return per date).
#' @param pre_window Pre-event benchmark window in trading-period offsets,
#'   `c(first, last)` with `first < last < 0`, e.g. `c(-200, -1)`: the pseudo-
#'   event days from `first` to `last` positions before the event.
#' @param method `"ols"` (default) or `"gls"`.
#' @param npc Principal components for the GLS covariance (default 100); must
#'   not exceed the pre-window length.
#' @param solver GLS solver: `"cholesky"` (default, most precise) or
#'   `"woodbury"` (faster for large cross-sections, slightly less precise).
#' @param min_return For GLS, firms whose absolute returns over the covariance
#'   window sum to less than this are excluded (drops (near-)constant histories
#'   that would make the covariance singular); default `0.01`.
#' @param verbose Print progress over the pre-event loop (the slow GLS path).
#' @return An object of class `fes_cse`: `$coefficients` (term, estimate,
#'   `p_cdf`, `p_parametric`, pre-event mean/sd), `$params` (named event-day
#'   coefficients), `$p_cdf`, `$p_parametric`, `$pre_betas` (L x k pseudo-event
#'   coefficients), `$n_obs` (firms per date, event first), `$diagnostics`,
#'   `$conventions`.
#' @references Cohn, J. B., Johnson, T. L., Liu, Z. & Wardlaw, M. I. (2026).
#'   Past is Prologue: Inference from the Cross Section of Returns Around an
#'   Event. *Journal of Financial Economics*.
#' @export
cross_section <- function(data, unit, time, ret, event_time, chars = NULL,
                          pre_window = c(-250, -1),
                          method = c("ols", "gls"), npc = 100L,
                          solver = c("cholesky", "woodbury"),
                          min_return = 0.01, verbose = FALSE) {
  method <- match.arg(method)
  solver <- match.arg(solver)
  gls <- method == "gls"
  if (solver == "woodbury" && !gls)
    stop("solver = 'woodbury' requires method = 'gls'")
  if (length(pre_window) != 2L || pre_window[1] >= pre_window[2] ||
      pre_window[2] >= 0)
    stop("`pre_window` must be c(first, last) with first < last < 0")

  data <- as.data.frame(data)
  needed <- c(unit, time, ret, chars)
  miss <- setdiff(needed, names(data))
  if (length(miss)) stop("columns not found in `data`: ", paste(miss, collapse = ", "))
  data <- data[needed]
  data <- data[stats::complete.cases(data), , drop = FALSE]
  if (!nrow(data)) stop("no complete observations in `data`")

  panel_id <- data[[unit]]
  time_id <- data[[time]]
  y_all <- as.numeric(data[[ret]])
  X_rhs <- if (length(chars)) {
    m <- as.matrix(data[chars]); storage.mode(m) <- "double"; m
  } else matrix(numeric(0), nrow(data), 0L)

  # rectangular (unit x time) lookups: row position and return per cell
  panels <- sort(unique(panel_id))
  times <- sort(unique(time_id))
  pcode <- match(panel_id, panels)
  tcode <- match(time_id, times)
  rect_row <- matrix(NA_integer_, length(panels), length(times))
  rect_row[cbind(pcode, tcode)] <- seq_along(panel_id)
  rect_y <- matrix(NA_real_, length(panels), length(times))
  rect_y[cbind(pcode, tcode)] <- y_all

  event_col <- match(event_time, times)
  if (is.na(event_col)) stop("`event_time` not found among panel times")
  pre_start_col <- event_col + pre_window[1]   # earliest pseudo-event day
  pre_end_col <- event_col + pre_window[2]     # latest pseudo-event day
  if (pre_start_col < 1L || event_col > length(times))
    stop("`pre_window` extends beyond the panel's time range")
  n_pre <- pre_end_col - pre_start_col + 1L
  if (gls && npc > n_pre)
    stop("`npc` (", npc, ") must be <= pre-window length (", n_pre, ")")
  # each date's covariance window sits at the same offsets relative to it, so
  # the earliest pseudo-event day needs data another -pre_window[1] periods back
  if (gls && (pre_start_col + pre_window[1]) < 1L)
    stop("`method = 'gls'` needs the panel to extend ", -pre_window[1],
         " periods before the first pseudo-event day; it does not")

  param_names <- c(chars, "(Intercept)")
  n_params <- length(param_names)

  run_date <- function(target_col) {
    valid <- !is.na(rect_row[, target_col])
    win <- NULL
    if (gls) {
      win <- (target_col + pre_window[1]):(target_col + pre_window[2])
      present <- !is.na(rect_row[, c(win, target_col), drop = FALSE])
      valid <- valid & (rowSums(present) == length(win) + 1L)
      valid <- valid &
        (rowSums(abs(rect_y[, win, drop = FALSE]), na.rm = TRUE) >= min_return)
    }
    mask <- which(valid)
    if (!length(mask)) return(list(beta = rep(NA_real_, n_params), n = 0L))
    if (gls && length(mask) < npc)
      stop("date offset ", target_col - event_col, " has ", length(mask),
           " valid firms, fewer than npc = ", npc)
    rows <- rect_row[mask, target_col]
    y <- y_all[rows]
    X <- if (ncol(X_rhs)) cbind(X_rhs[rows, , drop = FALSE], 1) else
      matrix(1, length(mask), 1L)
    beta <- if (gls) {
      dec <- .cse_pca(t(rect_y[mask, win, drop = FALSE]), npc)  # T_pre x N
      if (solver == "woodbury") .cse_gls_woodbury(y, X, dec$V, dec$lam, dec$d)
      else .cse_gls_chol(y, X, dec$V, dec$lam, dec$d)
    } else .cse_ols(y, X)
    list(beta = as.numeric(beta), n = length(mask))
  }

  ev <- run_date(event_col)
  if (!ev$n) stop("no valid firms on the event date")

  all_betas <- matrix(NA_real_, n_pre + 1L, n_params,
                      dimnames = list(NULL, param_names))
  all_n <- integer(n_pre + 1L)
  all_betas[1L, ] <- ev$beta
  all_n[1L] <- ev$n
  pre_cols <- pre_end_col:pre_start_col          # latest -> earliest
  for (i in seq_along(pre_cols)) {
    rd <- run_date(pre_cols[i])
    all_betas[i + 1L, ] <- rd$beta
    all_n[i + 1L] <- rd$n
    if (verbose && (i %% 25L == 0L || i == n_pre))
      message("  pseudo-event day ", i, " / ", n_pre)
  }

  sig <- .cse_significance(all_betas, ev$beta)
  method_label <- if (gls)
    if (solver == "woodbury") "gls (woodbury)" else "gls" else "ols"

  structure(list(
    coefficients = data.frame(
      term = param_names, estimate = ev$beta,
      p_cdf = unname(sig$p_cdf), p_parametric = unname(sig$p_parametric),
      pre_mean = unname(sig$pre_mean), pre_sd = unname(sig$pre_sd),
      row.names = NULL),
    params = stats::setNames(ev$beta, param_names),
    p_cdf = stats::setNames(sig$p_cdf, param_names),
    p_parametric = stats::setNames(sig$p_parametric, param_names),
    pre_betas = all_betas[-1L, , drop = FALSE],
    n_obs = all_n,
    diagnostics = list(n_obs_event = ev$n, n_pre_days = n_pre,
                       mean_units = mean(all_n)),
    method = method_label,
    conventions = list(event_time = event_time, pre_window = pre_window,
                       npc = if (gls) npc, solver = if (gls) solver,
                       min_return = if (gls) min_return),
    call = match.call()
  ), class = "fes_cse")
}

#' @export
print.fes_cse <- function(x, ...) {
  w <- x$conventions$pre_window
  cat("feventr cross-sectional event study: method '", x$method, "'\n", sep = "")
  cat(x$diagnostics$n_obs_event, " firms on the event day, ",
      x$diagnostics$n_pre_days, " pseudo-event days [", w[1], ", ", w[2],
      "]\n", sep = "")
  co <- x$coefficients
  cat(sprintf("%16s %12s %10s %12s\n", "term", "estimate", "cdf p", "param p"))
  for (i in seq_len(nrow(co)))
    cat(sprintf("%16s %12.5g %10.3f %12.3f\n", co$term[i], co$estimate[i],
                co$p_cdf[i], co$p_parametric[i]))
  invisible(x)
}

#' @export
summary.fes_cse <- function(object, ...) {
  out <- object$coefficients
  attr(out, "method") <- object$method
  attr(out, "n_pre_days") <- object$diagnostics$n_pre_days
  out
}

#' @export
coef.fes_cse <- function(object, ...) object$params

#' Plot the pre-event coefficient distribution of a cross-sectional fit
#'
#' Histogram of the pseudo-event-day slopes for one characteristic with the
#' event-day slope marked — the "past is prologue" benchmark: the event slope
#' is attributable to the event only if it lies in the tails of the pre-event
#' distribution.
#'
#' @param x An `fes_cse`.
#' @param which Characteristic to plot: a term name or its index in
#'   `coef(x)` (default the first non-intercept characteristic).
#' @param ... Passed to the underlying histogram call.
#' @export
plot.fes_cse <- function(x, which = NULL, ...) {
  terms <- names(x$params)
  if (is.null(which))
    which <- if (length(terms) > 1L) terms[1L] else terms[1L]
  j <- if (is.character(which)) match(which, terms) else as.integer(which)
  if (is.na(j) || j < 1L || j > length(terms))
    stop("`which` does not match a coefficient")
  pre <- x$pre_betas[, j]
  b <- x$params[[j]]
  op <- graphics::par(mar = c(4, 4, 2, 1))
  on.exit(graphics::par(op))
  graphics::hist(pre, breaks = "FD", col = "grey85", border = "white",
                 xlim = range(pre, b, na.rm = TRUE),
                 xlab = paste0("slope on ", terms[j], " (pre-event days)"),
                 main = paste0("cross-sectional benchmark: ", terms[j]), ...)
  graphics::abline(v = b, col = "firebrick", lwd = 2)
  graphics::legend("topright",
                   sprintf("event slope = %.4g\nCDF p = %.3f", b, x$p_cdf[[j]]),
                   text.col = "firebrick", bty = "n")
  invisible(x)
}
