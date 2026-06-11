# Beta diagnostics --------------------------------------------------------

#' Factor-loading comparison of treated and control units
#'
#' Estimates per-unit OLS factor loadings over the fit's estimation window
#' and tabulates them by group: treated, control (equal-weighted), and —
#' when the fit carries synthetic-control unit weights — the omega-weighted
#' control portfolio (the paper's Tables 3, 4 and 7 diagnostics).
#'
#' @param object An `fes_fit` from [event_study()] with `keep_data = TRUE`.
#' @param factors Data frame: a time column (named as in the original panel,
#'   or the first column) plus numeric factor columns. As in [event_study()],
#'   any excess-return adjustment is the caller's responsibility.
#' @param time Name of the time column in `factors`; default its first column.
#' @return A data.frame of class `fes_betas`: one row per group x statistic
#'   (mean, sd, median), one column per loading (alpha first).
#' @export
event_betas <- function(object, factors, time = NULL) {
  stopifnot(inherits(object, "fes_fit"))
  p <- object$panel
  if (is.null(p)) stop("fit was created with keep_data = FALSE")
  if (is.null(time)) time <- names(factors)[1]
  F <- align_factors(factors, time, p)
  pre <- seq_len(p$T0)
  X <- cbind(alpha = 1, F[pre, , drop = FALSE])
  B <- t(qr.solve(X, t(p$Y[, pre, drop = FALSE])))   # units x (1+K)
  colnames(B) <- c("alpha", colnames(F))
  tr <- rownames(B) %in% p$treated
  grp <- function(M, label, w = NULL) {
    stats <- if (is.null(w)) {
      rbind(mean = colMeans(M), sd = apply(M, 2, stats::sd),
            median = apply(M, 2, stats::median))
    } else {
      rbind(mean = apply(M, 2, stats::weighted.mean, w = w),
            sd = rep(NA_real_, ncol(M)), median = rep(NA_real_, ncol(M)))
    }
    data.frame(group = label, stat = rownames(stats), stats,
               row.names = NULL, check.names = FALSE)
  }
  out <- rbind(grp(B[tr, , drop = FALSE], "treated"),
               grp(B[!tr, , drop = FALSE], "control"))
  om <- object$weights$omega
  if (!is.null(om))
    out <- rbind(out, grp(B[!tr, , drop = FALSE][names(om), , drop = FALSE],
                          paste0("control_weighted_", object$method), w = om))
  attr(out, "n") <- c(treated = sum(tr), control = sum(!tr))
  class(out) <- c("fes_betas", "data.frame")
  out
}

#' @export
print.fes_betas <- function(x, digits = 3, ...) {
  n <- attr(x, "n")
  cat("Factor loadings over the estimation window (", n["treated"],
      " treated, ", n["control"], " control units)\n", sep = "")
  print.data.frame(x, digits = digits, row.names = FALSE)
  invisible(x)
}
