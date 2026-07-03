# S3 methods for fes_fit --------------------------------------------------

#' @export
print.fes_fit <- function(x, ...) {
  w <- x$conventions$window
  cat("feventr fit: method '", x$method, "', ", x$diagnostics$n_treated,
      " treated / ", x$diagnostics$n_donors, " donors\n", sep = "")
  cat("event window [", w[1], ", ", w[2], "], estimation window [",
      x$conventions$est_window[1], ", ", x$conventions$est_window[2],
      "], returns: ", x$conventions$returns, "\n", sep = "")
  cat("ATT (avg over window): ", formatC(x$att_avg, digits = 4, format = "f"),
      if (!is.null(x$att_avg_se))
        paste0(" (se ", formatC(x$att_avg_se, digits = 4, format = "f"), ")"),
      "\n", sep = "")
  if (identical(x$se$method, "conformal")) {
    if (anyNA(x$se$avg_ci))
      # the point-estimate constant is rejected by the joint test, so the CI
      # cannot be formed as an interval bracketing it (avg_p tests h0 = 0, a
      # different hypothesis, so it is deliberately not reported here)
      cat(format(100 * x$se$level), "% conformal CI (constant effect): ",
          "unavailable — the constant-effect model is rejected at the point ",
          "estimate\n", sep = "")
    else
      cat(format(100 * x$se$level), "% conformal CI (constant effect): [",
          formatC(x$se$avg_ci[1], digits = 4, format = "f"), ", ",
          formatC(x$se$avg_ci[2], digits = 4, format = "f"), "], joint p = ",
          formatC(x$se$avg_p, digits = 3, format = "f"), "\n", sep = "")
  }
  cat("Cumulative effect (", x$conventions$cumulate, ") at end of window: ",
      formatC(unname(x$car[length(x$car)]), digits = 4, format = "f"),
      "\n", sep = "")
  invisible(x)
}

#' @export
summary.fes_fit <- function(object, ...) {
  out <- data.frame(
    event_time = as.numeric(names(object$att)),
    att = unname(object$att),
    se = if (!is.null(object$se$att)) unname(object$se$att) else NA_real_,
    car = unname(object$car)
  )
  if (identical(object$se$method, "conformal")) {
    out$se <- NULL
    out$lower <- unname(object$se$ci[, 1])
    out$upper <- unname(object$se$ci[, 2])
    out$p <- unname(object$se$p)
  }
  attr(out, "method") <- object$method
  out
}

#' @export
coef.fes_fit <- function(object, cumulative = FALSE, ...) {
  if (cumulative) object$car else object$att
}

#' @export
vcov.fes_fit <- function(object, ...) {
  if (is.null(object$se$att)) stop("fit has no stored standard errors")
  if (!is.null(object$se$draws)) {
    v <- stats::cov(object$se$draws)
  } else {
    v <- diag(object$se$att^2, length(object$se$att))
  }
  dimnames(v) <- list(names(object$att), names(object$att))
  v
}

#' @export
confint.fes_fit <- function(object, parm, level = 0.95, ...) {
  if (identical(object$se$method, "conformal")) {
    if (!isTRUE(all.equal(level, object$se$level)))
      stop("conformal CIs are computed at fit time (level ", object$se$level,
           "); refit with `level = ", level, "`")
    return(object$se$ci)
  }
  if (is.null(object$se$att)) stop("fit has no stored standard errors")
  # Use the exact t critical value when the SEs carry degrees of freedom
  # (tstat: firm-day/pooled df; placebo: reps - 1); fall back to the normal
  # only when df is unavailable (e.g. gsynth bootstrap SEs).
  df <- object$se$df
  z <- if (!is.null(df) && is.finite(df)) stats::qt(1 - (1 - level) / 2, df = df)
       else stats::qnorm(1 - (1 - level) / 2)
  out <- cbind(object$att - z * object$se$att, object$att + z * object$se$att)
  colnames(out) <- paste(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2),
                                trim = TRUE), "%")
  out
}
