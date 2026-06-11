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
  if (is.null(object$se$att)) stop("fit has no stored standard errors")
  z <- stats::qnorm(1 - (1 - level) / 2)
  out <- cbind(object$att - z * object$se$att, object$att + z * object$se$att)
  colnames(out) <- paste(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2),
                                trim = TRUE), "%")
  out
}
