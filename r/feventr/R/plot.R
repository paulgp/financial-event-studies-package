# Plot methods (base graphics; deterministic for tinysnapshot) -------------

#' Plot an event-study fit
#'
#' @param x An `fes_fit`.
#' @param what `"att"` (ATT path with CI band, default), `"car"` (cumulative
#'   effect path), `"paths"` (cumulated treated vs synthetic return over the
#'   full loaded range), or `"weights"` (donor weights, synthetic methods).
#' @param ci Draw a confidence band when SEs are stored.
#' @param level Confidence level.
#' @param ... Passed to the underlying plot call.
#' @export
plot.fes_fit <- function(x, what = c("att", "car", "paths", "weights"),
                         ci = TRUE, level = 0.95, ...) {
  what <- match.arg(what)
  op <- graphics::par(mar = c(4, 4, 2, 1))
  on.exit(graphics::par(op))
  if (what == "att" || what == "car") {
    y <- if (what == "att") x$att else x$car
    t <- as.numeric(names(y))
    # band: stored conformal CI bounds if present (level fixed at fit time),
    # otherwise a normal band from the stored SEs
    has_ci <- !is.null(x$se$ci)
    band <- ci && what == "att" && (has_ci || !is.null(x$se$att))
    if (band) {
      if (has_ci) {
        lo <- x$se$ci[, 1]; hi <- x$se$ci[, 2]
      } else {
        z <- stats::qnorm(1 - (1 - level) / 2)
        lo <- y - z * x$se$att; hi <- y + z * x$se$att
      }
      band <- all(is.finite(lo)) && all(is.finite(hi))
    }
    ylim <- if (band) range(lo, hi, 0) else range(y, 0)
    plot(t, y, type = "n", ylim = ylim, xlab = "Event time",
         ylab = if (what == "att") "ATT" else
           paste0("Cumulative effect (", x$conventions$cumulate, ")"),
         main = paste0(x$method, if (what == "car") " - cumulative"), ...)
    graphics::abline(h = 0, col = "grey60")
    if (band)
      graphics::polygon(c(t, rev(t)), c(lo, rev(hi)),
                        col = grDevices::adjustcolor("steelblue", 0.25),
                        border = NA)
    graphics::lines(t, y, type = "b", pch = 19, col = "steelblue")
  } else if (what == "paths") {
    t <- as.numeric(names(x$paths$treated))
    cum <- function(r) cumprod(1 + r) - 1
    y1 <- cum(x$paths$treated)
    y0 <- cum(x$paths$synthetic)
    plot(t, y1, type = "l", lwd = 2, col = "firebrick",
         ylim = range(y1, y0), xlab = "Event time",
         ylab = "Cumulative return", main = x$method, ...)
    graphics::lines(t, y0, lwd = 2, lty = 2, col = "grey30")
    graphics::abline(v = x$conventions$window[1], col = "grey60", lty = 3)
    graphics::legend("topleft", c("treated", "synthetic"), lwd = 2,
                     lty = c(1, 2), col = c("firebrick", "grey30"), bty = "n")
  } else {
    w <- x$weights$omega
    if (is.null(w)) stop("fit has no donor weights (synthetic methods only)")
    w <- sort(w[w > 1e-6], decreasing = TRUE)
    if (length(w) > 30) w <- w[1:30]
    graphics::dotchart(rev(w), labels = rev(names(w)), pch = 19,
                       col = "steelblue", xlab = "Donor weight",
                       main = paste0(x$method, " - active donors"), ...)
  }
  invisible(x)
}

#' Plot cross-event average effects from a batch fit
#'
#' @param x An `fes_batch`.
#' @param what `"att"` or `"car"`.
#' @param ci Draw the cross-event confidence band.
#' @param level Confidence level.
#' @param ... Passed through.
#' @export
plot.fes_batch <- function(x, what = c("att", "car"), ci = TRUE,
                           level = 0.95, ...) {
  what <- match.arg(what)
  y <- if (what == "att") x$att else x$car
  s <- if (what == "att") x$se$att else x$se$car
  t <- as.numeric(names(y))
  op <- graphics::par(mar = c(4, 4, 2, 1))
  on.exit(graphics::par(op))
  band <- ci && !is.null(s)
  z <- stats::qnorm(1 - (1 - level) / 2)
  ylim <- if (band) range(y - z * s, y + z * s, 0) else range(y, 0)
  plot(t, y, type = "n", ylim = ylim, xlab = "Event time",
       ylab = if (what == "att") "ATT (cross-event avg)" else "Cumulative effect",
       main = paste0(x$method, " - ", sum(x$events$status == "ok"), " events"),
       ...)
  graphics::abline(h = 0, col = "grey60")
  if (band)
    graphics::polygon(c(t, rev(t)), c(y - z * s, rev(y + z * s)),
                      col = grDevices::adjustcolor("steelblue", 0.25),
                      border = NA)
  graphics::lines(t, y, type = "b", pch = 19, col = "steelblue")
  invisible(x)
}
