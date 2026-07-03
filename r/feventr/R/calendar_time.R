# Calendar-time portfolio estimator (Jaffe 1974; Mandelker 1974; Fama 1998) ----

#' Calendar-time portfolio event study
#'
#' For every calendar period, forms a portfolio of all units whose event
#' occurred within `window` event-time offsets, and regresses the portfolio
#' return on the supplied factors. The intercept is the average abnormal
#' return per calendar period — the Jaffe-Mandelker calendar-time portfolio
#' estimator advocated by Fama (1998) for event studies with many staggered
#' events, where overlapping event windows make event-time CARs
#' cross-sectionally dependent; aggregating to one portfolio return per
#' calendar period absorbs that dependence into the time series.
#'
#' @param data,unit,time,ret As in [event_study()]: one long panel covering
#'   every event's window.
#' @param events Data frame with columns `unit` and `event_time`, one row per
#'   event, as in [event_study_batch()]. A unit appearing in several
#'   overlapping windows enters the portfolio once per period.
#' @param window Portfolio membership window in event-time offsets, e.g.
#'   `c(0, 10)`: a unit is held from `window[1]` to `window[2]` trading
#'   periods after its event.
#' @param factors Data frame with the panel `time` column plus numeric factor
#'   columns covering every portfolio period (as in [event_study()]; the
#'   excess-return convention is the caller's responsibility). `NULL` fits an
#'   intercept-only model: alpha is the mean portfolio return, for use when
#'   `ret` is already an abnormal or excess return.
#' @param returns `"simple"` or `"log"`; required, no default.
#' @param align `"position"` (default; event/calendar times count positions
#'   among the panel's sorted unique times) or `"value"` (numeric time values
#'   used directly), as in [event_study()].
#' @param weight Optional column name in `data` holding portfolio weights
#'   (e.g. lagged market cap) for a value-weighted portfolio; default
#'   equal-weighted. Rows with a missing weight are excluded from the
#'   portfolio (as with missing returns).
#' @param min_units Calendar periods whose portfolio holds fewer units are
#'   dropped from the regression (Fama (1998) uses 10).
#' @param se `"ols"` (classical) or `"nw"` (Newey-West with Bartlett kernel).
#' @param lag Newey-West lag; default `floor(4 * (n/100)^(2/9))`.
#' @return An object of class `fes_caltime`: `$alpha`/`$alpha_se` (per-period
#'   abnormal return), `$coefficients` (regression table), `$portfolio`
#'   (calendar-period series: return, units held, abnormal return), `$nobs`.
#' @export
calendar_time <- function(data, unit, time, ret, events, window = c(0, 10),
                          factors = NULL, returns,
                          align = c("position", "value"), weight = NULL,
                          min_units = 1L, se = c("ols", "nw"), lag = NULL) {
  returns <- match.arg(returns, c("simple", "log"))
  align <- match.arg(align)
  se <- match.arg(se)
  events <- as.data.frame(events)
  stopifnot(all(c("unit", "event_time") %in% names(events)),
            length(window) == 2L, window[1] <= window[2])

  times <- sort(unique(data[[time]]))
  pos_of <- function(v) if (align == "value") as.numeric(v) else match(v, times)
  ev_pos <- pos_of(events$event_time)
  if (anyNA(ev_pos)) stop("some `event_time`s not found among panel times")

  # (unit, period) membership pairs, deduplicated across overlapping events
  offs <- seq.int(window[1], window[2])
  units_all <- unique(c(as.character(data[[unit]]), as.character(events$unit)))
  d_u <- match(as.character(data[[unit]]), units_all)
  d_p <- pos_of(data[[time]])
  m_u <- rep(match(as.character(events$unit), units_all), each = length(offs))
  m_p <- rep(ev_pos, each = length(offs)) + offs
  lo <- min(d_p, m_p)
  P <- max(d_p, m_p) - lo + 2
  if (!is.null(weight) && !weight %in% names(data))
    stop("`weight` column not found in `data`")
  keep <- (d_u * P + (d_p - lo)) %in% unique(m_u * P + (m_p - lo)) &
    !is.na(data[[ret]])
  # A single NA weight would make that calendar period's portfolio return NA
  # and propagate through qr.solve into every coefficient; drop NA-weight rows
  # (a unit's first listed period commonly lacks a lagged market cap).
  if (!is.null(weight)) keep <- keep & !is.na(data[[weight]])

  r <- data[[ret]][keep]
  p <- d_p[keep]
  w <- if (is.null(weight)) rep(1, length(r)) else data[[weight]][keep]
  if (!length(r)) stop("no panel observations fall in any event window")

  n_t <- drop(rowsum(rep(1L, length(p)), p))
  pr <- drop(rowsum(r * w, p)) / drop(rowsum(w, p))
  pos <- as.numeric(names(n_t))
  ok <- n_t >= min_units
  n_dropped <- sum(!ok)
  if (!any(ok)) stop("no calendar period holds `min_units` units")
  tvals <- if (align == "value") pos[ok] else times[pos[ok]]

  F <- NULL
  if (!is.null(factors)) {
    fdt <- as.data.frame(factors)
    if (!time %in% names(fdt))
      stop("`factors` must contain the panel time column '", time, "'")
    idx <- match(tvals, fdt[[time]])
    if (anyNA(idx)) stop("`factors` does not cover every portfolio period")
    F <- as.matrix(fdt[idx, setdiff(names(fdt), time), drop = FALSE])
    if (!is.numeric(F) || !ncol(F))
      stop("`factors` must have numeric factor columns")
  }

  y <- pr[ok]
  X <- cbind(alpha = rep(1, length(y)), F)
  n <- length(y)
  k <- ncol(X)
  if (n <= k) stop("fewer calendar periods than regression coefficients")
  B <- qr.solve(X, y)
  e <- y - as.vector(X %*% B)
  XX_inv <- solve(crossprod(X))
  V <- if (se == "ols") {
    XX_inv * sum(e^2) / (n - k)
  } else {
    if (is.null(lag)) lag <- floor(4 * (n / 100)^(2 / 9))
    Xe <- X * e
    S <- crossprod(Xe)
    for (l in seq_len(lag)) {
      G <- crossprod(Xe[seq_len(n - l), , drop = FALSE],
                     Xe[seq.int(l + 1L, n), , drop = FALSE])
      S <- S + (1 - l / (lag + 1)) * (G + t(G))
    }
    XX_inv %*% S %*% XX_inv
  }
  ses <- sqrt(diag(V))

  structure(list(
    alpha = unname(B[1]),
    alpha_se = unname(ses[1]),
    coefficients = data.frame(term = colnames(X), estimate = unname(B),
                              se = unname(ses),
                              t = unname(B / ses), row.names = NULL),
    portfolio = data.frame(time = tvals, ret = unname(y),
                           n_units = unname(n_t[ok]),
                           abnormal = unname(B[1] + e)),
    nobs = n,
    diagnostics = list(n_events = nrow(events), n_dropped_periods = n_dropped,
                       mean_units = mean(n_t[ok])),
    conventions = list(returns = returns, window = window, se = se,
                       lag = if (se == "nw") lag,
                       weight = if (is.null(weight)) "equal" else weight,
                       min_units = min_units),
    call = match.call()
  ), class = "fes_caltime")
}

#' @export
print.fes_caltime <- function(x, ...) {
  w <- x$conventions$window
  cat("feventr calendar-time portfolio: ", x$diagnostics$n_events,
      " events, window [", w[1], ", ", w[2], "]\n", sep = "")
  cat(x$nobs, " calendar periods (", x$conventions$weight,
      "-weighted, mean ", formatC(x$diagnostics$mean_units, digits = 1,
                                  format = "f"),
      " units/period)\n", sep = "")
  cat("alpha (abnormal return per period): ",
      formatC(x$alpha, digits = 4, format = "f"), " (se ",
      formatC(x$alpha_se, digits = 4, format = "f"), ", ",
      x$conventions$se, ")\n", sep = "")
  invisible(x)
}

#' @export
summary.fes_caltime <- function(object, ...) {
  out <- object$coefficients
  attr(out, "nobs") <- object$nobs
  out
}

#' @export
coef.fes_caltime <- function(object, ...) {
  stats::setNames(object$coefficients$estimate, object$coefficients$term)
}

#' Plot a calendar-time portfolio fit
#'
#' @param x An `fes_caltime`.
#' @param what `"car"` (default): cumulative abnormal portfolio return over
#'   calendar time, with the constant-alpha trend dashed. `"n_units"`:
#'   portfolio size over calendar time.
#' @param ... Passed to the underlying plot call.
#' @export
plot.fes_caltime <- function(x, what = c("car", "n_units"), ...) {
  what <- match.arg(what)
  op <- graphics::par(mar = c(4, 4, 2, 1))
  on.exit(graphics::par(op))
  t <- x$portfolio$time
  if (what == "car") {
    y <- cumsum(x$portfolio$abnormal)
    trend <- x$alpha * seq_along(y)
    plot(t, y, type = "l", lwd = 2, col = "steelblue",
         ylim = range(y, trend, 0), xlab = "Calendar time",
         ylab = "Cumulative abnormal return",
         main = "calendar-time portfolio", ...)
    graphics::abline(h = 0, col = "grey60")
    graphics::lines(t, trend, lty = 2, col = "grey30")
    graphics::legend("topleft", c("portfolio", "constant alpha"), lwd = c(2, 1),
                     lty = c(1, 2), col = c("steelblue", "grey30"), bty = "n")
  } else {
    plot(t, x$portfolio$n_units, type = "s", lwd = 2, col = "steelblue",
         ylim = c(0, max(x$portfolio$n_units)), xlab = "Calendar time",
         ylab = "Units in portfolio", main = "calendar-time portfolio", ...)
  }
  invisible(x)
}
