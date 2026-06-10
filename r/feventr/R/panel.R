# Panel layer: long data -> balanced event-window matrices ---------------------
#
# All R-side data handling (data.table, column lookup, time alignment) stops
# here. Estimator engines (R/estimators.R) see plain base matrices in the
# synthdid orientation: rows = units (controls 1..N0, treated below), columns
# = periods (estimation columns 1..T0, event-window columns after).

# Convert a long panel into the matrix form consumed by the engines.
#
# Event time is positional: the offset of each panel time from `event_time` in
# the sorted unique times of `data` (trading periods, not calendar arithmetic),
# so daily and monthly panels work identically. A gap between `est_window` and
# `window` (e.g. Geithner's excluded placebo days -30..-1) is allowed: gap
# columns are simply not loaded.
fes_panel <- function(data, unit, time, ret, treated, event_time,
                      window, est_window, donors = NULL) {
  if (est_window[1] > est_window[2] || window[1] > window[2])
    stop("windows must be increasing (start, end) pairs")
  if (est_window[2] >= window[1])
    stop("`est_window` must end before the event `window` starts")

  dt <- data.table::as.data.table(data)[, c(unit, time, ret), with = FALSE]
  data.table::setnames(dt, c("unit", "time", "ret"))
  times <- sort(unique(dt[["time"]]))
  pos0 <- match(event_time, times)
  if (is.na(pos0)) stop("`event_time` not found among panel times")
  offset <- seq_along(times) - pos0
  keep <- (offset >= est_window[1] & offset <= est_window[2]) |
    (offset >= window[1] & offset <= window[2])
  times_keep <- times[keep]
  dt <- dt[dt[["time"]] %in% times_keep]
  if (anyDuplicated(dt, by = c("unit", "time")))
    stop("duplicate unit-time rows in `data`")

  W <- data.table::dcast(dt, unit ~ time, value.var = "ret")
  units <- as.character(W[["unit"]])
  Y <- as.matrix(W[, -1L, with = FALSE])
  Y <- Y[, match(as.character(times_keep), colnames(Y)), drop = FALSE]
  rownames(Y) <- units

  treated <- as.character(treated)
  if (!any(units %in% treated)) stop("no treated units found in `data`")
  pool <- if (is.null(donors)) setdiff(units, treated)
          else setdiff(intersect(units, as.character(donors)), treated)

  complete <- rowSums(is.na(Y)) == 0L
  dropped <- data.frame(unit = character(), reason = character())
  drop_tr <- intersect(treated, units[!complete])
  if (length(drop_tr)) {
    warning(length(drop_tr), " treated unit(s) dropped: incomplete history")
    dropped <- rbind(dropped, data.frame(unit = drop_tr, reason = "treated: incomplete history"))
  }
  drop_co <- intersect(pool, units[!complete])
  if (length(drop_co))
    dropped <- rbind(dropped, data.frame(unit = drop_co, reason = "donor: incomplete history"))

  tr <- intersect(units[complete], treated)
  co <- intersect(units[complete], pool)
  if (!length(tr)) stop("no treated units with complete history over the windows")
  if (!length(co)) stop("no donor units with complete history over the windows")

  Y <- Y[c(co, tr), , drop = FALSE]
  list(Y = Y, N0 = length(co), T0 = sum(offset[keep] < window[1]),
       units = c(co, tr), treated = tr, times = offset[keep],
       time_values = times_keep, dropped = dropped)
}

# Align a user-supplied factor table (column `time` + numeric factor columns)
# to the panel's kept columns. Returns a T x K matrix.
align_factors <- function(factors, time, panel) {
  if (is.null(factors)) stop("`factors` is required for this method")
  fdt <- as.data.frame(factors)
  if (!time %in% names(fdt))
    stop("`factors` must contain the panel time column '", time, "'")
  idx <- match(panel$time_values, fdt[[time]])
  if (anyNA(idx))
    stop("`factors` does not cover every panel period in the windows")
  F <- as.matrix(fdt[idx, setdiff(names(fdt), time), drop = FALSE])
  if (!is.numeric(F) || !ncol(F)) stop("`factors` must have numeric factor columns")
  F
}
