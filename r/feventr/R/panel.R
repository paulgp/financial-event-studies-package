# Panel layer: long data -> balanced event-window matrices ---------------------
#
# All R-side data handling (data.table, column lookup, time alignment) stops
# here. Estimator engines (R/estimators.R) see plain base matrices in the
# synthdid orientation: rows = units (controls 1..N0, treated below), columns
# = periods (estimation columns 1..T0, event-window columns after).

# Convert a long panel into the matrix form consumed by the engines.
#
# Event time: with `align = "position"` (default) the offset of each panel
# time is its position from `event_time` in the sorted unique times of `data`
# (trading periods, not calendar arithmetic), so daily and monthly panels work
# identically. With `align = "value"` the numeric time values are used
# directly as offsets (time - event_time) — required when the time column is
# already an event-time index AND periods are missing from the data (e.g. the
# Geithner panels with the placebo days -30..-1 pre-deleted: positionally,
# dif = -256 would silently land at offset -226). A gap between `est_window`
# and `window` is allowed: gap columns are simply not loaded.
fes_panel <- function(data, unit, time, ret, treated, event_time,
                      window, est_window, donors = NULL,
                      align = c("position", "value")) {
  align <- match.arg(align)
  if (est_window[1] > est_window[2] || window[1] > window[2])
    stop("windows must be increasing (start, end) pairs")
  if (est_window[2] >= window[1])
    stop("`est_window` must end before the event `window` starts")

  miss <- setdiff(c(unit, time, ret), names(data))
  if (length(miss))
    stop("column(s) not found in `data`: ", paste(miss, collapse = ", "))
  tvals <- data[[time]]
  times <- sort(unique(tvals))
  pos0 <- match(event_time, times)
  if (is.na(pos0)) stop("`event_time` not found among panel times")
  offset <- if (align == "value") {
    if (!is.numeric(times)) stop("align = 'value' needs a numeric time column")
    as.numeric(times) - as.numeric(event_time)
  } else {
    seq_along(times) - pos0
  }
  in_win <- offset >= window[1] & offset <= window[2]
  keep <- (offset >= est_window[1] & offset <= est_window[2]) | in_win
  # The requested event window must be fully covered by the panel's periods.
  # Filtering `keep` alone would silently truncate the ATT path for an event
  # near a sample edge (conventions$window would still report the full window),
  # and a window with no periods at all makes `post` descend out of bounds in
  # event_study(). A partial *estimation* window is allowed (fewer pre-periods
  # still identifies; positional alignment on gapped data trims it by design).
  wo <- offset[in_win]
  if (!length(wo) || min(wo) > window[1] || max(wo) < window[2])
    stop("panel does not cover event window [", window[1], ", ", window[2],
         "] for event_time ", event_time, " (available offsets in window: ",
         if (length(wo)) paste0(min(wo), "..", max(wo)) else "none", ")")
  times_keep <- times[keep]
  # Trim to the loaded windows BEFORE materializing anything: copying the
  # whole input first meant every batch worker briefly held a full-panel
  # duplicate (48 years x 66.7M rows in the index-inclusion application),
  # which is what OOM-killed parallel runs. Only the kept rows of the three
  # needed columns are ever copied, and the time key becomes a small integer
  # period index so grouping/pivoting never runs on Date comparison paths.
  # match() on classed vectors (Date, POSIXct) routes through mtfrm(), i.e.
  # a full format() of every row — the ~50s-per-call cost on the 66.7M-row
  # panel. Matching the unclassed keys is the same join on the raw numbers.
  ti <- if (is.object(tvals)) match(unclass(tvals), unclass(times))
        else match(tvals, times)
  rows <- which(keep[ti])
  dt <- data.table::data.table(unit = data[[unit]][rows], ti = ti[rows],
                               ret = data[[ret]][rows])
  if (anyDuplicated(dt, by = c("unit", "ti")))
    stop("duplicate unit-time rows in `data`")

  W <- data.table::dcast(dt, unit ~ ti, value.var = "ret")
  units <- as.character(W[["unit"]])
  Y <- as.matrix(W[, -1L, with = FALSE])
  Y <- Y[, match(as.character(which(keep)), colnames(Y)), drop = FALSE]
  # engines and fit objects keep seeing the caller's time values in the
  # column names (sdid's lambda weights are named from them)
  colnames(Y) <- as.character(times_keep)
  rownames(Y) <- units

  treated <- as.character(treated)
  if (!any(units %in% treated)) stop("no treated units found in `data`")
  pool <- if (is.null(donors)) setdiff(units, treated)
          else setdiff(intersect(units, as.character(donors)), treated)

  complete <- rowSums(is.na(Y)) == 0L
  dropped <- data.frame(unit = character(), reason = character())
  # treated ids absent from the panel (typo, ticker-case mismatch) would never
  # reach the incomplete-history drop below and would silently shrink the
  # treated set; flag them explicitly and record the drop
  missing_tr <- setdiff(treated, units)
  if (length(missing_tr)) {
    warning(length(missing_tr), " treated unit(s) not found in `data`: ",
            paste(missing_tr, collapse = ", "))
    dropped <- rbind(dropped, data.frame(unit = missing_tr,
                                         reason = "treated: not in panel"))
  }
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
# to a target vector of time values. Returns a length(time_values) x K matrix.
# Shared by event_study() (panel columns) and calendar_time() (portfolio
# periods) so the coverage/format checks stay in one place.
align_factors <- function(factors, time, time_values) {
  if (is.null(factors)) stop("`factors` is required for this method")
  fdt <- as.data.frame(factors)
  if (!time %in% names(fdt))
    stop("`factors` must contain the panel time column '", time, "'")
  idx <- match(time_values, fdt[[time]])
  if (anyNA(idx))
    stop("`factors` does not cover every period in the windows")
  F <- as.matrix(fdt[idx, setdiff(names(fdt), time), drop = FALSE])
  if (!is.numeric(F) || !ncol(F)) stop("`factors` must have numeric factor columns")
  F
}
