# fes_panel: long data -> matrices

mk_long <- function(n = 8, T = 30) {
  expand.grid(id = paste0("u", seq_len(n)), t = seq_len(T),
              KEEP.OUT.ATTRS = FALSE) |>
    transform(ret = stats::rnorm(n * T, 0, 0.01))
}

set.seed(1)
d <- mk_long()
p <- feventr:::fes_panel(d, "id", "t", "ret", treated = c("u1", "u2"),
                         event_time = 25, window = c(0, 3),
                         est_window = c(-20, -5))
# gap between est_window and window is simply not loaded
expect_equal(p$times, c(-20:-5, 0:3))
expect_equal(p$T0, 16L)
expect_equal(p$N0, 6L)
# treated rows come last, in unit order
expect_equal(tail(p$units, 2), c("u1", "u2"))
expect_equal(dim(p$Y), c(8L, 20L))
# values land in the right cells
expect_equal(p$Y["u3", "20"], d$ret[d$id == "u3" & d$t == 20])

# incomplete donor history is dropped silently into `dropped`
d2 <- d[!(d$id == "u5" & d$t == 10), ]
p2 <- feventr:::fes_panel(d2, "id", "t", "ret", treated = c("u1", "u2"),
                          event_time = 25, window = c(0, 3),
                          est_window = c(-20, -5))
expect_equal(p2$N0, 5L)
expect_equal(p2$dropped$unit, "u5")

# incomplete treated history warns
d3 <- d[!(d$id == "u1" & d$t == 25), ]
expect_warning(
  feventr:::fes_panel(d3, "id", "t", "ret", treated = c("u1", "u2"),
                      event_time = 25, window = c(0, 3),
                      est_window = c(-20, -5)),
  "treated"
)

# errors: duplicates (within the loaded windows), overlapping windows,
# missing event time; duplicates outside the loaded windows are harmless
expect_error(
  feventr:::fes_panel(rbind(d, d[d$id == "u1" & d$t == 25, ]),
                      "id", "t", "ret", "u1", 25, c(0, 3), c(-20, -5)),
  "duplicate"
)
expect_silent(
  feventr:::fes_panel(rbind(d, d[d$id == "u3" & d$t == 1, ]),
                      "id", "t", "ret", c("u1", "u2"), 25, c(0, 3), c(-20, -5))
)
expect_error(
  feventr:::fes_panel(d, "id", "t", "ret", "u1", 25, c(0, 3), c(-20, 1)),
  "must end before"
)
expect_error(
  feventr:::fes_panel(d, "id", "t", "ret", "u1", 99, c(0, 3), c(-20, -5)),
  "event_time"
)

# requested event window must be fully covered: an event near the sample edge
# errors instead of silently truncating the ATT path (issue 1)
expect_error(
  feventr:::fes_panel(d, "id", "t", "ret", "u1", event_time = 28,
                      window = c(0, 10), est_window = c(-20, -5)),
  "does not cover event window"
)
# a window with no periods at all (panel ends on event day, window starts at 1)
# would otherwise make `post` descend out of bounds downstream
expect_error(
  feventr:::fes_panel(d, "id", "t", "ret", "u1", event_time = 30,
                      window = c(1, 10), est_window = c(-20, -5)),
  "does not cover event window"
)
# a partial *estimation* window is still allowed (positional alignment on
# gapped data trims it by design; see the align tests below)
expect_silent(
  feventr:::fes_panel(d, "id", "t", "ret", c("u1", "u2"), event_time = 25,
                      window = c(0, 3), est_window = c(-40, -5))
)

# a treated id absent from the panel (typo) is flagged, not silently dropped,
# and recorded in `dropped` (issue 9)
expect_warning(
  p_miss <- feventr:::fes_panel(d, "id", "t", "ret",
                               treated = c("u1", "u2", "nope"), event_time = 25,
                               window = c(0, 3), est_window = c(-20, -5)),
  "not found"
)
expect_true("nope" %in% p_miss$dropped$unit)
expect_equal(sort(p_miss$treated), c("u1", "u2"))

# donors= restricts the pool
p4 <- feventr:::fes_panel(d, "id", "t", "ret", treated = "u1",
                          event_time = 25, window = c(0, 3),
                          est_window = c(-20, -5),
                          donors = c("u3", "u4"))
expect_equal(p4$units[seq_len(p4$N0)], c("u3", "u4"))

# align = "value": when the time column is already an event-time index and
# some periods are pre-deleted from the data (the Geithner placebo-window
# trap), positional counting silently shifts the windows; value alignment
# must not. Build t = -20..3 with -4..-1 deleted:
d5 <- mk_long(6, 24)
d5$t <- d5$t - 21                      # t = -20..3
d5 <- d5[!(d5$t %in% -4:-1), ]
pv <- feventr:::fes_panel(d5, "id", "t", "ret", treated = "u1",
                          event_time = 0, window = c(0, 3),
                          est_window = c(-20, -5), align = "value")
expect_equal(pv$times, c(-20:-5, 0:3))
expect_equal(pv$T0, 16L)
# positional alignment on the same gapped data shifts: t = -20 sits 16
# positions before 0, so est_window c(-20,-5) loads only t -20..-9
pp <- feventr:::fes_panel(d5, "id", "t", "ret", treated = "u1",
                          event_time = 0, window = c(0, 3),
                          est_window = c(-20, -5), align = "position")
expect_equal(pp$T0, 12L)

# Date time keys (issues 21 + 25): the trimmed integer-index path must
# reproduce the numeric-key panel exactly, keep Date time_values, and keep
# the caller's time values as Y's column names
dts <- as.Date("2020-01-01") + 0:79
set.seed(4)
longd <- data.frame(id = rep(paste0("u", 1:12), each = 80),
                    d = rep(dts, times = 12), ret = rnorm(960, 0, 0.01))
pd <- feventr:::fes_panel(longd, "id", "d", "ret", treated = "u12",
                          event_time = dts[61], window = c(0, 5),
                          est_window = c(-40, -11))
longi <- transform(longd, d = as.integer(d))
pin <- feventr:::fes_panel(longi, "id", "d", "ret", treated = "u12",
                           event_time = as.integer(dts[61]), window = c(0, 5),
                           est_window = c(-40, -11))
expect_equivalent(pd$Y, pin$Y)
expect_equal(pd$N0, pin$N0)
expect_equal(pd$T0, pin$T0)
expect_equal(pd$times, pin$times)
expect_inherits(pd$time_values, "Date")
expect_equal(colnames(pd$Y), as.character(pd$time_values))

# duplicate unit-time rows inside the loaded windows are still a hard error
# after the trim-first change; duplicates outside them never were
dup_in <- rbind(longd, longd[longd$id == "u3" & longd$d == dts[62], ])
expect_error(feventr:::fes_panel(dup_in, "id", "d", "ret", treated = "u12",
                                 event_time = dts[61], window = c(0, 5),
                                 est_window = c(-40, -11)),
             pattern = "duplicate")
dup_out <- rbind(longd, longd[longd$id == "u3" & longd$d == dts[15], ])
expect_silent(feventr:::fes_panel(dup_out, "id", "d", "ret", treated = "u12",
                                  event_time = dts[61], window = c(0, 5),
                                  est_window = c(-40, -11)))

# missing columns get a direct error
expect_error(feventr:::fes_panel(longd, "id", "nope", "ret", treated = "u12",
                                 event_time = dts[61], window = c(0, 5),
                                 est_window = c(-40, -11)),
             pattern = "not found")
