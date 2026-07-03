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
