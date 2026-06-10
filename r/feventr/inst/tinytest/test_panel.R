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

# donors= restricts the pool
p4 <- feventr:::fes_panel(d, "id", "t", "ret", treated = "u1",
                          event_time = 25, window = c(0, 3),
                          est_window = c(-20, -5),
                          donors = c("u3", "u4"))
expect_equal(p4$units[seq_len(p4$N0)], c("u3", "u4"))
