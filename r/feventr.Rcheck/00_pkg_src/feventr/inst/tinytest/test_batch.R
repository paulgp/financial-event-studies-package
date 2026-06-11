# event_study_batch: many staggered events, cross-event SEs, donor exclusion

set.seed(5)
# build a staggered panel: 3 events at different dates, distinct treated units
n <- 60; T <- 120
ret <- matrix(rnorm(n * T, 0, 0.01), n, T)
ev <- data.frame(unit = as.character(c(1, 2, 11, 12, 21)),
                 event_time = c(60, 60, 80, 80, 100),
                 event = c("e1", "e1", "e2", "e2", "e3"),
                 bucket = c("x", "x", "x", "x", "y"))
tau <- 0.05
for (i in seq_len(nrow(ev)))
  ret[as.integer(ev$unit[i]), ev$event_time[i]] <- ret[as.integer(ev$unit[i]), ev$event_time[i]] + tau
long <- data.frame(id = rep(1:n, times = T), t = rep(1:T, each = n),
                   ret = as.vector(ret))

b <- feventr::event_study_batch(long, "id", "t", "ret", events = ev,
                                method = "mean", window = c(0, 5),
                                est_window = c(-40, -1), returns = "simple",
                                se = "cross")
expect_inherits(b, "fes_batch")
expect_equal(nrow(b$atts), 3L)
expect_equal(names(b$att), as.character(0:5))
# event-day cross-event mean recovers tau
expect_true(abs(b$att[["0"]] - tau) < 0.01)
# cross-event SE equals sd/sqrt(n) of the per-event paths
expect_equal(unname(b$se$att), unname(apply(b$atts, 2, sd) / sqrt(3)), tolerance = 1e-12)
# per-event ATT equals a single-event fit with the same donor pool: for e1
# (loaded range 20..65) the other events (t=80, t=100) are outside the range,
# so their units remain clean donors — only e1's own units are excluded
f1 <- feventr::event_study(long, "id", "t", "ret",
                           treated = c("1", "2"), event_time = 60,
                           donors = setdiff(as.character(1:n), c("1", "2")),
                           method = "mean", window = c(0, 5),
                           est_window = c(-40, -1), returns = "simple",
                           se = "none")
expect_equal(unname(b$atts["e1", ]), unname(f1$att), tolerance = 1e-12)

# exclusion binds when windows overlap: for e2 (range 40..85) unit 1/2's
# event at 60 contaminates, so they must NOT be donors of e2
f2 <- feventr::event_study(long, "id", "t", "ret",
                           treated = c("11", "12"), event_time = 80,
                           donors = setdiff(as.character(1:n),
                                            c("1", "2", "11", "12")),
                           method = "mean", window = c(0, 5),
                           est_window = c(-40, -1), returns = "simple",
                           se = "none")
expect_equal(unname(b$atts["e2", ]), unname(f2$att), tolerance = 1e-12)

# status + extra columns carried:
expect_equal(b$events$status, rep("ok", 3))
expect_true("bucket" %in% names(b$events))

# grouped summary
s <- summary(b, by = "bucket", horizon = 0)
expect_equal(sort(unique(s$group)), c("x", "y"))
expect_equal(s$n_events[s$group == "x"][1], 2L)

# methods
expect_stdout(print(b), pattern = "3/3 events")
expect_equal(coef(b), b$att)

# sc engine through batch
bsc <- feventr::event_study_batch(long, "id", "t", "ret", events = ev,
                                  method = "sc", window = c(0, 5),
                                  est_window = c(-40, -1), returns = "simple")
expect_equal(nrow(bsc$atts), 3L)
expect_true(all(is.finite(bsc$att)))

# plots run silently to a null device
tf <- tempfile(fileext = ".png")
grDevices::png(tf)
expect_silent(plot(b))
expect_silent(plot(bsc, what = "car"))
grDevices::dev.off()
unlink(tf)
