# calendar_time(): hand-built portfolio + lm() reference, effect recovery,
# weighting/min_units/dedup mechanics.

set.seed(31)
n_u <- 30; n_t <- 120
long <- data.frame(id = rep(sprintf("u%02d", 1:n_u), times = n_t),
                   t = rep(1:n_t, each = n_u),
                   ret = rnorm(n_u * n_t, 0, 0.01))
fac <- data.frame(t = 1:n_t, mkt = rnorm(n_t, 0, 0.008))
long$ret <- long$ret + 0.9 * fac$mkt[long$t]
ev <- data.frame(unit = sprintf("u%02d", 1:20),
                 event_time = c(31:40, 31:40) + rep(c(0, 30), each = 10))
tau <- 0.02
win <- c(0, 4)
in_win <- function(d, e) {
  hit <- rep(FALSE, nrow(d))
  for (i in seq_len(nrow(e))) {
    o <- d$t - e$event_time[i]
    hit <- hit | (d$id == e$unit[i] & o >= win[1] & o <= win[2])
  }
  hit
}
long$ret[in_win(long, ev)] <- long$ret[in_win(long, ev)] + tau

ct <- feventr::calendar_time(long, "id", "t", "ret", events = ev,
                             window = win, factors = fac, returns = "simple")

# reference: portfolio by hand, then lm()
member <- in_win(long, ev)
pr <- tapply(long$ret[member], long$t[member], mean)
ref <- lm(pr ~ mkt, data = data.frame(pr = as.vector(pr),
                                      mkt = fac$mkt[as.numeric(names(pr))]))
expect_equal(ct$alpha, unname(coef(ref)[1]), tolerance = 1e-12)
expect_equal(ct$coefficients$estimate[2], unname(coef(ref)[2]), tolerance = 1e-12)
expect_equal(ct$alpha_se, summary(ref)$coefficients[1, 2], tolerance = 1e-12)
expect_equal(ct$nobs, length(pr))

# alpha recovers the per-period effect
expect_true(abs(ct$alpha - tau) < 3 * ct$alpha_se)

# intercept-only model = mean portfolio return
ct0 <- feventr::calendar_time(long, "id", "t", "ret", events = ev,
                              window = win, returns = "simple")
expect_equal(ct0$alpha, mean(pr), tolerance = 1e-12)

# overlapping events for the same unit enter once per period
ev_dup <- rbind(ev, ev[1, ])
ct_dup <- feventr::calendar_time(long, "id", "t", "ret", events = ev_dup,
                                 window = win, factors = fac, returns = "simple")
expect_equal(ct_dup$alpha, ct$alpha, tolerance = 1e-12)

# value weighting: put all weight on one member unit
long$mcap <- ifelse(long$id == "u01", 1, 1e-12)
ct_vw <- feventr::calendar_time(long, "id", "t", "ret", events = ev[1, , drop = FALSE],
                                window = win, weight = "mcap", returns = "simple")
u1 <- long[long$id == "u01" & long$t %in% (ev$event_time[1] + win[1]:win[2]), ]
expect_equal(ct_vw$portfolio$ret, u1$ret[order(u1$t)], tolerance = 1e-9)

# min_units drops thin ramp-up/ramp-down periods (peak concurrency is 5:
# events arrive one per period and each unit is held 5 periods)
ct_min <- feventr::calendar_time(long, "id", "t", "ret", events = ev,
                                 window = win, factors = fac,
                                 returns = "simple", min_units = 4)
expect_true(all(ct_min$portfolio$n_units >= 4))
expect_equal(ct_min$nobs + ct_min$diagnostics$n_dropped_periods, ct$nobs)

# Newey-West with lag 0 = HC0, and equals classical up to the df factor
ct_nw <- feventr::calendar_time(long, "id", "t", "ret", events = ev,
                                window = win, factors = fac,
                                returns = "simple", se = "nw", lag = 0)
expect_equal(ct_nw$alpha, ct$alpha, tolerance = 1e-12)
expect_true(abs(ct_nw$alpha_se / ct$alpha_se - 1) < 0.25)

# align = "value" on a pre-gapped event-time index (positions would shift)
gap <- long[long$t != 35, ]
ct_val <- feventr::calendar_time(gap, "id", "t", "ret", events = ev,
                                 window = win, factors = fac,
                                 returns = "simple", align = "value")
expect_true(!35 %in% ct_val$portfolio$time)
expect_true(abs(ct_val$alpha - tau) < 4 * ct_val$alpha_se)

# methods
expect_stdout(print(ct), pattern = "calendar-time portfolio")
expect_equal(unname(coef(ct)["alpha"]), ct$alpha)
expect_equal(nrow(summary(ct)), 2L)
