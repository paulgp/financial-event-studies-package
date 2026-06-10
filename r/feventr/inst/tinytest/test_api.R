# API smoke tests: fit objects, methods, inference fields

set.seed(3)
sim <- feventr::simulate_events(n_units = 50, n_pre = 40, n_candidate = 1,
                                n_post = 4, seed = 21)
args <- list(data = sim$data, unit = "id", time = "t", ret = "ret",
             treated = sim$events$unit, event_time = sim$event_time,
             window = c(0, 4), est_window = c(-40, -1), returns = "simple")

for (m in c("mean", "did", "market", "factor")) {
  a <- args
  a$method <- m
  if (m %in% c("market", "factor")) {
    a$factors <- if (m == "market") data.frame(t = sim$factors$t, mkt = sim$factors$mktrf)
                 else sim$factors[, c("t", "mktrf", "smb")]
  }
  f <- do.call(feventr::event_study, a)
  expect_inherits(f, "fes_fit")
  expect_equal(length(f$att), 5L)
  expect_equal(names(f$att), as.character(0:4))
  expect_equal(length(f$car), 5L)
  expect_true(is.finite(f$att_avg))
  # default tstat inference populated
  expect_equal(f$se$method, "tstat")
  expect_equal(length(f$se$att), 5L)
  expect_true(is.finite(f$att_avg_se))
  # methods
  expect_silent(s <- summary(f))
  expect_equal(nrow(s), 5L)
  expect_equal(coef(f), f$att)
  expect_equal(coef(f, cumulative = TRUE), f$car)
  expect_equal(dim(vcov(f)), c(5L, 5L))
  ci <- confint(f)
  expect_true(all(ci[, 1] <= f$att & f$att <= ci[, 2]))
  expect_stdout(print(f), pattern = m)
}

# se = "none" skips inference
f0 <- do.call(feventr::event_study, c(args, list(method = "mean", se = "none")))
expect_null(f0$se)
expect_error(vcov(f0), "no stored")

# event day lands at the simulated effect: mean estimator at t=0 near tau
f1 <- do.call(feventr::event_study, c(args, list(method = "mean")))
expect_true(abs(f1$att[["0"]] - sim$tau) < 0.02)

# keep_data = FALSE drops the panel
f2 <- do.call(feventr::event_study, c(args, list(method = "mean", keep_data = FALSE)))
expect_null(f2$panel)

# synthetic methods are Phase 2
expect_error(do.call(feventr::event_study, c(args, list(method = "sc"))),
             "Phase 2")
