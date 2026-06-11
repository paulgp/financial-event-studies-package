## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6, fig.height = 4)

## -----------------------------------------------------------------------------
library(feventr)
sim <- simulate_events(n_units = 200, n_pre = 120, n_candidate = 1,
                       n_post = 10, tau = 0.03, seed = 42)
str(sim$events)

## -----------------------------------------------------------------------------
args <- list(data = sim$data, unit = "id", time = "t", ret = "ret",
             treated = sim$events$unit, event_time = sim$event_time,
             window = c(0, 10), est_window = c(-120, -1), returns = "simple")

fit_mean <- do.call(event_study, c(args, list(method = "mean")))
fit_sc   <- do.call(event_study, c(args, list(method = "sc", se = "placebo",
                                              reps = 50, seed = 1)))
fit_sc

## -----------------------------------------------------------------------------
head(summary(fit_sc))

## -----------------------------------------------------------------------------
fit_2f <- do.call(event_study, c(args, list(
  method = "factor", factors = sim$factors[, c("t", "mktrf", "smb")])))
round(fit_2f$att, 4)

## -----------------------------------------------------------------------------
plot(fit_sc)
plot(fit_sc, what = "paths")

## -----------------------------------------------------------------------------
event_betas(fit_sc, sim$factors[, c("t", "mktrf", "smb")])

## -----------------------------------------------------------------------------
ev <- data.frame(unit = c("1", "2", "10"), event_time = c(80, 80, 100),
                 event = c("a", "a", "b"))
b <- event_study_batch(sim$data, "id", "t", "ret", events = ev,
                       method = "sc", window = c(0, 5),
                       est_window = c(-60, -1), returns = "simple")
b
summary(b)

## -----------------------------------------------------------------------------
A <- matrix(rnorm(50 * 200), 50, 200)   # 200 donors, 50 pre-periods
w_true <- c(rep(0.25, 4), rep(0, 196))
b_ <- as.vector(A %*% w_true) + rnorm(50, sd = 0.05)
sol <- solve_simplex_ls(A, b_)
round(head(sort(sol$w, decreasing = TRUE), 6), 3)

