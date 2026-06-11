# Plot-method snapshot tests (tinysnapshot; SVG snapshots in _tinysnapshot/).
# Snapshots are generated on macOS; skip where font rendering differs.
if (!requireNamespace("tinysnapshot", quietly = TRUE) ||
    !requireNamespace("svglite", quietly = TRUE))
  exit_file("tinysnapshot/svglite not available")
if (nzchar(Sys.getenv("CI")) || Sys.info()[["sysname"]] != "Darwin")
  exit_file("snapshots rendered on macOS only")

# using() (not library()) so tinytest collects the snapshot expectations
using(tinysnapshot)
# the guards above already restrict to local macOS runs; don't let
# tinysnapshot's NOT_CRAN heuristic silently skip inside the test runner
options(tinysnapshot_plot_skip = FALSE, tinysnapshot_plot_review = FALSE)

sim <- feventr::simulate_events(n_units = 60, n_pre = 50, n_candidate = 1,
                                n_post = 6, tau = 0.03, seed = 77)
args <- list(data = sim$data, unit = "id", time = "t", ret = "ret",
             treated = sim$events$unit, event_time = sim$event_time,
             window = c(0, 6), est_window = c(-50, -1), returns = "simple")

f_mean <- do.call(feventr::event_study, c(args, list(method = "mean")))
f_sc <- do.call(feventr::event_study, c(args, list(method = "sc", se = "placebo",
                                                   reps = 25, seed = 4)))

expect_snapshot_plot(function() plot(f_mean), label = "fit_att_mean")

# conformal band: stored CI bounds, deterministic (no seed needed)
f_conf <- do.call(feventr::event_study, c(args, list(method = "sc",
                                                     se = "conformal")))
expect_snapshot_plot(function() plot(f_conf), label = "fit_att_sc_conformal")
expect_snapshot_plot(function() plot(f_sc, what = "car"), label = "fit_car_sc")
expect_snapshot_plot(function() plot(f_sc, what = "paths"), label = "fit_paths_sc")
expect_snapshot_plot(function() plot(f_sc, what = "weights"), label = "fit_weights_sc")

# batch plot
ev <- data.frame(unit = as.character(c(1, 2, 5, 6)),
                 event_time = c(60, 60, 70, 70), event = c("a", "a", "b", "b"))
set.seed(9)
long <- data.frame(id = rep(1:40, times = 100), t = rep(1:100, each = 40),
                   ret = rnorm(4000, 0, 0.01))
long$ret[long$id %in% c(1, 2) & long$t == 60] <-
  long$ret[long$id %in% c(1, 2) & long$t == 60] + 0.04
long$ret[long$id %in% c(5, 6) & long$t == 70] <-
  long$ret[long$id %in% c(5, 6) & long$t == 70] + 0.04
b <- feventr::event_study_batch(long, "id", "t", "ret", events = ev,
                                method = "mean", window = c(0, 4),
                                est_window = c(-30, -1), returns = "simple")
expect_snapshot_plot(function() plot(b), label = "batch_att_mean")

# calendar-time portfolio plots (reuse the staggered batch panel)
ct <- feventr::calendar_time(long, "id", "t", "ret", events = ev,
                             window = c(0, 4), returns = "simple")
expect_snapshot_plot(function() plot(ct), label = "caltime_car")
expect_snapshot_plot(function() plot(ct, what = "n_units"),
                     label = "caltime_n_units")
