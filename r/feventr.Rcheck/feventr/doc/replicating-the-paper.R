## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", eval = FALSE)

## -----------------------------------------------------------------------------
# sim <- simulate_events(selection = "both", seed = 1234)  # Panel D, sim 1

## -----------------------------------------------------------------------------
# fit <- event_study(banks, "index_ds", "dif", "ret",
#                    treated = schedule_connected, event_time = 0,
#                    method = "sdid",                 # or mean/did/market/...
#                    window = c(0, 10), est_window = c(-256, -31),
#                    returns = "simple", se = "placebo", reps = 100, seed = 123)

## -----------------------------------------------------------------------------
# fit <- event_study(cohort_panel, "permno", "event_date", "daret",
#                    treated = included_permnos, event_time = 0,
#                    method = "gsynth", force = "unit", r = c(1, 10),
#                    window = c(-100, 20), est_window = c(-280, -101),
#                    returns = "simple", se = "none")
# fit$att[["1"]]   # the announcement-day effect

## -----------------------------------------------------------------------------
# event_study(..., method = "gsynth", window = c(-30, 250),
#             est_window = c(-280, -31), returns = "simple", cumulate = "log")

