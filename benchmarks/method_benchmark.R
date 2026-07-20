# Benchmark: every estimator on simulated staggered panels, two stages.
#
# Stage 1 — estimation only: event_study_batch() over a grid of cohort
# counts x estimation-window lengths (per-event se = "none", cores = 1), the
# honest per-fit cost of the point estimates. Per-event cost is flat in the
# cohort count (the panel layer trims each event to its own windows before
# copying), so the README table reports the K = 50 column per t0.
#
# Stage 2 — with default inference: single-event event_study() fits at
# se = "auto" (t-stat for mean/did/market/factor, placebo 100 reps for
# sc/ridge/sdid, parametric bootstrap 1,000 for gsynth, weighted bootstrap
# 200 for apm), t0 = 250, same 500-donor DGP.
#
# DGP: 500 donors + one treated unit per cohort, two-factor returns, events
# spaced 12 trading periods apart, window c(0, 10).
# Output: benchmarks/method_benchmark_results.csv + a markdown table for the
# README on stdout.
suppressMessages(library(feventr))

n_donors <- 500L
methods <- c("mean", "did", "market", "factor", "sc", "ridge", "sdid",
             "gsynth", "apm")
if (!requireNamespace("gsynth", quietly = TRUE))
  methods <- setdiff(methods, "gsynth")
if (!requireNamespace("apm", quietly = TRUE))
  methods <- setdiff(methods, "apm")

mk_panel <- function(K, t0, seed) {
  set.seed(seed)
  spacing <- 12L
  d1 <- t0 + 11L
  ev_dates <- d1 + spacing * (0:(K - 1L))
  Tn <- ev_dates[K] + 10L
  N <- n_donors + K
  Fm <- cbind(mkt = rnorm(Tn, 3e-4, 0.01), smb = rnorm(Tn, 0, 0.006))
  L <- cbind(rnorm(N, 1, 0.3), rnorm(N, 0.3, 0.4))
  ret <- L %*% t(Fm) + matrix(rnorm(N * Tn, 0, 0.008), N, Tn)
  for (i in seq_len(K)) {
    tr <- n_donors + i
    ret[tr, ev_dates[i] + 0:10] <- ret[tr, ev_dates[i] + 0:10] + 0.01
  }
  list(long = data.frame(id = rep(seq_len(N), times = Tn),
                         t = rep(seq_len(Tn), each = N),
                         ret = as.vector(ret)),
       events = data.frame(unit = as.character(n_donors + seq_len(K)),
                           event_time = ev_dates),
       factors = data.frame(t = seq_len(Tn), mkt = Fm[, 1], smb = Fm[, 2]))
}

pick_factors <- function(m, p) {
  if (m == "market") p$factors[, c("t", "mkt")]
  else if (m == "factor") p$factors
  else NULL
}

# ---- stage 1: batch point estimates over the grid ---------------------------
grid <- expand.grid(K = c(10L, 50L), t0 = c(100L, 250L))
rows <- list()
for (i in seq_len(nrow(grid))) {
  K <- grid$K[i]; t0 <- grid$t0[i]
  p <- mk_panel(K, t0, seed = 1000 + i)
  for (m in methods) {
    tt <- system.time(
      b <- event_study_batch(p$long, "id", "t", "ret", events = p$events,
                             method = m, window = c(0, 10),
                             est_window = c(-(t0 + 10L), -11L),
                             returns = "simple", factors = pick_factors(m, p),
                             se = "cross", cores = 1L)
    )[["elapsed"]]
    ok <- sum(b$events$status == "ok")
    rows[[length(rows) + 1L]] <- data.frame(
      stage = "batch_point", method = m, n_events = K, t0 = t0,
      n_donors = n_donors, ok = ok, total_sec = round(tt, 2),
      sec_per_event = round(tt / K, 3))
    cat(sprintf("K=%d t0=%d %-7s %6.1fs (%.3fs/event, %d/%d ok)\n",
                K, t0, m, tt, tt / K, ok, K))
  }
}

# ---- stage 2: single-event fits with default inference ----------------------
`%||%` <- function(a, b) if (is.null(a)) b else a
p <- mk_panel(10L, 250L, seed = 1003)
ev1 <- p$events[1, ]
for (m in methods) {
  tt <- system.time(
    f <- event_study(p$long, "id", "t", "ret", treated = ev1$unit,
                     event_time = ev1$event_time, method = m,
                     window = c(0, 10), est_window = c(-260, -11),
                     returns = "simple", factors = pick_factors(m, p),
                     donors = as.character(seq_len(n_donors)),
                     se = "auto", seed = 1, keep_data = FALSE)
  )[["elapsed"]]
  rows[[length(rows) + 1L]] <- data.frame(
    stage = "single_event_se_auto", method = m, n_events = 1L, t0 = 250L,
    n_donors = n_donors, ok = 1L, total_sec = round(tt, 2),
    sec_per_event = round(tt, 3))
  cat(sprintf("se=auto t0=250 %-7s %7.1fs (%s)\n", m, tt,
              f$se$method %||% "none"))
}

res <- do.call(rbind, rows)
write.csv(res, "benchmarks/method_benchmark_results.csv", row.names = FALSE)

# markdown table for the README
val <- function(m, stage, K, t0v)
  res$sec_per_event[res$stage == stage & res$method == m &
                      res$n_events == K & res$t0 == t0v]
fmt <- function(v) if (v >= 10) sprintf("%.0fs", v) else sprintf("%.2fs", v)
cat("\n| `method =` | estimation only, t0=100 | estimation only, t0=250 |",
    "with `se = \"auto\"`, t0=250 |\n|---|---:|---:|---:|\n")
for (m in methods)
  cat(sprintf("| `%s` | %s | %s | %s |\n", m,
              fmt(val(m, "batch_point", 50L, 100L)),
              fmt(val(m, "batch_point", 50L, 250L)),
              fmt(val(m, "single_event_se_auto", 1L, 250L))))
