# Benchmark: every estimator through event_study_batch() on simulated
# staggered panels, over a grid of cohort counts x estimation-window lengths.
# Single core (cores = 1), se = "none" per event, so the numbers are the
# honest per-fit cost; batch mode parallelizes linearly with `cores` (gsynth
# derated to cores %/% 3 by default, see ?event_study_batch).
#
# DGP: 500 donors + one treated unit per cohort, two-factor returns, events
# spaced 12 trading periods apart, window c(0, 10), est_window of length t0.
# Output: benchmarks/method_benchmark_results.csv + a markdown table for the
# README on stdout.
suppressMessages(library(feventr))

n_donors <- 500L
methods <- c("mean", "did", "market", "factor", "sc", "ridge", "sdid",
             "gsynth", "cfm", "apm")
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

grid <- expand.grid(K = c(10L, 50L), t0 = c(100L, 250L))
rows <- list()
for (i in seq_len(nrow(grid))) {
  K <- grid$K[i]; t0 <- grid$t0[i]
  p <- mk_panel(K, t0, seed = 1000 + i)
  for (m in methods) {
    fx <- if (m == "market") p$factors[, c("t", "mkt")]
          else if (m == "factor") p$factors else NULL
    tt <- system.time(
      b <- event_study_batch(p$long, "id", "t", "ret", events = p$events,
                             method = m, window = c(0, 10),
                             est_window = c(-(t0 + 10L), -11L),
                             returns = "simple", factors = fx,
                             se = "cross", cores = 1L)
    )[["elapsed"]]
    ok <- sum(b$events$status == "ok")
    rows[[length(rows) + 1L]] <- data.frame(
      method = m, n_events = K, t0 = t0, n_donors = n_donors,
      ok = ok, total_sec = round(tt, 2),
      sec_per_event = round(tt / K, 3))
    cat(sprintf("K=%d t0=%d %-7s %6.1fs (%.3fs/event, %d/%d ok)\n",
                K, t0, m, tt, tt / K, ok, K))
  }
}
res <- do.call(rbind, rows)
write.csv(res, "benchmarks/method_benchmark_results.csv", row.names = FALSE)

# markdown table for the README: sec/event by config
cfg <- unique(res[, c("n_events", "t0")])
cat("\n| `method =` |",
    paste(sprintf("%d events, t0=%d", cfg$n_events, cfg$t0), collapse = " | "),
    "|\n|---|", paste(rep("---:", nrow(cfg)), collapse = "|"), "|\n", sep = "")
for (m in methods) {
  v <- vapply(seq_len(nrow(cfg)), function(j)
    res$sec_per_event[res$method == m & res$n_events == cfg$n_events[j] &
                        res$t0 == cfg$t0[j]], 0)
  cat(sprintf("| `%s` | %s |\n", m,
              paste(sprintf("%.2fs", v), collapse = " | ")))
}
