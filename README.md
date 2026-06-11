# feventr — Causal Inference for Financial Event Studies

`feventr` is a lean R package implementing the estimators studied in
Goldsmith-Pinkham & Lyu, *Causal Inference in Financial Event Studies*, and the
replication code that reproduces the paper's Tables 1–7 with it. Every
estimator consumes the same long panel (unit × date × return, plus treatment
assignment and event dates) and returns a common fit object: ATT path by event
time, cumulative effects, standard errors, diagnostics, and plot methods.

## Estimators

| `method =` | Estimator | Inference (`se = "auto"`) |
|---|---|---|
| `mean` / `did` | Difference in means / difference-in-differences | two-sample t |
| `market` | Market-adjusted abnormal returns (β pinned at 1) | one-sample t |
| `factor` | Factor-model abnormal returns (CAPM, FF3F, arbitrary user factors) | one-sample t |
| `sc` | Synthetic control, fast hybrid solver (below) | placebo (100 reps) |
| `ridge` | Ridge-augmented SC (augsynth-style closed form, LOO λ CV) | placebo |
| `sdid` | Synthetic difference-in-differences (synthdid algorithm port) | placebo |
| `gsynth` | Generalized synthetic control (wraps CRAN `gsynth`) | parametric bootstrap (1,000) |

Two first-class modes: `event_study()` (one event, many treated units) and
`event_study_batch()` (hundreds–thousands of events, each with its own event
date and donor pool; fit per event in parallel, cross-event SEs). A third,
complementary design via `calendar_time()`: the Jaffe–Mandelker calendar-time
portfolio estimator (Fama 1998) — each calendar period's portfolio holds every
unit within the event window, portfolio returns are regressed on factors, and
the alpha is the per-period abnormal return (classical or Newey–West SEs;
equal- or value-weighted). Beta diagnostics via `event_betas()`; the paper's
two-factor simulation DGP via `simulate_events()` (reproduces the published
simulation panels bit-for-bit given the same seeds).

## Quickstart

```r
library(feventr)
sim <- simulate_events(seed = 1234)              # 500 firms, FF-factor DGP
fit <- event_study(sim$data, unit = "id", time = "t", ret = "ret",
                   treated = sim$events$unit, event_time = sim$event_time,
                   method = "sc", window = c(0, 10), est_window = c(-239, -1),
                   returns = "simple")
summary(fit)
plot(fit)                                        # ATT path with placebo CIs
plot(fit, what = "paths")                        # treated vs synthetic cumret
event_betas(fit, sim$factors[, c("t", "mktrf", "smb")])
```

## The fast synthetic-control solver

Synthetic control at financial-panel scale (thousands of donor stocks) is the
binding constraint: total replication cost is per-fit time × n_events ×
(1 + inference reps). Standard implementations form the dense n₀ × n₀ donor
Gram matrix and hand it to a QP solver — quadratic-plus cost in the donor
count. `solve_simplex_ls()` instead runs Frank-Wolfe with exact line search
(the synthdid algorithm; O(n₀·t₀) per iteration, no Gram matrix) to seed the
sparse active donor set, polishes with an exact QP restricted to that support,
and verifies KKT optimality with a full gradient screen, re-admitting any
missed donors. The result matches the full QP at eps 1e-8 to ~1e-7 relative
(often better) at a small fraction of the cost:

| n₀ donors | t₀ = 100: full QP | hybrid | speedup | t₀ = 250: full QP | hybrid | speedup |
|---:|---:|---:|---:|---:|---:|---:|
| 500 | 0.45s | 0.02s | 21× | 0.08s | 0.06s | 1× |
| 2,000 | 2.17s | 0.08s | 29× | 2.07s | 0.19s | 11× |
| 5,000 | 23.5s | 0.20s | 121× | 21.0s | 0.49s | 42× |
| 10,000 | 186s | 0.42s | 446× | 492s | 0.98s | 504× |

(simulated 3-factor donor panels; `benchmarks/solver_benchmark.R`. On real
data — the Geithner all-CRSP window, 4,080 donors × 225 days — the hybrid is
123–131× faster than the full QP at +4.4e-9 relative objective, insensitive
to the support-size cutoff; `benchmarks/support_sensitivity_real.R`.)

## Replicating the paper

`replication/` reproduces Tables 1–7 against the paper's data repository
(licensed CRSP/Compustat/SDC content — paths configured via
`FEVENTR_PEAD_DIND`, never shipped here). Published point estimates were
transcribed from the PDF into `replication/targets/` and every reproduced
table is compared cell-by-cell; see `PLAN.md` for per-table verification
verdicts and documented deviations.

## Installation

```r
# from this repository
install.packages(c("data.table", "osqp", "gsynth"))
devtools::install_local("r/feventr")
```

## License and provenance

GPL-3. The Frank-Wolfe solver ports the algorithm of the
[synthdid](https://github.com/synth-inference/synthdid) package (Arkhangelsky,
Athey, Hirshberg, Imbens & Wager; dual BSD-3/GPL≥2); ridge augmentation
follows [augsynth](https://github.com/ebenmichael/augsynth) (Ben-Michael,
Feller & Rothstein; comparison benchmarks pinned at commit `982f650b`).
Bundled `ff_daily` factors are from the Kenneth R. French data library.
