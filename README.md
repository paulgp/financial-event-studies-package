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
| `apm` | Aggregated projection matrix: spectral counterfactual outcome means under general missingness (Lei & Ross, arXiv:2312.07520; GitHub [`apm`](https://github.com/brad-ross/apm)) | weighted bootstrap (200) |

Beyond `se = "auto"`: `se = "conformal"` (methods mean/did/sc/ridge/sdid) runs
Chernozhukov–Wüthrich–Zhu refit-under-the-null inference with *exact*
permutation distributions — per-period p-values and CIs via single-post
enumeration and a moving-block joint test for a constant effect — so it is
deterministic (no Monte Carlo, no seed) and fast: each null refit warm-starts
the solver from the previous solution (~4s for a full set of CIs at 2,000
donors). Placebo inference parallelizes with `cores =` (assignments are
pre-drawn, so results are identical for any core count).

Two first-class modes: `event_study()` (one event, many treated units) and
`event_study_batch()` (hundreds–thousands of events, each with its own event
date and donor pool; fit per event in parallel with per-event checkpointing,
cross-event SEs with optional event weights and per-event SE propagation —
see `?event_study_batch`). A third,
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

## Method timings

Seconds **per event** on simulated staggered panels (500 donors, one treated
unit per cohort, window `c(0, 10)`, single core), by estimation-window length
t₀. The first two columns are point estimates only — `event_study_batch()`
with its default per-event `se = "none"`, the batch design. The last column
is a single-event `event_study()` fit at the `se = "auto"` defaults (t-stat
for mean/did/market/factor, placebo with 100 reps for sc/ridge/sdid,
parametric bootstrap with 1,000 draws for gsynth, weighted bootstrap with
200 draws for apm):

| `method =` | estimation only, t₀=100 | estimation only, t₀=250 | with `se = "auto"`, t₀=250 |
|---|---:|---:|---:|
| `mean` | 0.02s | 0.03s | 0.02s |
| `did` | 0.02s | 0.04s | 0.03s |
| `market` | 0.02s | 0.03s | 0.03s |
| `factor` | 0.02s | 0.04s | 0.02s |
| `sc` | 0.04s | 0.09s | 5.5s |
| `ridge` | 0.05s | 0.15s | 6.5s |
| `sdid` | 1.43s | 3.62s | 330s |
| `gsynth` | 0.12s | 0.50s | 109s |
| `apm` | 0.09s | 0.27s | 25s |

(`benchmarks/method_benchmark.R`; the results CSV also carries a 10-event
grid showing per-event cost is flat in the number of cohorts — the panel
layer trims each event to its own windows before copying anything.)
Inference is free where it is closed-form (the t-stat methods) and costs
roughly reps × the refit time where it is resampled; placebo refits
parallelize with `cores =`, and conformal inference (sc-family) is a cheaper
deterministic alternative (~4s at 2,000 donors, above). Batch runs
parallelize linearly with `cores =` (gsynth derated to `cores %/% 3` workers
by default; see `?event_study_batch` for checkpointing, event weights, and
per-event SE propagation). Costs grow with the donor pool: on the S&P 500
index-inclusion application (~4,000–6,000 donors, 301-day windows, 635
cohorts) gsynth runs at ~1.8s/event and apm at ~3.0 core-seconds per event
end-to-end.

## Replicating the paper

`replication/` reproduces Tables 1–7 against the paper's data repository
(licensed CRSP/Compustat/SDC content — paths configured via
`FEVENTR_PEAD_DIND`, never shipped here). Published point estimates were
transcribed from the PDF into `replication/targets/` and every reproduced
table is compared cell-by-cell; see `PLAN.md` for per-table verification
verdicts and documented deviations.

### The M&A example end-to-end

The paper's largest application — Table 6's 14,847 acquirer deals — also
runs end-to-end through the package rather than consuming the paper's saved
fits. `replication/ma/ma_refit_full.R` rebuilds all 7,052 announcement-date
event panels from CRSP daily (4,000–4,300 complete-coverage donors × 531
trading days each; forked workers share the keyed 4 GB panel copy-on-write,
so there are no per-worker copies and no intermediate panel cache) and fits
every deal separately — the treated acquirer against a donor pool excluding
all same-day acquirers, `window = c(-30, 250)`, `est_window = c(-280, -31)`,
log CARs — checkpointed per announcement date and method (~75 minutes for
`sc` at 6 cores; the same runner produces the `apm` column).
`ma_refit_compare.R` then reproduces the Table 6 cells: the sc and apm
refits match the paper's saved per-deal gsynth CARs at deal-level
correlations 0.96 and 0.84, and both reproduce every subsample cell's sign
and ordering. `ma_refit_longrun.R` extends the horizon to +250 trading
days, reporting additive CATTs in simple returns — cumulated daily `att`,
the paper's estimand. Acquirers drift negative after the announcement-day
pop under every counterfactual, but how negative depends heavily on the
counterfactual: −11% (gsynth), −32% (sc), −36% (apm) at +250 days, with
100%-stock deals worse than 100%-cash under each
(`replication/output/ma_longrun.png`).

Three inference lessons come with the long horizon. First, deals cluster
on announcement dates and 250-day windows overlap across deals announced
within a year of each other, so naive cross-deal SEs are fiction: the
default bands are a circular block bootstrap over announcement time
(18-month blocks of announcement months, longer than the event window),
and the implied design effects are 2–11× the naive SEs at horizons
beyond a month. Second, the metric matters: Table 6's cells cumulate
log1p of realized and predicted *simple* returns, a convention that adds
a per-day Jensen term whenever the counterfactual is smoother than the
realization (log1p of a low-variance prediction loses none of the ~σ²/2
per day a noisy realization does). Shared across columns, it leaves the
paper's method comparison intact, but it lowers levels even at 3 days
(the published full-sample gsynth cell is 0.66; additive and
buy-and-hold agree on ~1.04) and compounds at +250 days — on the paper's
own saved gsynth fits, −23% log versus −11% additive — so the long-run
tables and figure here use the additive CATT (the CSV carries both
metrics). Third, even
the additive CATT is not bias-free at long horizons:
`ma_placebo_check.R` re-runs the identical pipeline on one date-matched
placebo non-acquirer per deal, and the SC placebo path drifts to −30% by
+250 days with the bid-ask-bounce fingerprint (Blume–Stambaugh):
−19bp/day before decimalization versus −4bp/day after 2001, increasing
in the unit's volatility, and already present in the out-of-sample
pre-announcement days (−30..−2), where no treatment exists — the
synthetic's measured simple returns are inflated because the fit loads
on volatile, wide-spread donors. Measured against the placebo benchmark,
the announcement effect is unambiguous — treated is outside the placebo
95% band at +1 and +21 in every subsample — while beyond a quarter the
treated path is statistically indistinguishable from its placebo
(`replication/output/ma_placebo_check.png`).

## Installation

```r
# from this repository
install.packages(c("data.table", "osqp", "gsynth"))
devtools::install_local("r/feventr")

# optional: method = "apm" (GitHub-only engine)
remotes::install_github("brad-ross/apm", subdir = "r")
```

## License and provenance

GPL-3. The Frank-Wolfe solver ports the algorithm of the
[synthdid](https://github.com/synth-inference/synthdid) package (Arkhangelsky,
Athey, Hirshberg, Imbens & Wager; dual BSD-3/GPL≥2); ridge augmentation
follows [augsynth](https://github.com/ebenmichael/augsynth) (Ben-Michael,
Feller & Rothstein; comparison benchmarks pinned at commit `982f650b`).
Bundled `ff_daily` factors are from the Kenneth R. French data library.
