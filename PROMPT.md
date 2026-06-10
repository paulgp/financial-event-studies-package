# Prompt: Build `feventr` — Causal Inference for Financial Event Studies

## Mission

Build an R package implementing the estimators from Goldsmith-Pinkham & Lyu,
*"Causal Inference in Financial Event Studies"*
(`/Users/psg24/Downloads/Event_Studies_and_Diff_in_Diff.pdf`), and use the package to
reproduce the paper's simulation and empirical tables (Tables 1–7). A Python port
is deferred — out of scope for now; just keep the estimator core free of
R-idiosyncratic state (plain matrices in, fit objects out) so a port stays feasible. Work through the
plan below autonomously, phase by phase, with the verification gates described in
"Validation protocol." Begin by entering plan mode and producing a concrete
implementation plan; once approved, execute it end to end.

## Design philosophy: minimum viable

This is a lean, research-grade package, not enterprise software. The model is
`fixest` (github.com/lrberge/fixest @ `1cb65a4`): a compact codebase with a minimal
dependency footprint and a single pragmatic test file (`tests/fixest_tests.R`, with
a small custom helper in `R/test_fun.R`) — not hundreds of unit tests. Concretely:

- **Implement only what Tables 1–7 need.** No speculative options, estimators, or
  abstraction layers; prefer a well-factored function over a class hierarchy. When
  in doubt, leave it out — features can be added when a table demands them.
- **Testing is gate-driven, not coverage-driven.** The test suite *is* the
  validation protocol: reference-solver agreement, table-target regression checks,
  plus a handful of API smoke tests. A few focused test files, not one per
  function. Do not chase coverage metrics.
- **Testing framework: `tinytest`, with `tinysnapshot` for plot methods.** Visual
  output is a core feature of an event-study package, so each fit object's plot
  method gets an SVG snapshot test with visual diffs — layout as in
  grantmcdermott/tinyplot: test files in `inst/tinytest/`, snapshots in
  `inst/tinytest/_tinysnapshot/`.
- **Few dependencies.** base/data.table + osqp, Rcpp only if profiling demands it;
  wrap CRAN `gsynth`/`synthdid` rather than reimplementing where allowed.

## Context and resources (read these first)

1. **The paper**: `/Users/psg24/Downloads/Event_Studies_and_Diff_in_Diff.pdf` (91 pp).
   Read the Simulations section and the four Empirical Examples in the Applications
   section. The tables to reproduce:
   - **Table 1** — Treatment-effect bias and coverage in simulations (two-factor
     structure with event-timing selection).
   - **Tables 2–3** — Geithner Treasury Secretary announcement (Acemoglu et al. 2016):
     ATT over an 11-day post window with estimators {difference in average, DiD,
     market-adjusted, CAPM, FF3F, synthetic control, synthetic DiD, gsynth}, two
     control samples (banks; all CRSP common stocks); plus treated/control beta table.
   - **Tables 4–5** — S&P 500 index inclusion (Greenwood & Sammon 2025 sample): beta
     distributions of included firms by decade; announcement-day treatment effects by
     decade with estimators {diff-in-means, market, CAPM, FF3F, gsynth}.
   - **Table 6** — M&A acquirer 3-day CARs by target type and payment method
     (market-adjusted mean vs. gsynth mean).
   - **Table 7** — Close merger contests: beta distributions of winners vs. losers.

2. **Existing replication code and data**: `~/Dropbox/PEAD_DinD` (read its `CLAUDE.md`).
   This is the paper's working repo: Stata + R, application subdirectories
   (`Geithner/`, `index_inclusion/`, `M&A/`, `KPSS/`, `BaGel/`, `simulations/`),
   `code/functions.r` (simulation DGP: FF betas + 100-PC factor structure from CRSP
   daily returns), and `data/` (~5.4GB raw CRSP daily returns CSV, ~26GB processed
   `.dta` files, FF factor files). The existing scripts are the ground truth for what
   each table computes — read the relevant scripts before implementing each
   application. The cleaned event-window datasets in each application directory are
   the inputs the replication layer should consume.

3. **Optimizer reference**: `~/repos/claude-container/projects/synthdid-package/synthdid`.
   Its `R/solver.R` implements a Frank-Wolfe solver with exact line search
   (`fw.step()`, `sc.weight.fw()`, `sc.weight.fw.covariates()`) for simplex-constrained
   least squares — the approach that made synthdid fast where general QP solvers
   (CVXR/ECOS, OSQP) were slow or fragile. `R/reference-solver.R`
   (`simplex.least.squares()` via CVXR) shows how to validate a fast solver against a
   reference. Use this as the template for the fast synthetic-control engine.

## What the package must provide

A coherent API where every estimator consumes the same panel input
(unit × date × return, plus treatment assignment and event dates, supporting both a
single clustered event and many staggered events) and returns a common fit object
(ATT path by event time, cumulative effects, SEs, diagnostics, plot methods).

Two distinct estimation modes, both first-class:

- **Single-event mode** (Geithner): one event date, many treated units, one fit.
- **Many-event batch mode** (index inclusion, M&A): hundreds–thousands of events,
  each with its own event date, treated unit(s), and donor pool; fit per event,
  then average ATTs across events with cross-event SEs. The PEAD_DinD repo runs
  these as batch jobs in 1,000-event chunks — the package should make this a single
  parallelized call. Note the runtime arithmetic: total cost = per-fit time ×
  n_events × (1 + inference reps); this, not any single fit, is why solver speed
  is binding.

Estimator inventory (all needed for the tables):

| Estimator | Notes |
|---|---|
| Diff-in-means / DiD | Trivial baselines, but same API |
| Factor-model abnormal returns | Market-adjusted, CAPM, FF3F (arbitrary user-supplied factors); pre-event estimation window (e.g. -250 to -101), project forward, CARs |
| Synthetic control | **Fast, scalable** — see below. Simplex weights matched on pre-event return paths |
| Augmented/ridge SC (augsynth-style) | Ridge augmentation has closed form given SC weights; choose lambda by CV as augsynth does (`R/ridge_lambda.R`) |
| Synthetic DiD | Port/wrap the Frank-Wolfe approach from the synthdid reference |
| Generalized synthetic control (gsynth/IFE) | Interactive fixed effects with cross-validated factor number; parametric bootstrap SEs |
| Inference | Placebo inference (Arkhangelsky et al. 2021 style; Table 2 uses 100 repetitions), parametric bootstrap (Table 2 gsynth: 1,000 draws), two-sample t; cross-event SEs for batch mode |
| Beta diagnostics | Treated-vs-control factor-loading tables (Tables 3, 4, 7) |
| Simulation DGP | Two-factor DGP with selection-on-timing from the paper's simulation section (Table 1: 50 simulations × {mean-diff, market, two-factor, gsynth}, reporting bias in pp, RMSE, 95% coverage), reproducible without CRSP data |

### Financial-panel details the implementation must get right

These conventions differ across applications and silently change the numbers —
extract each one from the paper and the corresponding PEAD_DinD script during
Phase 0 and record the choice per application in `PLAN.md`:

- **Return definition and cumulation**: simple vs. log returns, and whether CARs
  sum (log) or compound (simple). PEAD_DinD has both raw and log-return dataset
  variants.
- **What synthetic methods match on**: per-period returns vs. cumulated return
  paths in the pre-window. These give different weights; the paper's framing is
  "realized pre-event return paths" — confirm against the scripts.
- **Event time**: trading-day offsets, not calendar days; define day 0 per
  application (e.g. Geithner day 0 is a partial day: 3pm–close on 2008-11-21).
- **Donor eligibility**: complete pre-window return history required (balanced
  panel within each event window), share-code/exchange filters (the Stata cleaning
  scripts encode these: common shares, NYSE/AMEX/NASDAQ), exclusion of
  contemporaneously-treated units in batch mode.
- **Estimation vs. event windows** per application (e.g. index inclusion: loadings
  estimated on days -250 to -101 relative to announcement).
- **V / time-weighting matrix**: the solver API must accept a V matrix (augsynth
  compatibility), defaulting to identity.

### The key technical problem: synthetic control at scale

**Confirmed diagnosis** (benchmarked 2026-06-10 against the augsynth source,
github.com/ebenmichael/augsynth): augsynth's SCM solver (`synth_qp` in
`R/fit_synth.R`) forms the dense n0×n0 matrix `Pmat = X0 V X0'` over the donor pool
and hands it to OSQP at eps 1e-8 with n0+1 constraints. Cost is quadratic-plus in
the number of donors: on a 100-pre-period factor-structure problem, OSQP takes
~0.07s at n0=500, ~2.2s at n0=2,000, ~25s at n0=5,000 (and ~200MB for the dense P
at n0=5,000). The staggered `multisynth` path (`R/multi_synth_qp.R`) is worse — it
stacks one QP with n0×J variables across J treatment cohorts. The ridge
augmentation itself (`R/ridge.R`) is a cheap closed-form t0×t0 solve and is *not*
the bottleneck.

**Validated solver design** (prototype benchmarks on the same problems):

1. **Frank-Wolfe with exact line search** (synthdid `solver.R` pattern) — O(n0·t0)
   per iteration, never forms the n0×n0 Gram matrix. Linear scaling confirmed:
   2,000 iterations in ~1.1s at n0=5,000. Caveat: plain FW has a sublinear
   convergence tail, so its objective stalls slightly above the exact optimum.
2. **Support-restricted QP polish**: the exact SC solution is sparse (~50 active
   donors at n0=5,000 in the prototype). Run FW to identify the support, keep the
   top ~5·t0 donors by weight, then solve the exact OSQP problem restricted to that
   support. Prototype: objective within 1e-4 relative of the full-OSQP optimum in
   1.2s vs 27s — ~23× faster with effectively exact accuracy.

Implement the hybrid (FW → restricted QP) as the default solver, with pure-FW and
full-QP as options. Ridge augmentation in closed form on top. Core loops in
Rcpp/RcppArmadillo only if profiling demands it.

- Correctness gate: match a CVXR reference solver and `augsynth` itself to
  tolerance on problems small enough for both to finish.
- Speed gate: reproduce/extend the benchmark table (n0 ∈ {500, 2,000, 5,000,
  10,000} donors, t0 ∈ {100, 250}); confirm ≥10× at n0 ≥ 2,000 with objective
  within 1e-4 relative of the reference. Verify the support-restriction heuristic
  on real return data (not just simulated factor structure) — check sensitivity to
  the support-size cutoff.

## Repository layout

Work in `/Users/psg24/repos/financial-event-studies-package` (this repo). Suggested:

```
r/          # R package (e.g. `feventr`) — package code, tinytest suite, vignettes
replication/  # Scripts that call the package against ~/Dropbox/PEAD_DinD data
              # and emit Tables 1–7; one subdir per application; outputs to
              # replication/output/ as csv + tex
benchmarks/   # Optimizer benchmark scripts + results
PLAN.md       # The living plan: phases, status, decisions, deviations
```

**Hard constraint**: CRSP/Compustat/SDC data is licensed — none of it (raw or
derived unit-level) goes in the package or the git history. The package ships only
simulated example data from the built-in DGP. Replication scripts read from
`~/Dropbox/PEAD_DinD` via a configurable path and are expected to be runnable only
on machines with the data.

**Scope note**: the PEAD_DinD repo also contains `KPSS/`, `BaGel/`, and PEAD
earnings-announcement work — these are *not* among the paper's Tables 1–7 and are
out of scope unless Phase 0 recon finds a table depends on them.

**Dependencies and licensing**:
- `augsynth` is GitHub-only (`devtools::install_github("ebenmichael/augsynth")`) —
  pin the commit used for comparison benchmarks; `gsynth`, `synthdid`, `osqp`,
  `CVXR` are on CRAN.
- Code ported from synthdid (BSD-3/GPL≥2) and augsynth carries license obligations:
  pick a compatible license for the package (GPL-3 or MIT+attribution depending on
  what is ported vs. reimplemented from the math) and record provenance in the
  source headers.

## Phases

Develop the detailed plan yourself, but it must cover these phases in order, each
with its verification gate passed before moving on:

- **Phase 0 — Reconnaissance.** Read the paper's estimator definitions and table
  notes; read the PEAD_DinD scripts behind each table; read the synthdid solver.
  Write `PLAN.md` mapping each table → data inputs → estimators → expected numbers
  (transcribe the published point estimates from the PDF into a machine-readable
  `replication/targets/` file — these are the regression targets).
- **Phase 1 — Package core.** Panel data structures, factor-model estimators,
  diff baselines, simulation DGP. Gate: tinytest suite passes; factor-model CARs
  match a hand-rolled lm() computation; Table 1 simulation runs (small-scale smoke
  version).
- **Phase 2 — Fast synthetic engine.** FW solver, ridge augmentation, SDID,
  gsynth wrapper or IFE implementation (decide in the plan: wrapping the existing
  `gsynth` CRAN package is acceptable for v1; the fast SC engine is not optional).
  Gate: reference-solver agreement + benchmark table.
- **Phase 3 — Replication.** One application at a time: Geithner → index
  inclusion → M&A → close contests → Table 1 full-scale. Gate per application:
  reproduced table matches the PDF targets (see tolerances below).
- **Phase 4 — Polish.** Plot methods with tinysnapshot tests, two vignettes
  (quickstart on simulated data + "replicating the paper"), README with the
  benchmark table, R CMD check clean, GitHub Actions CI.

## Validation protocol

- **Tolerances**: point estimates within ±0.001 (returns in decimal) or ±0.1pp of
  published table values; SEs within 20% (placebo/bootstrap SEs are seed-dependent —
  fix seeds and document). Any cell outside tolerance gets investigated and either
  fixed or documented in `PLAN.md` with the reason (e.g. data vintage differences).
- **Never trust a single implementation**: every estimator gets at least one
  independent check (reference solver, existing CRAN package, or analytic
  result on a constructed example).
- After each phase, commit with a message describing what was verified. Do not
  rewrite history; do not push unless asked.

## Execution workflow

Run this as an orchestrated, multi-agent workflow — use the Workflow tool and
subagents; you are explicitly authorized to fan out agents wherever the structure
below calls for it.

- **Session shape**: enter plan mode for Phase 0 and present the plan once for
  approval. After approval, work through phases without asking for permission
  between steps. Stop and ask only when: (a) a replication target can't be matched
  after genuine investigation, (b) a scope decision changes the deliverable (e.g.
  IFE from scratch vs. wrapping gsynth turns out to matter for results), or
  (c) anything would require modifying files in `~/Dropbox/PEAD_DinD` (treat that
  directory as strictly read-only).
- **Fan out where work is independent**:
  - Phase 0: parallel reader agents — one per empirical application — each mapping
    its PEAD_DinD scripts/data to its paper table and returning structured notes;
    a separate agent transcribing table targets from the PDF.
  - Phase 3: the four applications are independent given the package — run each as
    its own agent (worktree isolation if they touch shared files), each ending with
    a target-vs-output comparison.
- **Verify adversarially**: after each replication table is produced, spawn an
  independent verifier agent that compares output against `replication/targets/`
  and tries to refute the match (units, sample, window, seed). A table is "done"
  only when the verifier signs off; record the verdict in `PLAN.md`.
- **Long-running jobs** (full Table 1 simulations, large-donor benchmarks,
  placebo/bootstrap inference) run as background tasks; downscale first to validate
  the pipeline, then scale up while continuing other work.
- **Checkpointing**: keep `PLAN.md` current — phases done, decisions, deviations,
  verifier verdicts. It is the resume point if the session restarts; on restart,
  re-read it and continue rather than re-deriving.
