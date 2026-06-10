# PLAN.md — feventr: Causal Inference for Financial Event Studies

Living plan for building the `feventr` R package and reproducing Tables 1–7 of
Goldsmith-Pinkham & Lyu, *Causal Inference in Financial Event Studies*
(`~/Downloads/Event_Studies_and_Diff_in_Diff.pdf`). Replication data lives in
`~/Dropbox/PEAD_DinD` (strictly read-only). Approved implementation plan:
`~/.claude/plans/iterative-mapping-zebra.md`.

## Phase status

| Phase | Status | Gate |
|---|---|---|
| 0 — Recon artifacts | **DONE 2026-06-10** | Targets transcribed + adversarially verified 7/7 (0 corrections); all 5 discrepancies resolved below |
| 1 — Package core | **DONE 2026-06-10** (`da4a6d9`) | factor CARs ≡ lm() at 1e-10; DiD ≡ feols TWFE; `simulate_events` **bit-exact** vs published sim panel (seed 1234); Table 1 smoke shows the selection-bias pattern; 103 tests OK |
| 2 — Fast synthetic engine | **DONE 2026-06-10** | hybrid/qp ≤1e-4 of CVXR (achieved ~1e-7); SC ≡ augsynth, SDID ≡ synthdid (τ̂ to 1e-6, ω/λ to 1e-4); ridge closed form verified; benchmark grid + real-data support check in `benchmarks/`; 138 tests OK |
| 3 — Replication T1–T7 | pending | ±0.001 / ±0.1pp points; SEs ±20%; verifier sign-off per table |
| 4 — Polish | pending | snapshots, vignettes, R CMD check, CI |

## Decisions

- gsynth: wrap CRAN `gsynth` 1.2.1 for v1. SDID: own FW solver (synthdid math, dual BSD-3/GPL≥2).
- License GPL-3; provenance headers cite synthdid and augsynth (pinned `982f650b926de089922b8501de07ec04aa0cd0f5`).
- API: `event_study()` / `event_study_batch()` with `method=`; pure-matrix `eng_*` engines; 5 exports.
- Toolchain verified: R 4.4.2; data.table 1.16.4, osqp 0.6.3.3, gsynth 1.2.1, synthdid 0.0.9,
  augsynth 0.2.0, CVXR 1.8.1, tinytest 1.4.1, tinysnapshot 0.2.0, haven, fixest.
- CAR conventions needed (from recon): arithmetic **sum**, **compound** (∏(1+r)−1 diff), and
  **log** (Σ[log(1+r₁)−log(1+r₀)] — Table 6 uses this). `cumulate ∈ {sum, compound, log}`.

## Solver design (Phase 2 outcome)

`solve_simplex_ls(A, b, V, method = c("hybrid", "fw", "qp"))`:
- **hybrid** (default): Frank-Wolfe with exact line search (synthdid port, O(n0·t0)/iter,
  capped at 300 iters with a relative objective-decrease stop) seeds the donor support
  (only donors pushed above the shrunken-uniform baseline), OSQP polishes on that support,
  then **KKT gradient screening** (inactive donors must have half-gradient ≥ the active-set
  level; one O(n0·t0) pass per round, ≤10 rounds) re-admits any donors FW missed and
  re-polishes. Result: objective ≤ full-OSQP at eps 1e-8 (rel diff ~1e-7, often negative),
  21–121× faster on the simulated grid at n0 ∈ {500..5,000} (full grid in
  `benchmarks/solver_benchmark_results.csv`).
- True SC solutions are sparse (40–170 active donors at n0 ≤ 10,000 in both simulated and
  real Geithner all-CRSP data), which is why support restriction works.
- eng_sdid ports synthdid's collapsed-form algorithm exactly (zeta regularization,
  demeaned FW, sparsify-then-refit, min.decrease stopping): τ̂ matches synthdid to 1e-6.

## Targets

`replication/targets/table{1..7}.csv`, schema
`table,panel,row,col,estimate,se,stars,n,units,page,notes`. Each transcribed from the PDF and
independently re-verified cell-by-cell (verifier verdicts: 7/7 verified, 0 corrections).
Known paper typos preserved verbatim: Table 4's second "Panel B" (should be C); Table 2 Panel A
Gsynth N=129,625 is a transposition of 129,165 (cannot be produced by any sample; see Geithner notes).

## Table → script → data map (verified ground truth)

### Table 1 — Simulations (two-factor selection DGP)
- **Source**: `code/simulations_selection_TL.R` (engine) + `code/make_tables_selection_TL.R` (table).
  NOT `code/functions.r`/`simulations.r` (older 100-PC DGP) and NOT `simulations_new.R` (post-publication).
- **DGP**: T=500 days (239 est + day 240 = event + 10 post reported as 11-day window), n=500 firms.
  ret = β_mkt·MktRF + β_smb·SMB + rf + ε; β ~ iid N(1, 0.3²) each; **ε ~ N(0, 0.01152636²)**
  (= sqrt(var(MktRF)+var(SMB)+2cov) on full 1926–2022 FF daily sample; paper's "0.1" is loose —
  binding rule: residual variance = systematic variance). Factors: FF daily rows (Mkt-RF, SMB, RF)
  sampled jointly **without replacement in random order** ("block" = whole-day rows, not contiguous).
  Treatment +0.03 on event day only. Assignment: Bernoulli(0.1), or logit(δ·β_smb), δ=log(0.1)/mean(β_smb).
  Timing: day 240, or **deterministic argmax of SMB** over candidate days 240..489 (the logit rbinom
  draw is dead code but consumes RNG — replay for bit-exact match). Seeds: `seq(1234, 1283)`, one per sim.
- **Estimators**: simple means (feols ~ -1+i(t)+i(t,treated), HC1); CAPM (pre-period exret~MKT);
  correct factor (raw ret~MKT+SMB, abnormal = ret−rf−fitted — note rf asymmetry); gsynth(ret~D,
  force="unit", r=c(1,40), CV=TRUE, se=TRUE, inference="parametric", nboots=200 default, no seed).
- **Metrics** (make_tables lines 73–153): All-periods E(Bias)=mean_sims(Σ_t bias)/11×100,
  MAD=mean(|Σ_t bias|)/11×100, **RMSE=sqrt(mean_sims(Σ_t bias²))×100 (no /11)**; treated period
  t=event: mean bias×100 + coverage; untreated /10. **"Coverage" = rejection rate of H0 (p<.05)**:
  power on event day, size on post days.

### Tables 2–3 — Geithner (single event, 2008-11-21)
- **Scripts**: `Geithner/code/Geithner_sdid.do` (Panel A cols 1,2,6,7 + Table 3 bank betas),
  `Geithner_sdid_allcrsp.do` (Panel B + Table 3 all-CRSP), `Geithner_estimate_tl.do` (cols 3–5
  both panels), `Geithner_estimate_tl.R` / `Geithner_estimate_allcrsp_tl.R` (col 8 gsynth).
  Published numbers live in `output/geithner-connections.txt` / `geithner-full-sample.txt`.
- **Data**: `Geithner/Data and Programs/Stata Files/cleaned_data_beforeest.dta` (545 banks × 267 d),
  `cleaned_data_beforesdid.dta` (545 × 237; 129,165 obs), `cleaned_data_allcrsp_beforeest.dta`
  (3,883 × 266), `cleaned_data_beforesdid_allcrsp.dta` (4,095 × 236; 966,420 obs).
- **Estimand**: ATT = **average daily return over days 0–10** (firm-day pooled), NOT a CAR.
  Day 0 = Nov 21 3–4pm TAQ return (market = .06876978). Treated: 15 schedule / 9 personal / 38 NY.
- **Windows**: betas (CAPM/FF3F, Tables 2 cols 3–5 + Table 3) on dif ∈ [−280,−31], ≥225 obs.
  Synthetic estimators on balanced pre **[−256,−31] banks / [−255,−31] all-CRSP**; placebo window
  [−30,−1] deleted everywhere; post [0,10].
- **Calls**: Stata `sdid ret index_ds dif treated, method(sc|sdid|did) vce(placebo) reps(10) seed(123)`
  (paper says 100 reps; reps(100) agrees to 3 decimals). gsynth: `force="none", CV=TRUE, r=c(0,5),
  se=TRUE, inference="parametric", nboots=1000` (no seed), placebo days filtered, est.avg reported.
  Cols 3–5 = **one-sample t on treated abnormal returns only** (identical across panels by construction).
- **Table 3**: Panel B weighted betas = sdid `e(omega)` unit weights, `asgen` weighted means.
- **Gotchas**: Panel B day-0 inconsistency — cols 1,2,6,7 use full-day CRSP day-0 return; cols 3–5,8
  use 3–4pm TAQ. Panel B has 4,095 firms for cols 1,2,6,7 vs 3,883 for col 8.

### Tables 4–5 — Index inclusion (792 firm-events, 635 anndates, 1976–2023)
- **Scripts**: curated pipeline `index_inclusion/code/gp_lyu_indexinc_replication/` steps 1–7.
  Table 4 ← `clean_index_dates_siblis.do` → `output/average_beta_index_inclusion{,_randcon}_siblis.csv`.
  Table 5 ← `index_include_carplots_siblis.do` "Table 5" block → `output/att_index_inclusion_siblis.xlsx`
  (hand-pasted collapse output; verified to the digit).
- **Data**: `include_event_date_siblis.dta` (events), `panel_ii_{1_272,...,500_635}.dta` (cohort
  panels, event_date −280..+20, donors shrcd 10/11/12 + exchcd 1–3, complete 301-day window,
  missing daret zero-filled), `sc_ii_siblis.dta` (gsynth Y.tr/Y.ct), beta .dta files,
  `permno_anndate_random_controls.dta` (**unseeded runiform draw — must consume saved file**).
- **Conventions**: day 0 = anndate moved to most recent PRIOR trading day; pre-1989-09 anndate :=
  effdate−1. **Table 5 effect measured at event_date == +1** (announcements after close). Betas on
  [−250,−101] (listwise missing). gsynth: pre [−280,−101], treated = include·1(event_date≥−100),
  `force="unit", r=c(1,10), se=FALSE`, per cohort; covers 613/635 cohorts (21 failed + #635 loop bug).
  ar_capm/ar_ff3f **subtract estimated alpha**; ar_sp = daret − sprtrn; ar_mean = treated−donor means.
- **No inference anywhere** in Tables 4–5; decade groups by year(anndate): 80–89, 90–99, 00–09, 10–20.

### Table 6 — M&A acquirer 3-day CARs (14,847 deals, 6,625 dates)
- **Scripts**: `M&A/code/sdc_ma_malmendier_gsynth.do` lines ~400–573 → `output/car_log_sc_summ_stats.tex`
  (verified exact). Upstream: `clean_sdc_deals.do`, `sdc_ma_malmendier_gsynth_batch.do/R` (panels +
  per-date gsynth), `sdc_ma_malmendier_canonical.do` (market-adjusted).
- **Data to consume**: `M&A/output/sl_m_deals_car_1_250_gsynth.dta` (deal-level, 14,847 rows) and/or
  `data/work/sdc_ma_sl_m_gsynth/sc_ma_1_7052.dta` (stacked Y.tr/Y.ct per deal-day) +
  `output/sl_m_deals_car_1_250_mkt.dta`. **Raw event panels (event_panel_i.dta) are NOT in Dropbox**
  (lived on RA machine) — rebuild via batch.do only if needed.
- **Conventions**: daily; day 0 = first trading day ≥ SDC dateann. Panels [−280,+250]. gsynth:
  treated = treat·1(event_date≥−30) → pre [−280,−31]; `force="unit", r=c(1,10), se=FALSE`; missing
  daret zero-filled. **Published cells = LOG CARs**: Σ_{[−1,+1]} [log(1+r_treat) − log(1+r_cf)];
  market-adjusted vs CRSP vwretd. Paper's "−250 to −100" loading window is **wrong vs code** ([−280,−31]).
  Deal filters: 100% cash or 100% stock, deal/acq mkt val ≥5%, completed, pheld<50, pctacq≥50.
- **No SEs in pipeline** (se=FALSE; table = estpost summarize). Cash(9,261)+stock(5,592)=14,853>14,847
  (6 deals double-coded).

### Table 7 — Close-contest winner/loser betas
- **Scripts**: `M&A/code/sdc_ma_malmendier_close_contest_betas.do` (pre_length=35 loop) +
  `sdc_ma_malmendier_close_contest_geo_catt.do` lines 146–229 (keep long_dur==1, Welch ttests,
  collapse mean/p50). Table hand-assembled in Excel — no tex export exists.
- **Data**: `M&A/data/MergerDataPGP.csv` (231 bidders, 111 contests),
  `data/work/ma_close_contest_permno_{capm,ff3f}_beta_prelen_12_35_wide.dta` (*35 columns),
  `close_contest_treatedret_adjt1_monthly.dta`.
- **Conventions**: **MONTHLY** (paper note's "daily" is a typo — proven by exact recomputation of all
  12 cells incl. HML Welch p=0.0497 → the single star). Betas: per-firm OLS of monthly excess
  delisting-adjusted simple returns (missing filled with vwretd) on FF3F_monthly (decimal), exactly
  35 obs, event months [−35,−1]. Sample: long_dur = FightLength ≥ median(175d) → 56 contests,
  119 firms (56 W / 63 L). Welch two-sided t-tests; fully deterministic, no seeds.
- KPSS/, BaGel/, PEAD earnings: **confirmed no Table 1–7 dependency**.

## Discrepancy resolutions (Phase 0 gate)

1. **Table 6 daily vs monthly** — RESOLVED: daily pipeline exists (`ma_sl_cleaned_2023_event_date.dta`
   is a daily trading-calendar map, not monthly); monthly panels belong to close contests.
2. **Table 1 DGP** — RESOLVED: `simulations_selection_TL.R` two-factor DGP (spec above); 100-PC
   `functions.r` DGP is the older, unpublished variant.
3. **Geithner windows** — RESOLVED: betas [−280,−31]; synthetic [−256,−31] (banks) / [−255,−31]
   (all-CRSP); placebo gap [−30,−1] excluded; post [0,10].
4. **SC matching basis** — RESOLVED: every synthetic method in every application matches on
   **per-period simple returns** (levels), never cumulated paths. `match_on="cumret"` stays as an
   API option but defaults to "ret".

## Replication-fidelity notes (for Phase 3 tolerance calls)

- Unseeded randomness upstream: Table 4 random controls (consume saved .dta), Geithner gsynth
  nboots=1000 (SE tolerance ±20% covers), Table 1 gsynth bootstrap (coverage cells).
- Hardcoded RA paths in .do files — replication layer reads Dropbox artifacts, never re-runs Stata.
- Where paper text and code disagree, **code wins** (documented above); cite in replication README.

## Verifier verdicts

- Phase 0 transcription: T1–T7 verified 7/7, 0 corrections (workflow `wf_306c0bf4-90f`, 2026-06-10).
- Phase 3 per-table verdicts: pending.
