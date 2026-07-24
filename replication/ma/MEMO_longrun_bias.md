# Memo: CAR conventions, placebo drift, and long-run bias in the M&A application

*2026-07-22 — feventr replication repo, `replication/ma/`. All numbers regenerate
from the scripts cited at the end; paper references are to the 2025-12-28 draft.*

## TL;DR

1. **Table 6 is in the log-CAR convention while the draft's stated estimand
   (Lemma 1) is the arithmetic ATT.** Under the stated estimand the 3-day cells
   are Market **1.23** / Gsynth **1.04**, not 0.8 / 0.7 — the announcement
   effect is ~50% larger. The arithmetic columns already exist in the saved
   output files (`car_arith_vwret_1`, `car_arith_sc_1`).
2. **A date-matched placebo (no deal, same pipeline) drifts to −30% additive
   CATT at +250 days under SC — but the *same placebo units* under gsynth
   drift +0.7%.** SC's long-run drift is design bias; gsynth's survives that
   test. The SC bias mechanism is *noise in measured prices*, not factor
   misspecification: it collapses 5× at decimalization, scales with unit
   volatility, and is present in pre-event days.
2a. **But a second confound — runup reversion (Ashenfelter dip) — takes out
   over half of gsynth's remaining long-run effect.** Placebos matched to
   their acquirers on estimation-window cumulative return revert −5.8%
   [−9.5, −2.4] by +250d under the gsynth design (random placebos: +0.7),
   because the unit FE carries a transient hot pre-mean forward; reversion
   is −12.7% for stock-deal-sized runups. The corrected long-run effect
   (treated − runup placebo): **−5.0% [−7.5, −2.4] full sample, −6.8 stock,
   −3.9 cash** — half the naive gsynth −10.8/−19.5/−5.5, still
   significantly negative, cash/stock ratio narrows from 3.5× to ~1.7×.
   This reconciles the designs: SC's placebo-calibrated estimate is −1.9
   [−4.0, +0.1], overlapping gsynth's corrected −5.0 — both point to a
   modest −2 to −5% one-year drift, far below the BHAR-literature
   magnitudes. The announcement effect is untouched by every correction:
   ~+1.1 [0.8, 1.4] under both designs, placebo-calibrated. The two-placebo
   decomposition (random = mechanical bias; runup-matched − random =
   selection reversion) is the general-purpose diagnostic; the close
   contest gets both implicitly because losers share runup and noise.
3. **The bias has a sign structure that explains Figure 8's fan**: the market
   counterfactual is biased positive (treated-side noise), SC negative
   (donor-side selection), gsynth mildest. At 252 days the market-adjusted CAR
   is −12.9% in logs but **+3.1%** in arithmetic — the sign of acquirer
   long-run "underperformance" vs the market depends on the cumulation
   convention alone.
4. **The close-contest (LaLonde) result gets a mechanical rationale**: contest
   losers match winners on noise variance, so winner-minus-loser differences
   the inflation out — exactly like treated-minus-placebo. Date-matched
   placebos generalize that calibration to any event study.

## 1. The three CAR conventions and what the code computes

`sdc_ma_malmendier_gsynth.do` builds, per deal-day, from gsynth fit on simple
returns (`daret ~ treated`, counterfactual `daret_sc` = Y.ct):

| variable | construction | 3-day mean | 252-day mean |
|---|---|---:|---:|
| `car_arith_sc` | Σ daret_treated − Σ daret_sc (additive CATT) | **1.044%** | **−10.8%** |
| `car_log_sc` | Σ log1p(daret_treated) − Σ log1p(daret_sc) | 0.658% | −23.5% |
| `car_sc` | exp(Σlog)−1 differences (BHAR) | 1.040% | −62.1% |

Market column (`sdc_ma_malmendier_canonical.do`): arithmetic **1.233%** /
log 0.825% at 3 days; arithmetic **+3.1%** / log −12.9% at 252 days.

- Table 6's published cells (0.8 / 0.7) are the **log** means. The browse
  comments at `sdc_ma_malmendier_gsynth.do:409-414` show why the log version
  was kept: the BHAR blows up at 252 days (deal-level range −84,451% to
  +1,824%; one synthetic compounded exp(6.74) ≈ 847×), and log tames both
  tails. But the log convention charges the noisy realized leg ~σ²/2 per day
  that it does not charge the smooth predicted leg — a 0.38pp haircut already
  at 3 days, growing linearly in horizon. At 3 days BHAR ≈ arithmetic (1.040
  vs 1.044), so compounding is irrelevant there; the entire published-vs-
  arithmetic gap is this Jensen term.
- The method comparison inside Table 6 is unaffected (the haircut sits on the
  common treated leg), but "difference less than 10 basis points" is 17-19bp
  for the full sample under either convention.

## 2. The placebo experiment and what it found

`ma_refit_full.R MA_REFIT_PLACEBO=1` fits one **date-matched placebo
non-acquirer** per deal (seeded, real acquirers excluded from every donor
pool) through the identical SC pipeline: 7,052 announcement-date panels,
complete-531-day-coverage universe, window −30..+250, est. window −280..−31.
Fitted units and donors come from the same complete-coverage pool, so
survivorship is symmetric by construction.

Pooled placebo path at +250 days (full sample): **−51.9% log / −30.8%
arithmetic**. The per-day arithmetic bias in `att`:

| slice | placebo mean daily att |
|---|---:|
| 1977–2000 | −19 bp/day (peak 1990s: −22) |
| 2001–2023 | −4.4 bp/day (2010s: −0.9) |
| pre-event days −30..−2 (era 1) | −22 bp/day |
| unit att-variance Q1 → Q4 (full sample) | −6 → −18 bp/day |
| unit att-variance Q1 → Q5 within 2001+ | −2.4 → −6.2 bp/day |

Against the placebo benchmark: the announcement effect is unambiguous
(treated outside the placebo 95% band at +1 **and** +21 in every subsample;
full-sample day +1: +0.83% treated vs −0.30% placebo, arithmetic), while
beyond a quarter the treated path is statistically indistinguishable from its
placebo (at +250: treated −32.1% vs placebo −30.2%). In the 2001+ subsample
the announcement effect is +1.0–1.1% under all three methods (placebo −0.1%)
and cash-merger SC lands exactly on its placebo (−10.9% both).

## 3a. The two-channel framing (the organizing frame for the long-run section)

Two distinct variance complications operate in long-run event studies — Jensen
terms at two different transformations. **The standing estimand convention is
arithmetic (additive CATT) precisely because it eliminates the unfixable
channel and leaves the fixable one.**

**Channel 1 — metric level (Lemma 1).** return → log(1+return) is concave:
E[log(1+R)] ≈ E[R] − Var(R)/2. Geometric/BHAR estimands need the
counterfactual matched on *return variance* as well as expected return —
impossible for a diversified portfolio against a single stock, and with no
theoretical guidance for constructing a variance-matched portfolio. Vanishes
entirely under arithmetic estimands (sums are linear). **No good fix; avoid
the channel by choosing the arithmetic estimand.**

**Channel 2 — measurement level.** price → return is convex in the price
level: with observed price = value·(1+ε), ε transitory,
E[measured simple return] ≈ E[true return] + Var(ε). Variance *in prices*
manufactures mean return with no change in value, sits inside each daily
mean, and so accumulates faithfully under addition — the arithmetic estimand
is exposed whenever the two legs differ in price-noise variance (bid-ask
bounce in volatile stocks). **Its sign is design-dependent** — it runs
against whichever leg is cleaner: treated stock vs noise-free index →
spurious positive drift (+3.1% at 252d, market arithmetic); noise-selected
synthetic vs random unit → spurious negative drift (−30% SC placebo). That
sign structure is Figure 8's fan. **Fixable while keeping the arithmetic
estimand**: the inflation is a stationary unit-level constant, so unit
intercepts/differencing (gsynth FE, SDID, DID), matched single-firm
benchmarks (control firms, contest losers, placebos), ABK noise-robust
weighting, or coarser frequency all remove it.

Caution for paper text: do NOT state the converse "logs escape channel 2"
unqualified. The telescoping immunity holds only for a buy-and-hold leg's own
cumulated log return; any counterfactual leg built from daily *simple*
returns (a rebalanced portfolio, a prediction that gets log1p'd) re-imports
the noise inflation each day. Logs trade channel 2 for channel 1 only in the
matched-single-firm case — where neither channel was a problem to begin with.

Empirical map: the ~21pp log-vs-arithmetic gap in the SC placebo at +250d is
channel 1; the −31pp surviving in arithmetic is channel 2 (plus its σ²/2
donor-selection cousin); the decimalization collapse dates channel 2's
microstructure core.

**Channel 2 has two members with different frequency scaling** (established
by the monthly 3-year build, `ma_monthly_longrun.R`): (i) *transitory-noise
inflation* Var(ε) is billed once per price touch, so monthly evaluation
cuts it ~21× — the member the decimalization test dates; (ii) the
*real-volatility convexity wedge* σ²_period/2 scales with the period
length (σ²_month ≈ 21·σ²_day), so per unit of calendar time it is
**frequency-invariant** — switching to monthly data does not touch it.
Evidence: the monthly-SC random placebo still drifts −22.5% by +36m and
the hybrid (daily-fit weights, monthly evaluation) −20.6% — the hybrid
fixes only member (i), which the complete-coverage screen had already made
small — while monthly gsynth's placebo is clean (−1.7 [−6.9, +4.1] at
+36m), because a unit intercept kills BOTH members (both are stationary
unit-level mean wedges). Practical rule: frequency choice mitigates the
noise member only; variance-selected counterfactuals need an intercept (or
factor averaging) at ANY frequency.

**Channel 2 decomposes further (the treated-cross-section question).** Piece
A, *selection amplification*: with one treated firm per event and n₀ ≫ T₀,
the fit must track one noisy path and recruits high-variance donors to
interpolate its wiggles — this is what makes SC's bias (−30% placebo) far
exceed a naive EW benchmark's. A treated cross-section of N firms shrinks
the target's idiosyncratic noise like 1/√N and kills the selection channel;
gsynth is effectively a many-firm design even at N=1 (its factors average
thousands of donors), which is why its placebo is clean. Piece B, *level
mismatch*: noise inflation sits in each firm's MEAN, and means do not
diversify — a portfolio of N treated firms still carries its constituents'
average inflation (CMTW's EW index: thousands of firms, still ~6%/yr). Many
firms converts the bias from an unestimable firm-specific quantity into a
stable composition quantity, matchable on observables (size/price/vol/
liquidity) or removable with an intercept — but does not eliminate it.
Two corollaries: (i) pooling across many single-treated EVENTS is NOT
equivalent — the selection tilt is systematic across events, so it survives
pooling (hence the pooled placebo at −13bp/day, not zero); pooling buys
precision, not bias reduction. (ii) Calendar-time portfolios kill piece A by
construction but keep piece B when equal-weighted — ABK's gross-return
weighting targets exactly that residual. Full recipe: aggregate the treated
side (kills A), then match/weight/intercept away composition inflation
(kills B), subject to the runup-persistence caveat for intercepts.

## 3b. Mechanism detail: noise in measured prices, not factor misspecification

A measured daily simple return satisfies E[measured r] ≈ geometric drift +
½·(real return variance) + (transitory-noise variance). The last term is
Blume–Stambaugh (1983) generalized by Asparouhova–Bessembinder–Kalcheva (JF
2013): *any* mean-reverting component of prices — bid-ask bounce, price
pressure, stale quotes — inflates measured arithmetic means by its full
variance. Two properties make it a counterfactual problem:

- **Linearity preserves it and diversification does not remove it.** A
  daily-rebalanced portfolio's expected return is the weighted average of its
  constituents' inflated means. (Logs are immune in cumulation — the noise
  telescopes to the endpoints.) This is the "Caveat Compounder" warning of
  Canina, Michaely, Thaler & Womack (JF 1998) — the daily CRSP EW index is an
  implicitly daily-rebalanced portfolio whose measured return is inflated
  ~6%/yr by bounce and nonsynchronous trading — transplanted to synthetic
  counterfactuals: Σw·r is exactly such a portfolio, and the SC fit is worse
  than the EW index (−19bp/day ≈ 48%/yr pre-2001 vs their ~2bp/day) because
  weight optimization *selects into* the bias rather than merely averaging
  over it.
- **SC selects on variance.** With ~4,000 donors and 250 pre-days, the
  simplex loads on donors whose realized noise co-moved with the unit's —
  disproportionately volatile, wide-spread names. A random unit has average
  noise; its synthetic has above-average noise; `att` drifts negative every
  day, in and out of sample.

The fingerprints above (decimalization step, volatility gradient, presence in
pre-event days, symmetry of survivorship) identify this channel and exclude
factor misspecification, which has no reason to halve at 2001Q2. Post-2001
the literal bounce is too small; the residual −4bp/day is the same formula
with the remaining variance sources (transitory components + the real-vol
σ²/2 wedge) — the within-2001+ volatility gradient confirms the variance link.

**Sign structure by counterfactual** (explains Figure 8's fan ordering):

| counterfactual | noise exposure | long-run bias sign |
|---|---|---|
| market index | treated-side only (index clean) | **positive** (arith +3.1% at 252d) |
| gsynth | factors average 1000s of donors | mild negative (−10.8%) |
| SC (simplex, n₀≫T₀) | donor-side, selected | **large negative** (placebo −30%) |
| APM | noisy counterfactual path | large negative (−35.7%) |

## 4. How this maps into the draft

1. **Table 6**: report `car_arith_*` (already saved) or state the log
   convention explicitly; reconciles the table with Lemma 1's "we focus on
   estimating the arithmetic ATT." Cells become 1.2 / 1.0; adjust the
   "<10bp" sentence (gap is ~19bp full-sample, still economically small).
2. **Lemma 1 companion (measurement)**: BHAR needs variance-matched
   counterfactuals; at daily frequency even the *arithmetic* ATT needs
   noise-variance-matched counterfactuals. Bias = per-day noise-variance
   differential × horizon; era-dependent (spreads), estimator-dependent
   (selection). Cites: Blume & Stambaugh (JFE 1983); Asparouhova,
   Bessembinder & Kalcheva (JFE 2010; JF 2013, "Noisy Prices and Inference
   Regarding Returns"); Canina, Michaely, Thaler & Womack (JF 1998),
   "Caveat Compounder: A Warning about Using the Daily CRSP Equal-Weighted
   Index to Compute Long-Run Excess Returns",
   doi:10.1111/0022-1082.165353 — the direct antecedent: our SC placebo
   drift is their daily-rebalancing bias with an optimized (and therefore
   noise-selected) portfolio in place of the EW index. The literature's
   response (Barber & Lyon JFE 1997 control firms; Lyon, Barber & Tsai JF
   1999 buy-and-hold benchmarks) anticipates the fix: a matched single-firm
   benchmark carries the same inflation and differences it out — the
   control-firm logic our placebo/close-contest/intercept results restate
   in causal-inference terms.
3. **"When does bias accumulate" section**: this channel is an order of
   magnitude larger at daily frequency (10–50%/yr for SC-style fits pre-2001)
   than the misspecification illustration (~1.8% at 3 years); worth a
   paragraph + the placebo as the diagnostic. Figure 7 already shows the
   symptom: SC pre-period counterfactual 0.19pp/day vs treated 0.14pp/day.
4. **Figure 8**: add the placebo path/band (or a 2001+ panel); the fan
   ordering is predicted by the sign table above, not only by factor
   misspecification.
5. **Close contest**: one sentence on *why* it works — losers match winners
   on noise variance, so the inflation differences out; date-matched placebos
   generalize the calibration when no contest exists.
6. **Inference**: deals cluster on dates and windows overlap; announcement-
   time block bootstrap (18-month blocks) gives design effects 2–11× naive
   cross-deal SEs at horizons ≥ +21.

**Prototype outcomes** (`ma_fix_check.R`; placebo drift, full sample,
additive CATT, 95% CI at +250d; per-day bias pre-2001 / post):

| design | +250d placebo drift | bp/day pre-2001 / post |
|---|---|---|
| SC baseline (fixed weights, no fix) | −30.2 [−39.8, −19.0] | −18.6 / −4.4 |
| SC + intercept (est-window demeaning) | +2.0 [−5.0, +9.3] | +0.1 / +2.2 |
| SC + ABK weighting (donor leg purged) | −3.3 [−9.8, +1.7] | −1.5 / −1.1 |
| Gsynth (unit FE + factor averaging) | +0.7 [−1.4, +3.1] | +0.3 / +0.3 |

Both estimand-preserving fixes eliminate the bias — a one-line demeaning or
a one-day gross-return reweighting each removes ~95% of a 30pp/yr drift.
One surprise refines the mechanism: ABK was predicted to land *positive*
(it purges only the donor leg, leaving the unit's own inflation), but it
lands at ≈0 — so a random complete-coverage stock's own noise inflation is
small (~1–2bp/day; the 531-day-presence screen drops the extreme-noise
tail), and the baseline SC bias was almost entirely **piece-A selection**
of high-noise donors rather than the level of typical-stock noise. The
placebo-validated menu for long-horizon SC at daily frequency is
therefore: intercept, ABK weighting, or switch to a factor-averaged
counterfactual — all three now demonstrated equivalent on the null.

## 5. Repro pointers (this repo)

| artifact | what |
|---|---|
| `ma/ma_refit_full.R` | 7,052-panel refit runner; `MA_REFIT_PLACEBO=1` for placebo fits |
| `ma/ma_refit_longrun.R` | pooled paths, both metrics, block bootstrap → `output/ma_longrun_paths.csv`, `output/ma_longrun.png` |
| `ma/ma_placebo_check.R` | SC placebo diagnostics → `output/ma_placebo_check.{csv,png}` |
| `ma/ma_placebo_gsynth.R` | gsynth placebo (same seeded units) vs published treated → `output/ma_placebo_gsynth.{csv,png}` |
| `ma/ma_runup_check.R` | runup-matched vs random gsynth placebos; reversion + corrected effect with CIs → `output/ma_runup_check.{csv,png}` |
| `ma/ma_placebo_calibrated.R` | placebo-calibrated CATT (treated − placebo, replicate-level CIs) → `output/ma_placebo_calibrated.{csv,png}` |
| `ma_refit_full.R` knobs | `MA_REFIT_PLACEBO=runup` (runup-matched placebos), `MA_REFIT_DEMEAN=1` (SC-with-intercept), `MA_REFIT_ABK=1` (prior-gross-return weighted counterfactual) |
| `ma/ma_fix_check.R` | four-placebo comparison of the fixes → `output/ma_fix_check.{csv,png}` |
| `ma/wrds_monthly_pull.R` | CRSP monthly 1974–2024 from WRDS (delist-adjusted) → Dropbox work dir |
| `ma/ma_refit_monthly.R` | monthly-frequency runner, −36..+36 months, 488 announcement-month cohorts, same knobs |
| `ma/ma_refit_hybrid.R` | daily-fit SC weights, monthly evaluation; placebo knobs |
| `ma/ma_monthly_longrun.R` | 3-year results, 9 columns → `output/ma_monthly_longrun.{csv,png}` |

**Three-year (literature-horizon) results, full sample, additive CATT [95% CI]:**
under the cleanest design (monthly gsynth), the raw treated drift at +36m is
−25.7 [−35.9, −15.5], of which the runup-matched placebo reverts −21.7
[−31.1, −12.9] — so the **corrected 3-year effect is −4.1 [−9.3, +2.0],
statistically zero** at exactly the horizon where Rau–Vermaelen report −4%
(significant) and Loughran–Vijh −25% for stock mergers. Monthly-SC corrected
is marginally negative (−5.8 [−10.4, −1.4]); the hybrid's corrected numbers
are wide and unstable (its own placebo drift is large, so the correction
differences two big numbers). The announcement quarter is +1.4 [0.4, 2.5]
under every design. Reading for the draft: at 3 years the classic long-run
merger underperformance decomposes as ≈0 mechanical (with an intercept
design) + ≈−22pp runup reversion + ≈−4pp residual deal effect (n.s.);
MMP's winner-loser −24% remains the honest upper bound among
quasi-experimental designs (contested-deal subsample).
| `ma/ma_longrun_placebo_fig.R` | combined figure → `output/ma_longrun_placebo.png` |
| `ma/ma_longrun_placebo_2001.R` | post-decimalization variant → `output/ma_longrun_placebo_2001.{csv,png}` |
| commits | `556d986` (bootstrap + placebo), `45e56d2` (additive CATT switch), `4c7f50a` (2001+ figure) |

Stata provenance: `M&A/code/sdc_ma_malmendier_gsynth.do` (~lines 127–137 log
pass, 300–313 BHAR pass, 400–461 all-three pass + save), market column in
`sdc_ma_malmendier_canonical.do`. gsynth is fit on simple returns
(`daret ~ treated`); `daret_sc` is a smooth prediction of the simple return,
which is why summing `log1p` of it re-introduces a Jensen term despite the
additive (CATT) intent.
