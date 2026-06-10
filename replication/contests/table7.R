# Table 7 — Close merger contest winner/loser betas (p. 52)
#
# Sample: 111 close merger contests (MergerDataPGP.csv bidders), restricted to
# long-duration contests (FightLength >= 175-day median): 56 contests,
# 119 firms (56 winners, 63 losers).
#
# Conventions per PLAN.md (code wins over the paper note's "daily" typo):
#   - MONTHLY frequency. Betas = per-firm OLS of monthly excess
#     delisting-adjusted simple returns (daret - rf; missing bidder returns
#     filled with vwretd upstream) on FF3F_monthly factors (decimal units),
#     exactly 35 observations, event months -35..-1.
#     Calendar map: ym = mofd(dateann_earliest) + event_date - 1 (pre-period),
#     as in sdc_ma_malmendier_close_contest_betas.do.
#   - Table = group means + medians by winner status and Welch two-sided
#     t-tests on the difference in means (loser - winner), as in
#     sdc_ma_malmendier_close_contest_geo_catt.do lines 146-229.
#     Star on the difference if Welch p < 0.10 (HML is the only starred cell).
#   - Fully deterministic; no seeds.
#
# Fidelity checks performed here:
#   1. Per-firm betas are cross-checked against the saved inputs to the
#      published table, ma_close_contest_permno_{capm,ff3f}_beta_prelen_12_35
#      _wide.dta (*35 columns) — exact diffs reported (float-precision match).
#   2. feventr::event_betas() is exercised on one contest panel and verified
#      against the per-firm lm() betas (API demonstration; Table 7 itself is a
#      per-firm diagnostics table that pools across contests, so the table is
#      assembled from per-firm OLS).
#
# Run from replication/: Rscript contests/table7.R
suppressMessages({library(feventr); library(haven); library(data.table)})
source("config.R")

work <- function(f) dind("M&A/data/work", f)

# ---------------------------------------------------------------------------
# 1. Load the bidders-only monthly panel and the FF3F monthly factors
# ---------------------------------------------------------------------------
panel <- as.data.table(read_dta(work("close_contest_treatedret_adjt1_monthly.dta")))
ff <- as.data.table(read_dta(work("FF3F_monthly.dta")))   # ym (Stata %tm), decimals

mofd <- function(x) (as.integer(format(x, "%Y")) - 1960L) * 12L +
  as.integer(format(x, "%m")) - 1L

pre <- panel[event_date >= -35 & event_date <= -1]
# Pre-period calendar month: contest clock anchored at the earliest bid date
pre[, ym := mofd(dateann_earliest) + event_date - 1L]
pre <- merge(pre, ff, by = "ym", all.x = TRUE)
stopifnot(!anyNA(pre$daret_adj), !anyNA(pre$mktrf))
# daret_adj == raw delisting-adjusted daret (vwretd-filled) in the pre-period;
# the t=1 contest-period compounding adjustment never enters months -35..-1.
pre[, exret := daret_adj - rf]

# ---------------------------------------------------------------------------
# 2. Per-firm CAPM and FF3F betas, exactly 35 monthly obs each
# ---------------------------------------------------------------------------
betas <- pre[, {
  stopifnot(.N == 35L)
  bc <- coef(lm(exret ~ mktrf))
  b3 <- coef(lm(exret ~ mktrf + smb + hml))
  .(winner = winner[1], long_dur = long_dur[1], n = .N,
    bmkt_capm = bc[["mktrf"]], bmkt_ff3f = b3[["mktrf"]],
    bsmb_ff3f = b3[["smb"]],   bhml_ff3f = b3[["hml"]])
}, by = .(contestnumber, permno)]
stopifnot(nrow(betas) == 231L)

# Cross-check vs the saved wide beta files (the published numbers' inputs)
wide <- merge(
  as.data.table(read_dta(work("ma_close_contest_permno_capm_beta_prelen_12_35_wide.dta")))[
    , .(contestnumber, permno, bmkt_capm35)],
  as.data.table(read_dta(work("ma_close_contest_permno_ff3f_beta_prelen_12_35_wide.dta")))[
    , .(contestnumber, permno, bmkt_ff3f35, bsmb_ff3f35, bhml_ff3f35)],
  by = c("contestnumber", "permno"))
chk <- merge(betas, wide, by = c("contestnumber", "permno"))
stopifnot(nrow(chk) == 231L)
xchk <- chk[, .(
  capm_mkt = max(abs(bmkt_capm - bmkt_capm35)),
  ff3f_mkt = max(abs(bmkt_ff3f - bmkt_ff3f35)),
  ff3f_smb = max(abs(bsmb_ff3f - bsmb_ff3f35)),
  ff3f_hml = max(abs(bhml_ff3f - bhml_ff3f35)))]
cat("Max |beta diff| vs saved *_prelen_12_35_wide.dta (published inputs):\n")
print(xchk)
stopifnot(unlist(xchk) < 1e-5)   # Stata float storage precision

# ---------------------------------------------------------------------------
# 3. feventr::event_betas() demonstration on one long-duration contest
# ---------------------------------------------------------------------------
demo_contest <- betas[long_dur == 1][order(contestnumber)][1, contestnumber]
cp <- panel[contestnumber == demo_contest & event_date >= -35 & event_date <= 1]
cp[, ym := fifelse(event_date <= 0, mofd(dateann_earliest) + event_date - 1L,
                   mofd(eff_date) + event_date - 1L)]
cp <- merge(cp, ff, by = "ym", all.x = TRUE)
cp[, exret := daret_adj - rf]
cf <- unique(cp[, .(event_date, mktrf, smb, hml)])[order(event_date)]
fit <- event_study(cp, unit = "permno", time = "event_date", ret = "exret",
                   treated = as.character(cp[winner == 1, unique(permno)]),
                   event_time = 0, window = c(0, 1), est_window = c(-35, -1),
                   method = "mean", returns = "simple", se = "none")
eb <- event_betas(fit, factors = cf, time = "event_date")
# Single winner per contest => treated mean beta == that firm's lm() beta
eb_w <- as.data.frame(eb)[eb$group == "treated" & eb$stat == "mean", ]
ref <- betas[contestnumber == demo_contest & winner == 1]
demo_diff <- max(abs(c(eb_w$mktrf - ref$bmkt_ff3f, eb_w$smb - ref$bsmb_ff3f,
                       eb_w$hml - ref$bhml_ff3f)))
cat(sprintf("event_betas() demo, contest %s: max |beta - lm beta| = %.3g\n",
            demo_contest, demo_diff))
stopifnot(demo_diff < 1e-10)

# ---------------------------------------------------------------------------
# 4. Table 7: long-duration sample, means/medians by winner, Welch t-tests
# ---------------------------------------------------------------------------
ld <- betas[long_dur == 1]
cat(sprintf("Long-duration sample: %d contests, %d firms (%d winners, %d losers)\n",
            uniqueN(ld$contestnumber), nrow(ld), sum(ld$winner == 1),
            sum(ld$winner == 0)))
stopifnot(uniqueN(ld$contestnumber) == 56L, nrow(ld) == 119L,
          sum(ld$winner == 1) == 56L, sum(ld$winner == 0) == 63L)

# Stata `ttest ..., by(winner) welch`: diff = mean(winner==0) - mean(winner==1)
# = loser - winner; Welch (1947) degrees of freedom.
welch <- function(x0, x1) {                # x0 = losers, x1 = winners
  n0 <- length(x0); n1 <- length(x1)
  w0 <- var(x0) / n0; w1 <- var(x1) / n1
  tstat <- (mean(x0) - mean(x1)) / sqrt(w0 + w1)
  df <- (w0 + w1)^2 / (w0^2 / (n0 + 1) + w1^2 / (n1 + 1)) - 2
  list(diff = mean(x0) - mean(x1), t = tstat, df = df,
       p = 2 * pt(-abs(tstat), df))
}

vars <- c("CAPM Beta" = "bmkt_capm", "FF3F Mkt Beta" = "bmkt_ff3f",
          "FF3F SMB Beta" = "bsmb_ff3f", "FF3F HML Beta" = "bhml_ff3f")
rows <- rbindlist(lapply(names(vars), function(lab) {
  v <- vars[[lab]]
  x1 <- ld[winner == 1][[v]]; x0 <- ld[winner == 0][[v]]
  wt <- welch(x0, x1)
  cat(sprintf("  %-14s Welch t = %6.3f, df = %6.2f, p = %.4f\n",
              lab, wt$t, wt$df, wt$p))
  data.table(row = lab,
             col = c("Winner Mean", "Winner Median", "Loser Mean",
                     "Loser Median", "Mean t-test Loser - Winner"),
             estimate = c(mean(x1), median(x1), mean(x0), median(x0), wt$diff),
             stars = c("", "", "", "", if (wt$p < 0.10) "*" else ""))
}))

out <- data.table(table = 7L, panel = "", rows, units = "beta")
fwrite(out, out_path("table7.csv"))

# ---------------------------------------------------------------------------
# 5. Comparison vs targets
# ---------------------------------------------------------------------------
tgt <- as.data.table(read_target(7))
cmp <- merge(tgt[, .(row, col, target = estimate,
                     target_stars = fifelse(is.na(stars), "", stars))],
             out[, .(row, col, reproduced = estimate, stars)],
             by = c("row", "col"), sort = FALSE)
cmp[, diff := reproduced - target]
# published values rounded to 3 dp => +/-0.0005 rounding slack on the 0.001 tol
cmp[, ok := abs(diff) <= 0.0015 & target_stars == stars]
fwrite(cmp, out_path("table7_comparison.csv"))

cat(sprintf("\nTable 7: %d/%d cells within tolerance (max |diff| = %.5f)\n",
            sum(cmp$ok), nrow(cmp), max(abs(cmp$diff))))
print(as.data.frame(cmp), digits = 4)
if (!all(cmp$ok)) stop("Table 7 comparison failed")
