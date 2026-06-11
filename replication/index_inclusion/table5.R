# Table 5 — Announcement-day effect of S&P 500 inclusion by decade
# (Goldsmith-Pinkham & Lyu, p. 44; units = percent). Effect measured at
# event_date == +1 (announcements come after the close; PLAN.md).
#
# Original: index_inclusion/code/gp_lyu_indexinc_replication/
# 6_index_include_carplots_siblis.do "Table 5" block ->
# output/att_index_inclusion_siblis.xlsx. Columns:
#  * Diff-in-Means: treated minus non-treated mean daret per announcement+1
#    date over all eligible CRSP firms (shrcd 10/11/12, exchcd 1-3), then
#    decade mean over dates (decade by year of the +1 date).
#  * Market: daret - sprtrn (CRSP dsp500), decade mean over firm-events.
#  * CAPM / FF3F: daret - rf - alpha - beta'factors, loadings from the
#    [-250, -101] per-event OLS (alpha subtracted, per PLAN.md).
#  * Gsynth: mean of (Y.tr - Y.ct) at +1. PROVISIONAL here: taken from the
#    saved sc_ii_siblis.dta artifact until table5_gsynth_full.R re-estimates
#    all 635 cohorts with feventr and overwrites the Gsynth rows.
# Decades by year(anndate) except Diff-in-Means (year of +1 date, as in the
# original code).
#
# Run from replication/: Rscript index_inclusion/table5.R

suppressMessages({library(haven); library(data.table)})
source("config.R")
source("index_inclusion/betas_common.R")

inp <- load_ii_inputs()
tw <- load_treated_returns(inp)

betas_f <- feventr_cache("ii_treated_betas.rds")
tr_betas <- if (file.exists(betas_f)) readRDS(betas_f) else
  compute_treated_betas(inp, tw)

dec_lab <- c("1980-1989", "1990-1999", "2000-2009", "2010-2020")

# ---- treated firm-event day +1 abnormal returns -----------------------------
sp <- as.data.table(read_dta(ii_work("sp500_sprtrn.dta")))
ev1 <- tw[event_date == 1]
ev1 <- merge(ev1, sp, by = "date", all.x = TRUE)
ev1 <- merge(ev1, inp$ff, by = "date", all.x = TRUE)
ev1 <- merge(ev1, tr_betas[, .(index, alpha_capm, beta_capm, alpha_ff3f,
                               bmkt_ff3f, bsmb_ff3f, bhml_ff3f)],
             by = "index", all.x = TRUE)
ev1[, ar_sp := daret - sprtrn]
ev1[, ar_capm := daret - rf - alpha_capm - beta_capm * mktrf]
ev1[, ar_ff3f := daret - rf - alpha_ff3f - bmkt_ff3f * mktrf -
      bsmb_ff3f * smb - bhml_ff3f * hml]
ev1[, group := decade_group(year(anndate))]
cat(sprintf("day+1 treated coverage: %d/%d firm-events with daret\n",
            sum(!is.na(ev1$daret)), nrow(ev1)))

ar_evt <- ev1[!is.na(group),
              .(Market = mean(ar_sp, na.rm = TRUE),
                CAPM = mean(ar_capm, na.rm = TRUE),
                FF3F = mean(ar_ff3f, na.rm = TRUE)), by = group]

# ---- Diff-in-Means: per-date treated vs all eligible CRSP firms -------------
# (original block filters shrcd/exchcd for BOTH groups, so use the iiwindow
# file alone — shrcd-ineligible treated firms are excluded here as published)
d1 <- unique(tw[event_date == 1, .(permno, date)])
day1 <- inp$ret[date %in% unique(d1$date)]
day1[, treated := paste(permno, date) %chin% paste(d1$permno, d1$date)]
bydate <- day1[, .(m = mean(daret, na.rm = TRUE)), by = .(date, treated)]
bydate <- dcast(bydate, date ~ treated, value.var = "m")
setnames(bydate, c("FALSE", "TRUE"), c("daret0", "daret1"))
bydate[, ar_mean := daret1 - daret0]
bydate[, group := decade_group(year(date))]
ar_dim <- bydate[!is.na(group) & !is.na(ar_mean),
                 .(`Diff-in-Means` = mean(ar_mean)), by = group]

# ---- Gsynth (provisional, from the saved sc_ii_siblis.dta artifact) ---------
sc <- as.data.table(read_dta(ii_work("sc_ii_siblis.dta")))
sc1 <- sc[event_date == 1]
sc1[, ar_sc := daret_treated - daret_sc]
sc1[, group := decade_group(year(anndate))]
ar_sc <- sc1[!is.na(group), .(Gsynth = mean(ar_sc)), by = group]
cat(sprintf("gsynth artifact: %d cohorts / %d firm-events at +1\n",
            uniqueN(sc1$index_anndate), nrow(sc1)))

# ---- assemble (percent) ------------------------------------------------------
tab <- Reduce(function(a, b) merge(a, b, by = "group"),
              list(ar_dim, ar_evt, ar_sc))
tab <- melt(tab, id.vars = "group", variable.name = "col",
            value.name = "estimate")
tab[, `:=`(row = dec_lab[group], estimate = 100 * estimate,
           provisional = fifelse(col == "Gsynth",
                                 "saved sc_ii_siblis.dta artifact", ""))]
col_order <- c("Diff-in-Means", "Market", "CAPM", "FF3F", "Gsynth")
tab <- tab[order(group, match(col, col_order)),
           .(row, col, estimate, provisional)]
write.csv(tab, out_path("table5.csv"), row.names = FALSE)

# ---- compare against published targets --------------------------------------
tg <- as.data.table(read_target(5))
cmp <- merge(tg[, .(row, col, estimate_pub = estimate)],
             tab[, .(row, col, estimate_rep = estimate, provisional)],
             by = c("row", "col"))
cmp[, diff := estimate_rep - estimate_pub]
cmp[, ok := abs(diff) <= 0.1 + 0.005]   # +/-0.1pp plus 2-dp rounding slack
cmp <- cmp[order(row, match(col, col_order))]
print(cmp, digits = 4)
cat(sprintf("\nTable 5: %d/%d cells within tolerance (Gsynth provisional)\n",
            sum(cmp$ok), nrow(cmp)))
write.csv(cmp, out_path("table5_comparison.csv"), row.names = FALSE)
