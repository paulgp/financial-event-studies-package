# Table 4 — Beta distribution of S&P 500 index-inclusion firms by decade
# (Goldsmith-Pinkham & Lyu, p. 43). 792 firm-events / 635 announcement dates.
#
# Original: index_inclusion/code/gp_lyu_indexinc_replication/
# 1_clean_index_dates_siblis.do -> output/average_beta_index_inclusion
# {,_randcon}_siblis.csv. Per-firm OLS of (daret - rf) on mktrf (CAPM) and on
# mktrf+smb+hml (FF3F) over event days [-250, -101], listwise-dropping missing
# daret; mean/sd by decade of year(anndate) for treated firms and for the
# randomly drawn controls (saved unseeded draw consumed verbatim per PLAN.md).
#
# Run from replication/: Rscript index_inclusion/table4.R

suppressMessages({library(haven); library(data.table)})
source("config.R")
source("index_inclusion/betas_common.R")

inp <- load_ii_inputs()
cat(sprintf("events: %d firm-events / %d anndates; returns rows: %s\n",
            uniqueN(inp$ev$index), uniqueN(inp$ev$anndate),
            format(nrow(inp$ret), big.mark = ",")))

tr <- compute_treated_betas(inp)
rc <- compute_control_betas(inp)

# ---- cross-check against the saved beta artifacts ---------------------------
chk <- function(mine, saved_file, keys, vars) {
  sv <- as.data.table(read_dta(ii_work(saved_file)))
  m <- merge(mine, sv, by = keys, suffixes = c("_my", "_sv"))
  for (v in vars) {
    d <- abs(m[[paste0(v, "_my")]] - m[[paste0(v, "_sv")]])
    cat(sprintf("  %-28s %-10s matched %d/%d  max|diff| %.2e  both-NA %d\n",
                saved_file, v, sum(!is.na(d)), nrow(sv), max(d, na.rm = TRUE),
                sum(is.na(m[[paste0(v, "_my")]]) & is.na(m[[paste0(v, "_sv")]]))))
  }
  invisible(m)
}
cat("Cross-check vs saved Stata betas:\n")
chk(tr, "index_inclusion_siblis_capm_beta.dta", "index",
    c("alpha_capm", "beta_capm"))
chk(tr, "index_inclusion_siblis_ff3f_beta.dta", "index",
    c("alpha_ff3f", "bmkt_ff3f", "bsmb_ff3f", "bhml_ff3f"))
chk(rc, "index_inclusion_siblis_random_control_capm_beta.dta",
    c("permno", "anndate"), c("alpha_capm", "beta_capm"))
chk(rc, "index_inclusion_siblis_random_control_ff3f_beta.dta",
    c("permno", "anndate"), c("alpha_ff3f", "bmkt_ff3f", "bsmb_ff3f", "bhml_ff3f"))

# Persist treated betas for table5.R (cache, outside the git repo).
saveRDS(tr, feventr_cache("ii_treated_betas.rds"))

# ---- decade mean/sd (Stata tabstat, listwise within each factor model) ------
tr[, group := decade_group(year(anndate))]
rc[, group := decade_group(year(anndate))]

summarise <- function(dt, who) {
  capm <- dt[!is.na(group) & !is.na(alpha_capm) & !is.na(beta_capm),
             .(row = "CAPM Beta", mean = mean(beta_capm), sd = sd(beta_capm)),
             by = group]
  f <- dt[!is.na(group) & !is.na(alpha_ff3f) & !is.na(bmkt_ff3f) &
            !is.na(bsmb_ff3f) & !is.na(bhml_ff3f)]
  ff <- rbind(
    f[, .(row = "FF3F Mkt Beta", mean = mean(bmkt_ff3f), sd = sd(bmkt_ff3f)), by = group],
    f[, .(row = "FF3F SMB Beta", mean = mean(bsmb_ff3f), sd = sd(bsmb_ff3f)), by = group],
    f[, .(row = "FF3F HML Beta", mean = mean(bhml_ff3f), sd = sd(bhml_ff3f)), by = group])
  out <- rbind(capm, ff)
  out[, who := who]
  out
}
summ <- rbind(summarise(tr, "Treated"), summarise(rc, "Random Control"))
tab <- melt(summ, id.vars = c("group", "row", "who"),
            measure.vars = c("mean", "sd"), variable.name = "stat")
tab[, panel := decade_labels[group]]
tab[, col := paste(who, fifelse(stat == "mean", "Mean", "Std"))]
tab <- tab[, .(panel, row, col, estimate = value)]

row_order <- c("CAPM Beta", "FF3F Mkt Beta", "FF3F SMB Beta", "FF3F HML Beta")
col_order <- c("Treated Mean", "Treated Std", "Random Control Mean",
               "Random Control Std")
tab <- tab[order(match(panel, decade_labels), match(row, row_order),
                 match(col, col_order))]
write.csv(tab, out_path("table4.csv"), row.names = FALSE)

# ---- compare against published targets --------------------------------------
tg <- as.data.table(read_target(4))
cmp <- merge(tg[, .(panel, row, col, estimate_pub = estimate)],
             tab[, .(panel, row, col, estimate_rep = estimate)],
             by = c("panel", "row", "col"))
cmp[, diff := estimate_rep - estimate_pub]
cmp[, ok := abs(diff) <= 0.001 + 5e-4]  # +/-0.001 plus 3-dp rounding slack
cmp <- cmp[order(match(panel, decade_labels), match(row, row_order),
                 match(col, col_order))]
print(cmp, digits = 4, nrows = 100)
cat(sprintf("\nTable 4: %d/%d cells within tolerance\n", sum(cmp$ok), nrow(cmp)))
write.csv(cmp, out_path("table4_comparison.csv"), row.names = FALSE)
