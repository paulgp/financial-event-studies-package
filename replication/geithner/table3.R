# Table 3 — Treated and control betas, Geithner announcement
#
# The published betas are estimated on dif in [-280, -31] (unbalanced
# eligibility panel, >=225 obs). The cleaned panels only start at -256, so
# the -280-window betas are taken from the loadings the data files carry
# (beta/intercept, bmkt_3f/..., and the firm_*_beta_allcrsp companions) —
# these reproduce the published Panel A means exactly, which also pins down
# the window as the only difference from the package's own event_betas()
# estimates on the available balanced window (reported for comparison).
# Panel B applies feventr's SC/SDID omega weights to the same -280 betas
# (equivalent to the original's Stata sdid e(omega) weights once the panels
# are aligned correctly — align = "value", since the placebo days -30..-1
# are pre-deleted from these files; see PLAN.md).
#
# Run from replication/: Rscript geithner/table3.R
suppressMessages({library(feventr); library(haven); library(data.table)})
source("config.R")
stata <- function(f) dind("Geithner/Data and Programs/Stata Files", f)

bank <- as.data.table(read_dta(stata("cleaned_data_beforesdid.dta")))
bank_est <- as.data.table(read_dta(stata("cleaned_data_beforeest.dta")))[dif >= -256]
crsp <- as.data.table(read_dta(stata("cleaned_data_beforesdid_allcrsp.dta")))

# -280-window loadings carried by the data files
b_bank <- unique(bank_est[, .(index_ds, beta, bmkt_3f, bsmb_3f, bhml_3f)])
b_capm_ac <- as.data.table(read_dta(stata("firm_capm_beta_allcrsp.dta")))
b_ff_ac <- as.data.table(read_dta(stata("firm_ff3f_beta_allcrsp.dta")))
setnames(b_capm_ac, intersect(c("beta"), names(b_capm_ac)), "beta", skip_absent = TRUE)
b_ac <- merge(b_capm_ac[, .(index_ds, beta)],
              b_ff_ac[, .(index_ds, bmkt_3f, bsmb_3f, bhml_3f)], by = "index_ds")

tr_bank <- as.character(unique(bank[geithnersched == 1, index_ds]))
tr_crsp <- as.character(unique(crsp[geithnersched == 1, index_ds]))

rows <- list()
add <- function(row, col, est) rows[[length(rows) + 1]] <<- data.frame(
  row = row, col = col, estimate = est)
mlab <- c(beta = "CAPM Beta", bmkt_3f = "FF3F Market Beta",
          bsmb_3f = "FF3F Size Beta", bhml_3f = "FF3F Value Beta")

# ---- Panel A: simple averages (published -280 window) ------------------------
for (v in names(mlab)) {
  add(mlab[v], "Treated", mean(b_bank[index_ds %in% as.numeric(tr_bank)][[v]]))
  add(mlab[v], "Control", mean(b_bank[!index_ds %in% as.numeric(tr_bank)][[v]]))
  # control set = the all-CRSP PANEL's control firms (the original merges the
  # beta files into the panel), not every firm in the beta companion files
  ac_controls <- setdiff(unique(crsp$index_ds), as.numeric(tr_crsp))
  add(mlab[v], "Control (All CRSP)",
      mean(b_ac[index_ds %in% ac_controls][[v]], na.rm = TRUE))
}

# ---- Panel B: omega-weighted control betas (feventr SC/SDID weights) ---------
wfit <- function(d, tr, pre, m)
  event_study(d, "index_ds", "dif", "ret", treated = tr, event_time = 0,
              method = m, window = c(0, 10), est_window = c(pre, -31),
              returns = "simple", se = "none", align = "value")
for (spec in list(list(d = bank, tr = tr_bank, pre = -256, B = b_bank,
                       lab = "Bank Controls"),
                  list(d = crsp, tr = tr_crsp, pre = -255, B = b_ac,
                       lab = "All CRSP Controls"))) {
  for (m in c("sc", "sdid")) {
    om <- wfit(spec$d, spec$tr, spec$pre, m)$weights$omega
    bb <- spec$B[match(as.numeric(names(om)), index_ds)]
    keep <- !is.na(bb$beta)
    for (v in names(mlab))
      add(mlab[v], paste(spec$lab, toupper(m)),
          weighted.mean(bb[[v]][keep], om[keep]))
    cat(spec$lab, m, "done\n")
  }
}

out <- do.call(rbind, rows)
write.csv(out, out_path("table3.csv"), row.names = FALSE)

# package's own event_betas on the available balanced window, for the record
mkt <- unique(bank_est[, .(dif, market_return)])
fb <- event_study(bank, "index_ds", "dif", "ret", treated = tr_bank,
                  event_time = 0, method = "sc", window = c(0, 10),
                  est_window = c(-256, -31), returns = "simple", se = "none",
                  align = "value")
eb <- event_betas(fb, mkt, time = "dif")
cat("\nevent_betas() CAPM on balanced window (-256..-31): treated",
    round(eb[eb$group == "treated" & eb$stat == "mean", "market_return"], 3),
    "control", round(eb[eb$group == "control" & eb$stat == "mean", "market_return"], 3),
    " (published -280 window: 1.427 / 0.825 — window-vintage difference)\n")

# ---- compare against the published targets ----------------------------------
tg <- read_target(3)
tg$col <- sub("^Control \\(banks\\)$", "Control", tg$col)
cmp <- merge(tg[!is.na(tg$estimate), c("row", "col", "estimate")], out,
             by = c("row", "col"), suffixes = c("_pub", "_rep"))
cmp$diff <- cmp$estimate_rep - cmp$estimate_pub
cmp$ok <- abs(cmp$diff) <= 0.0105   # betas published to 3 dp; |diff| <= 0.01 + rounding
print(cmp[order(cmp$row, cmp$col), ], digits = 3, row.names = FALSE)
cat(sprintf("\nWithin +/-0.01: %d/%d (Panel A exact-match expected; Panel B SMB/HML\n",
            sum(cmp$ok), nrow(cmp)),
    "carry the documented Stata-vs-R omega difference)\n")
write.csv(cmp, out_path("table3_comparison.csv"), row.names = FALSE)
