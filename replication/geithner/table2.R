# Table 2 — ATT of Treasury Secretary Announcement (Geithner, Nov 21 2008)
#
# Panel A (bank controls):    cols 1,2,6,7 from cleaned_data_beforesdid.dta
#                             (545 firms x 237 days; pre -256..-31, post 0..10,
#                             placebo days -30..-1 already deleted);
#                             cols 3,4,5 + 8 from cleaned_data_beforeest.dta.
# Panel B (all-CRSP controls): same from the _allcrsp files (pre -255..-31);
#                             cols 3-5 are one-sample t-tests on treated firms
#                             only, hence identical across panels (as published).
# Day 0 is the 3pm-4pm TAQ return (already substituted in the bank/beforeest
# files; the beforesdid_allcrsp file deliberately keeps the full-day CRSP
# day-0 return, matching the published Panel B cols 1,2,6,7 — see PLAN.md).
#
# Conventions per PLAN.md: ATT = average daily return over days 0..10 (no
# cumulation); factor betas estimated in-package on the balanced pre window
# (-256/-255..-31; published used -280..-31 on the unbalanced panel — checked
# against tolerance below); placebo inference 100 reps; gsynth parametric
# bootstrap 1,000 draws, force = "none", r = c(0, 5).
#
# Run from replication/: Rscript geithner/table2.R
suppressMessages({library(feventr); library(haven); library(data.table)})
source("config.R")
stata <- function(f) dind("Geithner/Data and Programs/Stata Files", f)

groups <- c("Schedule connections" = "geithnersched",
            "Personal connections" = "geithner2",
            "New York connections" = "ny")

`%||%` <- function(a, b) if (is.null(a)) b else a

# Panel B's Market/CAPM/FF3F are one-sample t-tests on the treated banks only
# (no control sample), so they are identical across panels; the all-CRSP
# beforeest file carries no factor columns. `factor_rows` reuses Panel A's.
fit_panel <- function(dsdid_file, dest_file, pre_start, factor_rows = NULL) {
  dsdid <- as.data.table(read_dta(stata(dsdid_file)))
  dest <- as.data.table(read_dta(stata(dest_file)))
  dest <- dest[dif >= pre_start]                      # balanced trim as published
  has_factors <- all(c("mktrf", "rf") %in% names(dest))
  if (has_factors) {
    dest[, exret := ret - rf]
    mkt <- unique(dest[, .(dif, market_return)])
    ff <- unique(dest[, .(dif, mktrf, smb, hml)])
  }
  rows <- list()
  for (g in names(groups)) {
    v <- groups[[g]]
    tr_sdid <- as.character(unique(dsdid[get(v) == 1, index_ds]))
    tr_est <- as.character(unique(dest[get(v) == 1, index_ds]))
    base_sdid <- list(data = dsdid, unit = "index_ds", time = "dif", ret = "ret",
                      treated = tr_sdid, event_time = 0, window = c(0, 10),
                      est_window = c(pre_start, -31), returns = "simple",
                      align = "value")
    base_est <- base_sdid
    base_est$data <- dest          # NB: not modifyList — it recurses into data.frames
    base_est$treated <- tr_est
    f <- list(
      Average = do.call(event_study, c(base_sdid, list(method = "mean", se = "tstat"))),
      DID     = do.call(event_study, c(base_sdid, list(method = "did", se = "placebo",
                                                       reps = 100, seed = 123))),
      SC      = do.call(event_study, c(base_sdid, list(method = "sc", se = "placebo",
                                                       reps = 100, seed = 123))),
      SDID    = do.call(event_study, c(base_sdid, list(method = "sdid", se = "placebo",
                                                       reps = 100, seed = 123))),
      Gsynth  = do.call(event_study, c(base_est, list(method = "gsynth", force = "none",
                                                      r = c(0, 5), se = "bootstrap",
                                                      reps = 1000)))
    )
    if (has_factors) {
      f$Market <- do.call(event_study, c(base_est, list(method = "market",
                                                        factors = mkt, se = "tstat")))
      f$CAPM <- do.call(event_study, c(base_est, list(method = "factor",
                                                      factors = mkt, se = "tstat")))
      base_x <- base_est
      base_x$ret <- "exret"
      f$FF3F <- do.call(event_study, c(base_x, list(method = "factor", factors = ff,
                                                    se = "tstat")))
    }
    rows[[g]] <- data.frame(row = g, col = names(f),
                            estimate = sapply(f, `[[`, "att_avg"),
                            se = sapply(f, function(x) x$att_avg_se %||% NA_real_),
                            row.names = NULL)
    cat(sprintf("  %s done (%d treated, %d donors)\n", g, length(tr_sdid),
                f$Average$diagnostics$n_donors))
  }
  tab <- do.call(rbind, rows)
  if (!is.null(factor_rows)) tab <- rbind(tab, factor_rows)
  tab
}

set.seed(2026)   # gsynth parametric bootstrap (unseeded in the original)
cat("Panel A: bank controls\n")
panA <- fit_panel("cleaned_data_beforesdid.dta", "cleaned_data_beforeest.dta", -256)
cat("Panel B: all CRSP controls\n")
panB <- fit_panel("cleaned_data_beforesdid_allcrsp.dta",
                  "cleaned_data_allcrsp_beforeest.dta", -255,
                  factor_rows = panA[panA$col %in% c("Market", "CAPM", "FF3F"), ])

out <- rbind(cbind(panel = "A: Bank Controls", panA),
             cbind(panel = "B: All Firm Controls", panB))
write.csv(out, out_path("table2.csv"), row.names = FALSE)

# ---- compare against the published targets ----------------------------------
tg <- read_target(2)
tg <- tg[tg$row %in% names(groups) & !is.na(tg$estimate), ]
cmp <- merge(tg[, c("panel", "row", "col", "estimate", "se")], out,
             by = c("panel", "row", "col"), suffixes = c("_pub", "_rep"))
cmp$diff <- cmp$estimate_rep - cmp$estimate_pub
cmp$pt_ok <- abs(cmp$diff) <= 0.001
cmp$se_ok <- is.na(cmp$se_pub) | is.na(cmp$se_rep) |
  abs(cmp$se_rep - cmp$se_pub) <= 0.2 * cmp$se_pub + 5e-4  # published rounded to 3dp
print(cmp[order(cmp$panel, cmp$row, cmp$col),
          c("panel", "row", "col", "estimate_pub", "estimate_rep", "diff",
            "pt_ok", "se_pub", "se_rep", "se_ok")], digits = 3)
cat(sprintf("\nPoints within +/-0.001: %d/%d | SEs within 20%%: %d/%d\n",
            sum(cmp$pt_ok), nrow(cmp), sum(cmp$se_ok), nrow(cmp)))
write.csv(cmp, out_path("table2_comparison.csv"), row.names = FALSE)
