# NOTE: cfm was subsequently disabled as a public method (issue 31; its
# systematic-effect estimand smears short-lived effects — this script is the
# evidence). Cached outputs remain valid; to refit cfm, re-enable it in
# event_study()/event_study_batch() as described in R/event_study.R.
#
# Table 2 companion — the two latent-factor additions (cfm, apm) on the
# Geithner event, same groups/windows/files as table2.R. Both methods use the
# beforeest files, like the Gsynth column (day 0 is the 3pm-4pm TAQ return in
# Panel A; the all-CRSP beforeest file carries no factor columns, which these
# methods do not need).
#
# cfm: Bai & Wang (arXiv:2606.29691), analytic SEs (auto), Ahn-Horenstein
#      factor count over 1..5 (auto).
# apm: Lei & Ross (arXiv:2312.07520), needs the GitHub apm package;
#      multinomial weighted bootstrap (200 reps, seed 123).
#
# Writes output/table2_new_methods.csv (CFM/APM rows in table2.csv's schema,
# plus the selected r) and prints the combined estimator table.
#
# Run from replication/: Rscript geithner/table2_new_methods.R
suppressMessages({library(feventr); library(haven); library(data.table)})
source("config.R")
stata <- function(f) dind("Geithner/Data and Programs/Stata Files", f)

groups <- c("Schedule connections" = "geithnersched",
            "Personal connections" = "geithner2",
            "New York connections" = "ny")

fit_new <- function(dest_file, pre_start) {
  dest <- as.data.table(read_dta(stata(dest_file)))
  dest <- dest[dif >= pre_start]
  rows <- list()
  for (g in names(groups)) {
    v <- groups[[g]]
    tr <- as.character(unique(dest[get(v) == 1, index_ds]))
    base <- list(data = dest, unit = "index_ds", time = "dif", ret = "ret",
                 treated = tr, event_time = 0, window = c(0, 10),
                 est_window = c(pre_start, -31), returns = "simple",
                 align = "value")
    f_cfm <- do.call(event_study, c(base, list(method = "cfm")))
    f_apm <- do.call(event_study, c(base, list(method = "apm", seed = 123)))
    rows[[g]] <- data.frame(
      row = g, col = c("CFM", "APM"),
      estimate = c(f_cfm$att_avg, f_apm$att_avg),
      se = c(f_cfm$att_avg_se, f_apm$att_avg_se),
      r = c(f_cfm$diagnostics$info$r, f_apm$diagnostics$info$r))
    cat(sprintf("  %s: cfm %.4f (%.4f) | apm %.4f (%.4f)\n", g,
                f_cfm$att_avg, f_cfm$att_avg_se, f_apm$att_avg, f_apm$att_avg_se))
  }
  do.call(rbind, rows)
}

cat("Panel A: bank controls\n")
panA <- fit_new("cleaned_data_beforeest.dta", -256)
cat("Panel B: all CRSP controls\n")
panB <- fit_new("cleaned_data_allcrsp_beforeest.dta", -255)
new <- rbind(cbind(panel = "A: Bank Controls", panA),
             cbind(panel = "B: All Firm Controls", panB))
write.csv(new, out_path("table2_new_methods.csv"), row.names = FALSE)

# combined view against the replicated Table 2 columns
old <- read.csv(out_path("table2.csv"))
keep <- old[old$col %in% c("Average", "DID", "SC", "SDID", "Gsynth", "CAPM"), ]
comb <- rbind(keep[, c("panel", "row", "col", "estimate", "se")],
              new[, c("panel", "row", "col", "estimate", "se")])
comb$col <- factor(comb$col, levels = c("Average", "DID", "SC", "SDID",
                                        "Gsynth", "CAPM", "CFM", "APM"))
wide <- reshape(comb[, c("panel", "row", "col", "estimate")],
                idvar = c("panel", "row"), timevar = "col", direction = "wide")
names(wide) <- sub("estimate[.]", "", names(wide))
cat("\nATT (avg daily return, days 0-10):\n")
print(wide, digits = 3, row.names = FALSE)
