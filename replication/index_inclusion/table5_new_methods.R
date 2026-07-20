# Table 5 companion — CFM and APM columns for the index-inclusion decade
# table, aggregated from the cohort fits in cfm_out/ and apm_out/ (built by
# cfm_apm_full.R) exactly as table5_gsynth_full.R stage 3 aggregates the
# Gsynth column: the day +1 effect (announcements after close), weighted by
# each cohort's treated count within decade. Reported over the
# published-vintage cohorts (those in sc_ii_siblis.dta, as the paper's
# Gsynth column) with the all-successful-cohorts variant alongside.
#
# Reading the CFM column: under the gsynth-style single treated regime
# (-100..+20), cfm's path varies over t only through the factors, so a
# single-day estimand like Table 5's day +1 effect is ~0 by construction —
# the announcement spike is exactly the non-systematic component cfm
# excludes. Its meaningful index-inclusion quantity is the window-average
# systematic shift (see cfm_apm_compare.R). APM imputes day-specific
# counterfactual means like Gsynth and is directly comparable here.
#
# Writes output/table5_new_methods.csv and prints the combined decade table.
#
# Run from replication/: Rscript index_inclusion/table5_new_methods.R
suppressMessages({library(haven); library(data.table)})
source("config.R")
source("index_inclusion/betas_common.R")

read_day1 <- function(dir) {
  fs <- list.files(file.path("index_inclusion", dir),
                   pattern = "^cohort_[0-9]+[.]csv$", full.names = TRUE)
  if (!length(fs)) stop("no cohort fits in index_inclusion/", dir,
                        " — run cfm_apm_full.R first")
  rbindlist(lapply(fs, fread), fill = TRUE)[event_date == 1L]
}

ev <- as.data.table(read_dta(ii_work("include_event_date_siblis.dta")))
ann <- unique(ev[, .(anndate)])[order(anndate)][, index_anndate := .I]
sc <- as.data.table(read_dta(ii_work("sc_ii_siblis.dta")))
pub_cohorts <- unique(sc$index_anndate)

agg <- function(d) d[!is.na(group),
                     .(estimate = 100 * weighted.mean(att, n_treat),
                       n_events = sum(n_treat), n_cohorts = .N), by = group]
dec_lab <- c("1980-1989", "1990-1999", "2000-2009", "2010-2020")

rows <- list()
for (m in c("CFM", "APM")) {
  d <- read_day1(paste0(tolower(m), "_out"))
  d <- merge(d, ann, by = "index_anndate")
  d[, group := decade_group(year(anndate))]
  a_pub <- agg(d[index_anndate %in% pub_cohorts])
  a_all <- agg(d)
  rows[[m]] <- data.frame(row = dec_lab[a_pub$group], col = m,
                          estimate = a_pub$estimate,
                          estimate_allcohorts = a_all$estimate[
                            match(a_pub$group, a_all$group)],
                          n_events = a_pub$n_events,
                          n_cohorts = a_pub$n_cohorts)
}
new <- do.call(rbind, rows)
write.csv(new[order(new$row, new$col), ], out_path("table5_new_methods.csv"),
          row.names = FALSE)

# combined decade table against the replicated Table 5 columns
old <- read.csv(out_path("table5.csv"))
comb <- rbind(old[, c("row", "col", "estimate")],
              new[, c("row", "col", "estimate")])
comb$col <- factor(comb$col, levels = c("Diff-in-Means", "Market", "CAPM",
                                        "FF3F", "Gsynth", "CFM", "APM"))
wide <- reshape(comb, idvar = "row", timevar = "col", direction = "wide")
names(wide) <- sub("estimate[.]", "", names(wide))
cat("Table 5 (day +1 effect x100, n_treat-weighted within decade):\n")
print(wide[order(wide$row), ], digits = 3, row.names = FALSE)
