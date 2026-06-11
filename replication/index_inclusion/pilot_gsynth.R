# PILOT — feventr gsynth vs the saved sc_ii_siblis.dta artifact, cohorts 1-15.
#
# Validates the feventr call (method="gsynth", force="unit", r=c(1,10),
# window=c(-100,20), est_window=c(-280,-101), se="none") against the original
# per-cohort gsynth Y.tr/Y.ct before launching the full 635-cohort run
# (table5_gsynth_full.R). Writes the same checkpoint CSVs the full run uses,
# plus output/pilot_gsynth_agreement.csv (per-cohort correlation of the
# effect path over [-100,20] and the day +1 effect difference).
#
# Run from replication/: Rscript index_inclusion/pilot_gsynth.R

suppressMessages({library(haven); library(data.table); library(parallel)})
source("config.R")
source("index_inclusion/betas_common.R")
source("index_inclusion/gsynth_cohort.R")

outdir <- "index_inclusion/gsynth_out"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
cohorts <- 1:15

st <- mclapply(cohorts, run_cohort, outdir = outdir, mc.cores = 5)
for (i in seq_along(cohorts))
  cat(sprintf("cohort %d: %s\n", cohorts[i], st[[i]]))

# ---- agreement vs the saved artifact ----------------------------------------
sc <- as.data.table(read_dta(ii_work("sc_ii_siblis.dta")))
agree <- rbindlist(lapply(cohorts, function(i) {
  f <- file.path(outdir, sprintf("cohort_%d.csv", i))
  if (!file.exists(f)) return(data.table(index_anndate = i, status = "failed"))
  mine <- fread(f)
  sv <- sc[index_anndate == i,
           .(eff_saved = mean(daret_treated) - mean(daret_sc),
             tr_saved = mean(daret_treated), n_sv = .N), by = event_date]
  if (nrow(sv) == 0) return(data.table(index_anndate = i, status = "no-saved"))
  m <- merge(mine, sv, by = "event_date")
  w <- m[event_date >= -100]
  data.table(
    index_anndate = i, status = "ok",
    n_treat = mine$n_treat[1], r_cv = mine$r_cv[1],
    cor_eff = cor(w$att, w$eff_saved),
    mae_eff = mean(abs(w$att - w$eff_saved)),
    cor_tr_path = cor(m$treated_mean, m$tr_saved),
    max_tr_diff = max(abs(m$treated_mean - m$tr_saved)),
    day1_mine = w[event_date == 1, att],
    day1_saved = w[event_date == 1, eff_saved],
    day1_absdiff = abs(w[event_date == 1, att - eff_saved]))
}), fill = TRUE)
print(agree, digits = 4)
cat(sprintf("\nPilot: %d/%d cohorts fit; median path corr %.4f; MAE(day+1 effect) %.6f\n",
            sum(agree$status == "ok", na.rm = TRUE), length(cohorts),
            median(agree$cor_eff, na.rm = TRUE),
            mean(agree$day1_absdiff, na.rm = TRUE)))
write.csv(agree, out_path("pilot_gsynth_agreement.csv"), row.names = FALSE)
