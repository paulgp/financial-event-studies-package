# M&A refit comparison — assemble the per-deal CARs from ma_refit_full.R,
# validate against the published deal-level gsynth/market CARs, and produce
# Table 6-style cells (cross-deal mean of the [-1,+1] log CAR, in percent,
# by subsample) for each refit method alongside the published columns.
#
# Writes output/ma_refit_deals.csv (deal level, all methods merged) and
# output/table6_refit_cells.csv.
#
# Run from replication/: Rscript ma/ma_refit_compare.R
source("config.R")
suppressMessages({library(haven); library(data.table)})
ma_work <- function(...) dind("M&A", "data", "work", ...)
ma_out  <- function(...) dind("M&A", "output", ...)

read_method <- function(m) {
  fs <- list.files(file.path("ma", "ma_refit_out", m),
                   pattern = "^cohort_[0-9]+[.]csv$", full.names = TRUE)
  if (!length(fs)) return(NULL)
  d <- rbindlist(lapply(fs, fread), fill = TRUE)
  ok <- d[status == "ok"]
  cat(sprintf("%s: %d cohorts, %d deal fits ok, %d skipped/failed\n",
              m, length(fs), uniqueN(ok, by = c("permno", "date_index")),
              nrow(d[status != "ok"])))
  # [-1,+h] log CARs are path differences against event_date == -2
  ok[, .(car_refit = car_log[event_date == 1L] - car_log[event_date == -2L],
         car_refit_250 = car_log[event_date == 250L] -
           car_log[event_date == -2L],
         r = r[1L], pre_rmse = pre_rmse[1L]),
     by = .(permno, date_index)][, permno := as.numeric(permno)][]
}

methods <- Filter(function(m) dir.exists(file.path("ma", "ma_refit_out", m)),
                  c("sc", "apm", "gsynth"))
refits <- lapply(setNames(methods, methods), read_method)
refits <- Filter(Negate(is.null), refits)

di <- as.data.table(read_dta(
  ma_work("sdc_ma_details_sl_m_cleaned_dateindex_2023.dta")))
saved <- as.data.table(read_dta(
  ma_out("sl_m_deals_car_1_250_gsynth.dta"),
  col_select = c("group", "car_log_sc_1", "car_log_vwretd_1")))
det <- as.data.table(read_dta(
  ma_work("sdc_ma_details_sl_m_cleaned_2023.dta"),
  col_select = c("master_deal_no", "permno", "tpublic", "pct_cash",
                 "pct_stk")))

deals <- merge(di[, .(master_deal_no, permno, group, date_index)],
               saved, by = "group")
deals <- merge(deals, det, by = c("master_deal_no", "permno"))
for (m in names(refits))
  deals <- merge(deals,
                 setnames(copy(refits[[m]]),
                          c("car_refit", "r", "pre_rmse"),
                          paste0(c("car_", "r_", "pre_rmse_"), m)),
                 by = c("permno", "date_index"), all.x = TRUE)
fwrite(deals, out_path("ma_refit_deals.csv"))

cat("\nDeal-level agreement with the published gsynth CARs (car_log_sc_1):\n")
for (m in names(refits)) {
  v <- deals[[paste0("car_", m)]]
  ok <- !is.na(v)
  cat(sprintf("  %-6s N=%5d  cor=%.3f (spearman %.3f)  mean diff=%+.4f\n",
              m, sum(ok), cor(v[ok], deals$car_log_sc_1[ok]),
              cor(v[ok], deals$car_log_sc_1[ok], method = "spearman"),
              mean(v[ok] - deals$car_log_sc_1[ok])))
}

## Table 6-style cells: cross-deal means of the [-1,+1] log CAR, percent
cells <- function(d, cols) {
  sub <- list("Full sample"     = rep(TRUE, nrow(d)),
              "Public targets"  = d$tpublic == "Public",
              "Private targets" = d$tpublic == "Priv.",
              "Other targets"   = d$tpublic == "Sub.",
              "Cash merger"     = d$pct_cash == 100,
              "Stock merger"    = d$pct_stk == 100)
  rbindlist(lapply(names(sub), function(s) {
    i <- which(sub[[s]])
    rbindlist(lapply(names(cols), function(cn) {
      v <- d[[cols[[cn]]]][i]
      data.table(row = cn, col = s, estimate_pct = 100 * mean(v, na.rm = TRUE),
                 n = sum(!is.na(v[])))
    }))
  }))
}
cols <- c("Market (published)" = "car_log_vwretd_1",
          "Gsynth (published)" = "car_log_sc_1")
for (m in names(refits)) cols[paste0(toupper(m), " (refit)")] <- paste0("car_", m)
tab <- cells(deals, as.list(cols))
write.csv(tab, out_path("table6_refit_cells.csv"), row.names = FALSE)
wide <- dcast(tab, row ~ col, value.var = "estimate_pct")
setcolorder(wide, c("row", "Full sample", "Public targets", "Private targets",
                    "Other targets", "Cash merger", "Stock merger"))
cat("\nTable 6-style cells (mean [-1,+1] log CAR, %):\n")
print(wide, digits = 3, row.names = FALSE)
