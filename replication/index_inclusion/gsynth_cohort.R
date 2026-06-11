# Shared per-cohort gsynth fit for the index-inclusion Table 5 Gsynth column.
#
# Original: 4_index_inclusion_siblis_gsynth.R — per announcement-date cohort,
# gsynth(daret ~ treated) with treated = include * 1(event_date >= -100),
# force = "unit", r = c(1, 10) (CV), se = FALSE, on the zero-filled panel.
# feventr equivalent: method = "gsynth", window = c(-100, 20),
# est_window = c(-280, -101), returns = "simple", se = "none".
#
# fit_cohort() reads cohort_<i>.csv.gz from the cache (written by
# extract_cohorts.py) and returns a small per-cohort data.table:
# index_anndate, event_date (-280..20), treated_mean, synthetic_mean, att
# (window days only), n_treat, n_donors, r_cv — aggregate estimates only, no
# licensed firm-level data.

suppressMessages({library(feventr); library(data.table)})

fit_cohort <- function(i, cache = feventr_cache("cohorts")) {
  f <- file.path(cache, sprintf("cohort_%d.csv.gz", i))
  if (!file.exists(f)) stop("missing cohort extract: ", f)
  d <- fread(f)
  d[, permno := as.character(permno)]
  tr <- unique(d[include == 1, permno])
  if (length(tr) == 0L) stop("zero treated firms in panel for cohort ", i)
  fit <- event_study(
    data = d, unit = "permno", time = "event_date", ret = "daret",
    treated = tr, event_time = 0, method = "gsynth", force = "unit",
    r = c(1, 10), window = c(-100, 20), est_window = c(-280, -101),
    returns = "simple", se = "none", keep_data = FALSE)
  out <- data.table(index_anndate = i,
                    event_date = as.integer(names(fit$paths$treated)),
                    treated_mean = as.numeric(fit$paths$treated),
                    synthetic_mean = as.numeric(fit$paths$synthetic))
  att <- data.table(event_date = as.integer(names(fit$att)),
                    att = as.numeric(fit$att))
  out <- merge(out, att, by = "event_date", all.x = TRUE)
  out[, `:=`(n_treat = length(tr),
             n_donors = fit$diagnostics$n_donors %||% NA_integer_,
             r_cv = unname(fit$diagnostics$info$r %||% NA_integer_))]
  setcolorder(out, "index_anndate")
  out[]
}

`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a

# Checkpointed runner: writes gsynth_out/cohort_<i>.csv (or _failed marker)
# and skips work already done, so the full run is restartable.
run_cohort <- function(i, outdir, cache = feventr_cache("cohorts")) {
  done <- file.path(outdir, sprintf("cohort_%d.csv", i))
  fail <- file.path(outdir, sprintf("cohort_%d_failed.txt", i))
  if (file.exists(done)) return("done")
  if (file.exists(fail)) return("failed-prev")
  if (!file.exists(file.path(cache, sprintf("cohort_%d.csv.gz", i))))
    return("no-extract")   # transient: extraction stage incomplete, no marker
  res <- tryCatch(fit_cohort(i, cache), error = function(e) e)
  if (inherits(res, "error")) {
    writeLines(conditionMessage(res), fail)
    return(paste("failed:", conditionMessage(res)))
  }
  fwrite(res, done)
  "done"
}
