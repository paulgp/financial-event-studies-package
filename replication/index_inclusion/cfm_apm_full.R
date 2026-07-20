# NOTE: cfm was subsequently disabled as a public method (issue 31; its
# systematic-effect estimand smears short-lived effects — this script is the
# evidence). Cached outputs remain valid; to refit cfm, re-enable it in
# event_study()/event_study_batch() as described in R/event_study.R.
#
# cfm and apm over all 635 index-inclusion cohorts, mirroring the Gsynth
# full run's conventions (window=c(-100,20), est_window=c(-280,-101),
# returns="simple", se="none"; factor count by Ahn-Horenstein over 1..5) and
# its checkpointed layout: one CSV per cohort in index_inclusion/cfm_out/ and
# apm_out/, already-done cohorts skipped, failures recorded as
# cohort_<i>_failed.txt markers.
#
# Reads the per-cohort cache built by table5_gsynth_full.R stage 1 (run that
# first if the cache is empty). Going through the per-cohort files also keeps
# each worker's memory at one cohort panel rather than the 66.7M-row master
# panel (see issue #21). apm needs the GitHub apm package (Remotes in
# DESCRIPTION); each engine call pins its own internal threads, and workers
# pin data.table to one thread.
#
# Run from replication/:
#   nohup Rscript index_inclusion/cfm_apm_full.R > output/cfm_apm_full.log 2>&1 &
suppressMessages({library(feventr); library(data.table); library(parallel)})
source("config.R")
source("index_inclusion/betas_common.R")

cache <- feventr_cache("cohorts")
if (!length(list.files(cache, pattern = "^cohort_.*[.]csv[.]gz$")))
  stop("cohort cache is empty at ", cache,
       " — run table5_gsynth_full.R stage 1 first")
`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a

fit_one <- function(i, method, outdir) {
  done <- file.path(outdir, sprintf("cohort_%d.csv", i))
  fail <- file.path(outdir, sprintf("cohort_%d_failed.txt", i))
  if (file.exists(done) || file.exists(fail)) return("skip")
  f <- file.path(cache, sprintf("cohort_%d.csv.gz", i))
  if (!file.exists(f)) return("no-extract")
  res <- tryCatch({
    d <- fread(f)
    d[, permno := as.character(permno)]
    tr <- unique(d[include == 1, permno])
    fit <- event_study(data = d, unit = "permno", time = "event_date",
                       ret = "daret", treated = tr, event_time = 0,
                       method = method, window = c(-100, 20),
                       est_window = c(-280, -101), returns = "simple",
                       se = "none", keep_data = FALSE)
    out <- data.table(index_anndate = i,
                      event_date = as.integer(names(fit$paths$treated)),
                      treated_mean = as.numeric(fit$paths$treated),
                      synthetic_mean = as.numeric(fit$paths$synthetic))
    att <- data.table(event_date = as.integer(names(fit$att)),
                      att = as.numeric(fit$att))
    out <- merge(out, att, by = "event_date", all.x = TRUE)
    out[, `:=`(n_treat = length(tr),
               n_donors = fit$diagnostics$n_donors %||% NA_integer_,
               r = unname(fit$diagnostics$info$r %||% NA_integer_))]
    setcolorder(out, "index_anndate")
    out
  }, error = function(e) e)
  if (inherits(res, "error")) {
    writeLines(conditionMessage(res), fail)
    return("failed")
  }
  fwrite(res, done)
  "done"
}

for (method in c("cfm", "apm")) {
  outdir <- file.path("index_inclusion", paste0(method, "_out"))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  t0 <- proc.time()[3]
  st <- mclapply(1:635, function(i) {
    data.table::setDTthreads(1L)
    fit_one(i, method, outdir)
  }, mc.cores = 6, mc.preschedule = TRUE)
  nulls <- sum(vapply(st, is.null, TRUE))
  tab <- table(unlist(st))
  cat(sprintf("%s: %s%s | %.1f min\n", method,
              paste(names(tab), tab, collapse = ", "),
              if (nulls) paste0(", killed-fork ", nulls) else "",
              (proc.time()[3] - t0) / 60))
}
