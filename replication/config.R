# Shared configuration for the replication scripts.
#
# Scripts are run with the working directory set to replication/ (e.g.
# `cd replication && Rscript geithner/table2.R`). They read cleaned datasets
# from the paper's working repo — licensed CRSP/Compustat/SDC content that is
# NEVER copied into this repository or its git history. Set FEVENTR_PEAD_DIND
# to point elsewhere if needed.

pead_dind <- Sys.getenv("FEVENTR_PEAD_DIND", path.expand("~/Dropbox/PEAD_DinD"))
if (!dir.exists(pead_dind))
  stop("PEAD_DinD data directory not found at ", pead_dind,
       " — set FEVENTR_PEAD_DIND")

dind <- function(...) file.path(pead_dind, ...)
read_target <- function(n) utils::read.csv(sprintf("targets/table%d.csv", n))
out_path <- function(f) { dir.create("output", showWarnings = FALSE); file.path("output", f) }
