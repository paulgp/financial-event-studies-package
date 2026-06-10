# Support-restriction sensitivity on REAL return data (Phase 2 gate c).
# Donor matrix: the Geithner all-CRSP pre-event window (PEAD_DinD,
# cleaned_data_beforesdid_allcrsp.dta: ~4,095 firms x 225 pre days), target =
# the schedule-connected treated mean path. Checks that the hybrid solver's
# objective is insensitive to the support-size cutoff on data with a real
# (not simulated) factor structure.
# Output: benchmarks/support_sensitivity_results.csv
suppressMessages({library(feventr); library(haven); library(data.table)})

dta <- "~/Dropbox/PEAD_DinD/Geithner/Data and Programs/Stata Files/cleaned_data_beforesdid_allcrsp.dta"
d <- as.data.table(read_dta(path.expand(dta)))
pre <- d[dif >= -255 & dif <= -31]
W <- dcast(pre, index_ds ~ dif, value.var = "ret")
ids <- W$index_ds
M <- as.matrix(W[, -1])
keep <- rowSums(is.na(M)) == 0
M <- M[keep, ]; ids <- ids[keep]
treated_ids <- unique(d[geithnersched == 1, index_ds])
tr <- ids %in% treated_ids
A <- t(M[!tr, ])                 # t0 x n0 donor pre paths
b <- colMeans(M[tr, , drop = FALSE])
cat(sprintf("real-data problem: n0=%d donors, t0=%d pre days, %d treated\n",
            ncol(A), nrow(A), sum(tr)))

t_qp <- system.time(s_qp <- solve_simplex_ls(A, b, method = "qp"))[["elapsed"]]
cat(sprintf("full QP: %.2fs, obj %.6e, support %d\n",
            t_qp, s_qp$objective, length(s_qp$support)))

t0 <- nrow(A)
rows <- lapply(c(1, 2, 5, 10) * t0, function(k) {
  t_hy <- system.time(s <- solve_simplex_ls(A, b, method = "hybrid",
                                            support_size = k))[["elapsed"]]
  rel <- (s$objective - s_qp$objective) / s_qp$objective
  cat(sprintf("hybrid support<=%4d (%2dx t0): %.2fs (%4.0fx), rel obj %+.2e, support used %d\n",
              k, k / t0, t_hy, t_qp / t_hy, rel, length(s$support)))
  data.frame(support_cap = k, cap_x_t0 = k / t0, sec = t_hy,
             speedup = t_qp / t_hy, rel_obj = rel,
             support_used = length(s$support))
})
res <- do.call(rbind, rows)
res$qp_sec <- t_qp; res$n0 <- ncol(A); res$t0 <- t0
write.csv(res, "support_sensitivity_results.csv", row.names = FALSE)
cat("\nGATE:", if (all(abs(res$rel_obj[res$cap_x_t0 >= 5]) < 1e-4)) "PASS" else "FAIL", "\n")
