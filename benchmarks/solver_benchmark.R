# Benchmark: feventr::solve_simplex_ls hybrid vs pure FW vs full OSQP
# Gate (Phase 2): >= 10x speedup at n0 >= 2000 with objective within 1e-4
# relative of the full-QP optimum.
# Output: benchmarks/solver_benchmark_results.csv
suppressMessages(library(feventr))
set.seed(42)

mk_problem <- function(n0, t0, k = 3) {
  F_ <- matrix(rnorm(t0 * k), t0, k)
  L <- matrix(abs(rnorm(n0 * k)), n0, k)
  A <- F_ %*% t(L) + matrix(rnorm(t0 * n0, sd = 0.5), t0, n0)
  w_true <- rep(0, n0); w_true[1:5] <- 0.2
  list(A = A, b = as.vector(A %*% w_true + rnorm(t0, sd = 0.1)))
}

grid <- expand.grid(n0 = c(500, 2000, 5000, 10000), t0 = c(100, 250))
rows <- list()
for (i in seq_len(nrow(grid))) {
  n0 <- grid$n0[i]; t0 <- grid$t0[i]
  p <- mk_problem(n0, t0)
  t_qp <- system.time(s_qp <- solve_simplex_ls(p$A, p$b, method = "qp"))[["elapsed"]]
  t_hy <- system.time(s_hy <- solve_simplex_ls(p$A, p$b, method = "hybrid"))[["elapsed"]]
  t_fw <- system.time(s_fw <- solve_simplex_ls(p$A, p$b, method = "fw"))[["elapsed"]]
  rel_hy <- (s_hy$objective - s_qp$objective) / s_qp$objective
  rel_fw <- (s_fw$objective - s_qp$objective) / s_qp$objective
  rows[[i]] <- data.frame(n0 = n0, t0 = t0,
                          qp_sec = t_qp, hybrid_sec = t_hy, fw_sec = t_fw,
                          speedup_hybrid = t_qp / t_hy,
                          rel_obj_hybrid = rel_hy, rel_obj_fw = rel_fw,
                          support_hybrid = length(s_hy$support),
                          fw_iters = s_hy$iterations)
  cat(sprintf("n0=%5d t0=%3d | qp %7.2fs | hybrid %6.2fs (%4.0fx, rel obj %+.1e, support %d) | fw %6.2fs (rel %+.1e)\n",
              n0, t0, t_qp, t_hy, t_qp / t_hy, rel_hy,
              length(s_hy$support), t_fw, rel_fw))
}
res <- do.call(rbind, rows)
write.csv(res, "solver_benchmark_results.csv", row.names = FALSE)
ok <- with(res, all(speedup_hybrid[n0 >= 2000] >= 10) && all(abs(rel_obj_hybrid) < 1e-4))
cat("\nGATE:", if (ok) "PASS" else "FAIL", "\n")
