# CVXR reference solver --------------------------------------------------
#
# Exact (interior-point) solution of the simplex least-squares problem, used
# only to validate solve_simplex_ls() in tests and benchmarks (CVXR is in
# Suggests). Pattern follows synthdid/R/reference-solver.R.

ref_simplex_ls <- function(A, b, V = NULL, eta = 0) {
  if (!requireNamespace("CVXR", quietly = TRUE))
    stop("ref_simplex_ls requires the CVXR package")
  if (!is.null(V)) {
    if (is.matrix(V)) {
      R <- chol(V)
      A <- R %*% A
      b <- as.vector(R %*% b)
    } else {
      A <- A * sqrt(V)
      b <- b * sqrt(V)
    }
  }
  w <- CVXR::Variable(ncol(A))
  objective <- CVXR::Minimize(
    CVXR::sum_squares(A %*% w - b) + eta * CVXR::sum_squares(w))
  prob <- CVXR::Problem(objective, list(sum(w) == 1, w >= 0))
  cvxr_solve <- get("solve", envir = asNamespace("CVXR"))
  sol <- cvxr_solve(prob)
  as.vector(suppressWarnings(sol$getValue(w)))
}
