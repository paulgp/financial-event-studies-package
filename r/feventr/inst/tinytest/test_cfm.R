# Causal factor model engine (Bai & Wang 2026) vs a hand-rolled reference.
# cfm is implemented but disabled as a public method (issue 31): its
# systematic-effect estimand smears short-lived effects across the event
# window (see replication table1/table5_new_methods.R). The engine stays in
# the codebase and stays verified here at the internal level.

set.seed(42)
sim <- feventr::simulate_events(n_units = 60, n_pre = 50, n_candidate = 1,
                                n_post = 5, tau = 0.02, seed = 99)
ev <- sim$event_time
treated_ids <- sim$events$unit

p <- feventr:::fes_panel(sim$data, "id", "t", "ret", treated = treated_ids,
                         event_time = ev, window = c(0, 5),
                         est_window = c(-50, -1))
eng <- feventr:::eng_cfm(p$Y, p$N0, p$T0, r = 2, se = TRUE)

# --- point estimates: PCA on demeaned donors, pre/post lm on (1, factors) ----
wide <- reshape(sim$data, idvar = "id", timevar = "t", direction = "wide")
times <- (ev - 50):(ev + 5)
Tn <- length(times)
ipre <- 1:50
ipost <- 51:56
donor_ids <- setdiff(unique(sim$data$id), as.integer(treated_ids))
Yd <- as.matrix(wide[match(donor_ids, wide$id), paste0("ret.", times)])
Xc <- Yd - rowMeans(Yd)
egr <- eigen(crossprod(Xc), symmetric = TRUE)
Fh <- egr$vectors[, 1:2, drop = FALSE] * sqrt(Tn)
trt_path <- colMeans(as.matrix(
  wide[match(as.integer(treated_ids), wide$id), paste0("ret.", times),
       drop = FALSE]))
m0 <- lm(trt_path[ipre] ~ Fh[ipre, ])
m1 <- lm(trt_path[ipost] ~ Fh[ipost, ])
tau_ref <- as.vector(cbind(1, Fh[ipost, ]) %*% (coef(m1) - coef(m0)))
expect_equivalent(unname(eng$tau), tau_ref, tolerance = 1e-8)

# contract identity and stored loadings
expect_equivalent(trt_path[ipost] - eng$y0hat[ipost], eng$tau,
                  tolerance = 1e-10)
expect_equal(dim(eng$weights$beta), c(2L, 3L))
expect_equal(rownames(eng$weights$beta), c("pre", "post"))
expect_equal(eng$info$r, 2L)

# --- analytic SEs: HC1 blocks + factor-estimation term (Prop 1 / Lemma 4) ----
Zp <- cbind(1, Fh[ipre, ])
Zq <- cbind(1, Fh[ipost, ])
e0 <- unname(resid(m0))
e1 <- unname(resid(m1))
hc <- function(Z, e) {
  Zi <- solve(crossprod(Z))
  Zi %*% crossprod(Z * e) %*% Zi * nrow(Z) / (nrow(Z) - ncol(Z))
}
C <- hc(Zp, e0) + hc(Zq, e1)
Lam <- Xc %*% Fh / Tn
Ec <- Xc - Lam %*% t(Fh)
dl <- (coef(m1) - coef(m0))[-1]
N0d <- nrow(Yd)
G <- as.vector(Lam %*% solve(crossprod(Lam) / N0d, dl))
nf <- N0d * (N0d - 4L)
se_ref <- sqrt(rowSums((Zq %*% C) * Zq) +
                 colSums((G * Ec[, ipost])^2) / nf)
expect_equivalent(unname(eng$info$se$att), se_ref, tolerance = 1e-8)

zbar <- colMeans(Zq)
ebar <- rowMeans(Ec[, ipost])
avg_ref <- sqrt(sum(zbar * as.vector(C %*% zbar)) + sum((G * ebar)^2) / nf)
expect_equivalent(eng$info$se$avg, avg_ref, tolerance = 1e-8)

# --- eigenvalue-ratio factor-count selection and guards ----------------------
eng_ah <- feventr:::eng_cfm(p$Y, p$N0, p$T0, r = c(0, 5), se = FALSE)
expect_true(eng_ah$info$r %in% 1:5)
expect_error(feventr:::eng_cfm(p$Y, p$N0, p$T0, r = 20),
             pattern = "between 1 and")

# --- disabled as a public option ---------------------------------------------
expect_error(feventr::event_study(sim$data, "id", "t", "ret",
                                  treated = treated_ids, event_time = ev,
                                  method = "cfm", window = c(0, 5),
                                  est_window = c(-50, -1),
                                  returns = "simple"),
             pattern = "arg")
expect_error(feventr::event_study(sim$data, "id", "t", "ret",
                                  treated = treated_ids, event_time = ev,
                                  method = "mean", window = c(0, 5),
                                  est_window = c(-50, -1), returns = "simple",
                                  se = "analytic"),
             pattern = "arg")
