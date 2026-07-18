# Causal factor model (Bai & Wang 2026) vs a hand-rolled PCA + lm reference

set.seed(42)
sim <- feventr::simulate_events(n_units = 60, n_pre = 50, n_candidate = 1,
                                n_post = 5, tau = 0.02, seed = 99)
ev <- sim$event_time
treated_ids <- sim$events$unit

fit <- feventr::event_study(sim$data, "id", "t", "ret",
                            treated = treated_ids, event_time = ev,
                            method = "cfm", window = c(0, 5),
                            est_window = c(-50, -1), returns = "simple", r = 2)

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
expect_equivalent(unname(fit$att), tau_ref, tolerance = 1e-8)

# implied counterfactual keeps the realized idiosyncratic shock, so CARs
# cumulate tau* exactly under the compound convention
expect_equivalent(unname(fit$car),
                  feventr:::car_from_paths(trt_path[ipost],
                                           trt_path[ipost] - tau_ref,
                                           "compound"),
                  tolerance = 1e-8)

# pre/post loadings stored for diagnostics
expect_equal(dim(fit$weights$beta), c(2L, 3L))
expect_equal(rownames(fit$weights$beta), c("pre", "post"))
expect_equal(fit$diagnostics$info$r, 2L)

# --- analytic SEs: HC1 blocks + factor-estimation term (Prop 1 / Lemma 4) ----
expect_equal(fit$se$method, "analytic")
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
expect_equivalent(unname(fit$se$att), se_ref, tolerance = 1e-8)

zbar <- colMeans(Zq)
ebar <- rowMeans(Ec[, ipost])
avg_ref <- sqrt(sum(zbar * as.vector(C %*% zbar)) + sum((G * ebar)^2) / nf)
expect_equivalent(fit$att_avg_se, avg_ref, tolerance = 1e-8)

# --- factor-count selection and refit freezing -------------------------------
fit_ah <- feventr::event_study(sim$data, "id", "t", "ret",
                               treated = treated_ids, event_time = ev,
                               method = "cfm", window = c(0, 5),
                               est_window = c(-50, -1), returns = "simple")
expect_true(fit_ah$diagnostics$info$r %in% 1:5)

# placebo inference reruns the engine on donors (with the selected r frozen)
fit_pl <- feventr::event_study(sim$data, "id", "t", "ret",
                               treated = treated_ids, event_time = ev,
                               method = "cfm", window = c(0, 5),
                               est_window = c(-50, -1), returns = "simple",
                               se = "placebo", reps = 19, seed = 1)
expect_equal(fit_pl$se$method, "placebo")
expect_equal(length(fit_pl$se$att), 6L)
expect_equivalent(unname(fit_pl$att), unname(fit_ah$att), tolerance = 1e-10)

# --- guards ------------------------------------------------------------------
expect_error(feventr::event_study(sim$data, "id", "t", "ret",
                                  treated = treated_ids, event_time = ev,
                                  method = "cfm", window = c(0, 5),
                                  est_window = c(-50, -1), returns = "simple",
                                  se = "conformal"),
             pattern = "conformal")
expect_error(feventr::event_study(sim$data, "id", "t", "ret",
                                  treated = treated_ids, event_time = ev,
                                  method = "cfm", window = c(0, 5),
                                  est_window = c(-50, -1), returns = "simple",
                                  se = "tstat"),
             pattern = "tstat")
expect_error(feventr::event_study(sim$data, "id", "t", "ret",
                                  treated = treated_ids, event_time = ev,
                                  method = "mean", window = c(0, 5),
                                  est_window = c(-50, -1), returns = "simple",
                                  se = "analytic"),
             pattern = "analytic")
expect_error(feventr::event_study(sim$data, "id", "t", "ret",
                                  treated = treated_ids, event_time = ev,
                                  method = "cfm", window = c(0, 5),
                                  est_window = c(-50, -1), returns = "simple",
                                  r = 20),
             pattern = "between 1 and")

# --- batch mode --------------------------------------------------------------
b <- feventr::event_study_batch(sim$data, "id", "t", "ret",
                                events = data.frame(unit = treated_ids,
                                                    event_time = ev),
                                method = "cfm", window = c(0, 5),
                                est_window = c(-50, -1), returns = "simple",
                                r = 2)
expect_equivalent(unname(b$att), unname(fit$att), tolerance = 1e-10)
