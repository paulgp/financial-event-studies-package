# Engine correctness against independent references (lm, fixest::feols)

set.seed(42)
sim <- feventr::simulate_events(n_units = 60, n_pre = 50, n_candidate = 1,
                                n_post = 5, tau = 0.02, seed = 99)
ev <- sim$event_time
treated_ids <- sim$events$unit

# --- factor model vs hand-rolled lm() per treated unit -----------------------
fit <- feventr::event_study(sim$data, "id", "t", "ret",
                            treated = treated_ids, event_time = ev,
                            method = "factor", window = c(0, 5),
                            est_window = c(-50, -1), returns = "simple",
                            factors = sim$factors[, c("t", "mktrf", "smb")])

wide <- reshape(sim$data, idvar = "id", timevar = "t", direction = "wide")
Fm <- as.matrix(sim$factors[, c("mktrf", "smb")])
pre_t <- (ev - 50):(ev - 1)
post_t <- ev:(ev + 5)
ar_ref <- sapply(post_t, function(s) {
  mean(sapply(treated_ids, function(u) {
    y <- unlist(wide[wide$id == u, paste0("ret.", c(pre_t, s))])
    m <- lm(y[seq_along(pre_t)] ~ Fm[pre_t, ])
    y[length(y)] - sum(c(1, Fm[s, ]) * coef(m))
  }))
})
expect_equivalent(unname(fit$att), ar_ref, tolerance = 1e-10)

# per-unit loadings stored for diagnostics
expect_equal(dim(fit$weights$beta), c(length(treated_ids), 3L))
u1 <- treated_ids[1]
m1 <- lm(unlist(wide[wide$id == u1, paste0("ret.", pre_t)]) ~ Fm[pre_t, ])
expect_equivalent(unname(fit$weights$beta[u1, ]), unname(coef(m1)),
                  tolerance = 1e-10)

# --- market-adjusted: beta pinned at 1, no intercept --------------------------
mkt <- data.frame(t = sim$factors$t, mkt = sim$factors$mktrf)
fitm <- feventr::event_study(sim$data, "id", "t", "ret",
                             treated = treated_ids, event_time = ev,
                             method = "market", window = c(0, 5),
                             est_window = c(-50, -1), returns = "simple",
                             factors = mkt)
trt_mean <- sapply(post_t, function(s)
  mean(sim$data$ret[sim$data$t == s & sim$data$id %in% as.integer(treated_ids)]))
expect_equivalent(unname(fitm$att), trt_mean - sim$factors$mktrf[post_t],
                  tolerance = 1e-10)

# --- diff-in-means and DiD vs fixest::feols ----------------------------------
if (requireNamespace("fixest", quietly = TRUE)) {
  dd <- sim$data
  dd$treated <- dd$id %in% as.integer(treated_ids)
  dd <- dd[dd$t >= ev - 50 & dd$t <= ev + 5, ]
  dd$post <- dd$t >= ev

  # per-period diff-in-means = i(t, treated) coefficients on the post panel
  fm <- fixest::feols(ret ~ -1 + i(t) + i(t, treated), data = dd[dd$post, ])
  fit_mean <- feventr::event_study(sim$data, "id", "t", "ret",
                                   treated = treated_ids, event_time = ev,
                                   method = "mean", window = c(0, 5),
                                   est_window = c(-50, -1), returns = "simple")
  expect_equivalent(unname(fit_mean$att),
                    unname(coef(fm)[grepl("treated", names(coef(fm)))]),
                    tolerance = 1e-10)

  # DiD average over window = twfe treated x post coefficient
  dd$tp <- dd$treated & dd$post
  fm2 <- fixest::feols(ret ~ tp | id + t, data = dd)
  fit_did <- feventr::event_study(sim$data, "id", "t", "ret",
                                  treated = treated_ids, event_time = ev,
                                  method = "did", window = c(0, 5),
                                  est_window = c(-50, -1), returns = "simple")
  expect_equivalent(fit_did$att_avg, unname(coef(fm2)), tolerance = 1e-10)
}

# --- CAR conventions on a constructed example --------------------------------
trt <- c(0.10, -0.05, 0.02)
syn <- c(0.04, 0.01, -0.01)
expect_equal(feventr:::car_from_paths(trt, syn, "sum"), cumsum(trt - syn))
expect_equal(feventr:::car_from_paths(trt, syn, "compound"),
             cumprod(1 + trt) - cumprod(1 + syn))
expect_equal(feventr:::car_from_paths(trt, syn, "log"),
             cumsum(log(1 + trt) - log(1 + syn)))

# 'auto' picks compound for simple returns, sum for log returns
expect_equal(fit$conventions$cumulate, "compound")
fitl <- feventr::event_study(sim$data, "id", "t", "ret",
                             treated = treated_ids, event_time = ev,
                             method = "mean", window = c(0, 5),
                             est_window = c(-50, -1), returns = "log")
expect_equal(fitl$conventions$cumulate, "sum")
expect_equivalent(fitl$car, cumsum(fitl$att))
