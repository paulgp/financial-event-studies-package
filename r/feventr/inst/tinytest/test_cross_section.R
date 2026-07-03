# cross_section(): OLS vs lm, GLS (cholesky vs direct vs woodbury), the
# pre-event benchmarking statistics, guards, and S3 methods.

set.seed(7)
N <- 120; T <- 360
x <- rnorm(N)
z <- rnorm(N)
fac <- rnorm(T, 0, 0.02)
# returns load on a common factor through x, plus idiosyncratic noise, plus a
# true event-day relationship with x; balanced panel (no missing)
ret <- outer(fac, 0.5 * x) + matrix(rnorm(N * T, 0, 0.01), T, N)
ev <- 340L
ret[ev, ] <- ret[ev, ] + 0.03 * x
long <- data.frame(id = rep(seq_len(N), each = T), t = rep(seq_len(T), times = N),
                   ret = as.vector(ret),                 # firm-major
                   x = rep(x, each = T), z = rep(z, each = T))

# --- OLS: event and pre-event slopes equal lm() on those cross-sections ------
fo <- feventr::cross_section(long, "id", "t", "ret", event_time = ev,
                             chars = c("x", "z"), pre_window = c(-100, -1),
                             method = "ols")
expect_inherits(fo, "fes_cse")
expect_equal(names(fo$params), c("x", "z", "(Intercept)"))
expect_equal(dim(fo$pre_betas), c(100L, 3L))
expect_equal(fo$n_obs[1], N)                              # balanced: all firms

d0 <- long[long$t == ev, ]
expect_equal(unname(fo$params[c("x", "z")]),
             unname(coef(lm(ret ~ x + z, d0))[c("x", "z")]), tolerance = 1e-10)
# a pre-event date: offset -20 -> t = ev-20 -> pre_betas row 20 (latest first)
d20 <- long[long$t == ev - 20L, ]
expect_equal(unname(fo$pre_betas[20L, c("x", "z")]),
             unname(coef(lm(ret ~ x + z, d20))[c("x", "z")]), tolerance = 1e-10)

# --- intercept-only model = mean return per date -----------------------------
fm <- feventr::cross_section(long, "id", "t", "ret", event_time = ev,
                             chars = NULL, pre_window = c(-100, -1))
expect_equal(unname(fm$params[["(Intercept)"]]), mean(d0$ret), tolerance = 1e-12)

# --- GLS: cholesky == direct Omega^-1 formula == woodbury --------------------
fg <- feventr::cross_section(long, "id", "t", "ret", event_time = ev,
                             chars = "x", pre_window = c(-100, -1),
                             method = "gls", npc = 15L)
fw <- feventr::cross_section(long, "id", "t", "ret", event_time = ev,
                             chars = "x", pre_window = c(-100, -1),
                             method = "gls", npc = 15L, solver = "woodbury")
# direct GLS on the event date, independent of the whitening/woodbury paths
win <- (ev - 100L):(ev - 1L)
dec <- feventr:::.cse_pca(ret[win, ], 15L)
Om <- dec$V %*% diag(dec$lam) %*% t(dec$V) + diag(dec$d)
Xd <- cbind(x, 1); Oi <- solve(Om)
b_direct <- solve(t(Xd) %*% Oi %*% Xd, t(Xd) %*% Oi %*% ret[ev, ])
expect_equal(unname(fg$params[["x"]]), b_direct[1], tolerance = 1e-8)
expect_equal(unname(fg$params), unname(fw$params), tolerance = 1e-8)

# rolling GLS window: a pre-event date's beta uses its OWN preceding window
tgt <- ev - 30L
win_t <- (tgt - 100L):(tgt - 1L)
dec_t <- feventr:::.cse_pca(ret[win_t, ], 15L)
Om_t <- dec_t$V %*% diag(dec_t$lam) %*% t(dec_t$V) + diag(dec_t$d)
Oi_t <- solve(Om_t)
b_t <- solve(t(Xd) %*% Oi_t %*% Xd, t(Xd) %*% Oi_t %*% ret[tgt, ])
expect_equal(unname(fg$pre_betas[30L, "x"]), b_t[1], tolerance = 1e-8)

# --- significance statistics recomputed by hand ------------------------------
all_b <- rbind(fo$params, fo$pre_betas)
pm <- colMeans(fo$pre_betas)
psd <- sqrt(apply(fo$pre_betas, 2L, var))
Lp1 <- nrow(all_b)
p_cdf_hand <- colSums(abs(sweep(all_b, 2, pm)) >=
                        matrix(abs(fo$params - pm), Lp1, ncol(all_b), byrow = TRUE)) / Lp1
expect_equal(unname(fo$p_cdf), unname(p_cdf_hand), tolerance = 1e-12)
zc <- abs(fo$params - pm) / (psd * sqrt(Lp1 / (Lp1 - 1)))
p_par_hand <- 2 * pt(zc, df = Lp1 - 2, lower.tail = FALSE)
expect_equal(unname(fo$p_parametric), unname(p_par_hand), tolerance = 1e-12)
# p-values live on [0, 1]; the empirical CDF p is a multiple of 1/(L+1)
expect_true(all(fo$p_cdf >= 0 & fo$p_cdf <= 1))
expect_true(all(abs(fo$p_cdf * Lp1 - round(fo$p_cdf * Lp1)) < 1e-9))

# the event-day slope recovers the true positive relationship with x
expect_true(fo$params[["x"]] > 0)

# --- guards ------------------------------------------------------------------
expect_error(feventr::cross_section(long, "id", "t", "ret", event_time = ev,
                                    chars = "x", pre_window = c(-1, -100)),
             "pre_window")
expect_error(feventr::cross_section(long, "id", "t", "ret", event_time = ev,
                                    chars = "x", pre_window = c(-5, 3)),
             "pre_window")
expect_error(feventr::cross_section(long, "id", "t", "ret", event_time = ev,
                                    chars = "x", pre_window = c(-10, -1),
                                    method = "gls", npc = 50L),
             "npc")
# GLS needs history another -pre_window[1] before the first pseudo-event day:
# event at t=150 with pre_window -100 needs data back to t=-50 (absent)
expect_error(feventr::cross_section(long, "id", "t", "ret", event_time = 150L,
                                    chars = "x", pre_window = c(-100, -1),
                                    method = "gls", npc = 15L),
             "gls")
expect_error(feventr::cross_section(long, "id", "t", "ret", event_time = 99999L,
                                    chars = "x", pre_window = c(-100, -1)),
             "not found")
expect_error(feventr::cross_section(long, "id", "t", "ret", event_time = ev,
                                    chars = "nope", pre_window = c(-100, -1)),
             "not found")
expect_error(feventr::cross_section(long, "id", "t", "ret", event_time = ev,
                                    chars = "x", pre_window = c(-100, -1),
                                    method = "ols", solver = "woodbury"),
             "woodbury")

# --- S3 methods --------------------------------------------------------------
expect_stdout(print(fo), pattern = "cross-sectional")
s <- summary(fo)
expect_equal(nrow(s), 3L)
expect_true(all(c("estimate", "p_cdf", "p_parametric") %in% names(s)))
expect_equal(coef(fo), fo$params)
tf <- tempfile(fileext = ".png")
grDevices::png(tf)
expect_silent(plot(fo, which = "x"))
expect_silent(plot(fm))          # intercept-only
grDevices::dev.off()
unlink(tf)
