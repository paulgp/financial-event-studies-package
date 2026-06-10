# simulate_events(): DGP correctness and bit-exact agreement with the
# published simulation panels where the PEAD_DinD Dropbox is available.

sim <- feventr::simulate_events(seed = 7)
expect_equal(nrow(sim$data), 500 * 500)
expect_equal(sim$event_time, 240)
# treated share is binomial around 10%
expect_true(abs(mean(sim$betas$treated) - 0.1) < 0.06)
# loadings have the right moments
expect_true(abs(mean(sim$betas$b_mkt) - 1) < 0.05)
expect_true(abs(sd(sim$betas$b_smb) - 0.3) < 0.05)

# reproducibility
sim2 <- feventr::simulate_events(seed = 7)
expect_identical(sim$data$ret, sim2$data$ret)

# timing selection picks the max-SMB candidate day
simt <- feventr::simulate_events(selection = "timing", seed = 11)
cand <- 240:489
expect_equal(simt$event_time,
             cand[which.max(simt$factors$smb[cand])])

# assignment selection treats low-SMB-loading firms more often
sima <- feventr::simulate_events(selection = "assignment", seed = 13)
expect_true(mean(sima$betas$b_smb[sima$betas$treated]) <
              mean(sima$betas$b_smb[!sima$betas$treated]))

# --- bit-exact replication of the published simulated panels -----------------
# The published runs saved their panels for the assignment-selection config;
# seed for sim_idx i is 1234 + i - 1.
ref_csv <- path.expand(file.path(
  "~/Dropbox/PEAD_DinD/output/simulations/selection_2factors",
  "simul_data_500_10_240_one-shot-treatment_TRUE_FALSE_1.csv"))
if (file.exists(ref_csv)) {
  ref <- utils::read.csv(ref_csv)
  mine <- feventr::simulate_events(selection = "assignment", seed = 1234)
  d <- merge(mine$data, mine$betas, by = "id")
  d <- merge(d, mine$factors, by = "t")
  key <- order(ref$id, ref$t)
  keym <- order(d$id, d$t)
  expect_equal(d$ret[keym], ref$ret[key], tolerance = 1e-12)
  expect_equal(d$b_smb[keym], ref$b_smb[key], tolerance = 1e-12)
  expect_equal(d$treated[keym], ref$treated[key] == "TRUE" | ref$treated[key] == TRUE)
  expect_equal(unique(ref$treatment_period), mine$event_time)
} else {
  exit_file("PEAD_DinD reference panels not available")
}
