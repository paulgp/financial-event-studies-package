# Build data/ff_daily.rda from the Ken French data library daily factors file.
# Source: F-F_Research_Data_Factors_daily.CSV (Fama/French 3 factors, daily),
# https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html
# The copy read here is the one used by the paper's simulations
# (PEAD_DinD/data/raw), spanning 1926-07-01 .. 2022-11-30, so the bundled
# factors reproduce the published Table 1 DGP exactly.
raw <- read.csv("~/Dropbox/PEAD_DinD/data/raw/F-F_Research_Data_factors_daily.CSV")
names(raw) <- c("date", "mktrf", "smb", "hml", "rf")
ff_daily <- data.frame(
  date  = as.Date(as.character(raw$date), format = "%Y%m%d"),
  mktrf = as.numeric(raw$mktrf) / 100,
  smb   = as.numeric(raw$smb) / 100,
  hml   = as.numeric(raw$hml) / 100,
  rf    = as.numeric(raw$rf) / 100
)
stopifnot(!anyNA(ff_daily), nrow(ff_daily) == 25378)
save(ff_daily, file = "data/ff_daily.rda", compress = "xz")
