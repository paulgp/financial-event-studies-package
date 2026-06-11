# Shared helpers for the index-inclusion replication (Tables 4-5).
#
# Per-firm CAPM / FF3F betas over event days [-250, -101] (the original
# clean_index_dates_siblis.do keeps `event_date >= -250 & event_date < -100`
# and listwise-drops missing daret before the per-event OLS of
# exret = daret - rf on mktrf [+ smb + hml]).
#
# Returns sources:
#  * data/work/permno_date_ret_in_iiwindow.dta — every CRSP firm (shrcd
#    10/11/12, exchcd 1-3) on every date inside any event window, with daret
#    missingness preserved. Used for random controls (always eligible) and as
#    the primary treated source.
#  * treated_rows.csv.gz (extracted from the cohort panels by
#    extract_treated.py into the cache OUTSIDE the repo) — the panels keep
#    treated firms regardless of share code, so this covers the 38 treated
#    firm-events (REITs etc.) absent from the iiwindow file; `ret` retains
#    true missingness while `daret` was zero-filled, so daret is set NA where
#    ret is NA. The union mirrors the original unfiltered merge against
#    crsp_daily_raw_2023_cleaned; coverage is checked against the saved beta
#    artifacts in table4.R.
#
# Sourced from scripts run with cwd = replication/ (after source("config.R")).

suppressMessages({
  library(haven)
  library(data.table)
})

ii_work <- function(f) dind("index_inclusion/data/work", f)
feventr_cache <- function(...) {
  file.path(Sys.getenv("FEVENTR_CACHE",
                       path.expand("~/.cache/feventr_replication")), ...)
}

load_ii_inputs <- function() {
  ev <- as.data.table(read_dta(ii_work("include_event_date_siblis.dta")))
  ret <- as.data.table(read_dta(ii_work("permno_date_ret_in_iiwindow.dta")))
  ff <- as.data.table(read_dta(ii_work("FF5F_daily.dta")))
  ff <- ff[, .(date, mktrf, smb, hml, rf)]
  setkey(ret, permno, date)
  list(ev = ev, ret = ret, ff = ff)
}

# Per-group OLS of exret on the given factor columns; returns one row per id
# with alpha + betas (NA when the regression cannot be run, mirroring Stata's
# `capture reg ... continue`).
fit_betas <- function(dt, id_col, xvars) {
  prefix <- if (length(xvars) == 1L) "capm" else "ff3f"
  est <- dt[!is.na(daret), {
    X <- cbind(1, as.matrix(.SD[, xvars, with = FALSE]))
    y <- daret - rf
    if (.N >= ncol(X) && qr(X)$rank == ncol(X)) {
      b <- qr.solve(X, y)
      as.list(b)
    } else {
      as.list(rep(NA_real_, ncol(X)))
    }
  }, by = id_col]
  setnames(est, setdiff(names(est), id_col), c("alpha", paste0("b_", xvars)))
  est
}

# Treated firm-day returns over the full loaded window: iiwindow primary,
# panel-extracted rows fill shrcd-ineligible treated firms.
load_treated_returns <- function(inp) {
  win <- merge(inp$ev, inp$ret, by = c("permno", "date"), all.x = TRUE)
  pf <- feventr_cache("treated_rows.csv.gz")
  if (file.exists(pf)) {
    pr <- fread(pf)
    pr <- pr[, .(permno, anndate = as.Date(anndate), event_date,
                 daret_panel = fifelse(is.na(ret), NA_real_, daret))]
    win <- merge(win, pr, by = c("permno", "anndate", "event_date"),
                 all.x = TRUE)
    win[, daret := fcoalesce(daret, daret_panel)][, daret_panel := NULL]
  } else {
    warning("treated_rows.csv.gz not in cache; run extract_treated.py — ",
            "shrcd-ineligible treated firms will be missing")
  }
  win
}

# Treated betas: event map (792 firm-events) x returns, [-250, -101].
compute_treated_betas <- function(inp, tw = load_treated_returns(inp)) {
  win <- tw[event_date >= -250 & event_date < -100]
  win <- merge(win, inp$ff, by = "date")
  capm <- fit_betas(win, "index", "mktrf")
  ff3f <- fit_betas(win, "index", c("mktrf", "smb", "hml"))
  setnames(capm, c("alpha", "b_mktrf"), c("alpha_capm", "beta_capm"))
  setnames(ff3f, c("alpha", "b_mktrf", "b_smb", "b_hml"),
           c("alpha_ff3f", "bmkt_ff3f", "bsmb_ff3f", "bhml_ff3f"))
  hdr <- unique(inp$ev[, .(index, permno, anndate, effdate)])
  Reduce(function(a, b) merge(a, b, by = "index", all.x = TRUE),
         list(hdr, capm, ff3f))
}

# Random-control betas: MUST consume the saved unseeded draw
# permno_anndate_random_controls.dta (PLAN.md). Controls share each
# announcement date's event-window calendar.
compute_control_betas <- function(inp) {
  rc <- as.data.table(read_dta(ii_work("permno_anndate_random_controls.dta")))
  cal <- unique(inp$ev[, .(date, anndate, event_date)])
  win <- merge(cal, rc, by = "anndate", allow.cartesian = TRUE)
  win <- win[event_date >= -250 & event_date < -100]
  win <- merge(win, inp$ret, by = c("permno", "date"), all.x = TRUE)
  win <- merge(win, inp$ff, by = "date")
  win[, idx := .GRP, by = .(permno, anndate)]
  capm <- fit_betas(win, "idx", "mktrf")
  ff3f <- fit_betas(win, "idx", c("mktrf", "smb", "hml"))
  setnames(capm, c("alpha", "b_mktrf"), c("alpha_capm", "beta_capm"))
  setnames(ff3f, c("alpha", "b_mktrf", "b_smb", "b_hml"),
           c("alpha_ff3f", "bmkt_ff3f", "bsmb_ff3f", "bhml_ff3f"))
  hdr <- unique(win[, .(idx, permno, anndate)])
  Reduce(function(a, b) merge(a, b, by = "idx", all.x = TRUE),
         list(hdr, capm, ff3f))
}

# Decade groups by year(anndate): 1 = 80-89, ..., 4 = 10-20 (NA outside).
decade_group <- function(y) {
  fcase(y >= 1980 & y <= 1989, 1L, y >= 1990 & y <= 1999, 2L,
        y >= 2000 & y <= 2009, 3L, y >= 2010 & y <= 2020, 4L,
        default = NA_integer_)
}
decade_labels <- c("Panel A: 1980-1989", "Panel B: 1990-1999",
                   "Panel B: 2000-2009",  # paper prints a second "Panel B"
                   "Panel D: 2010-2020")
