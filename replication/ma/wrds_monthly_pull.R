# Pull CRSP monthly returns from WRDS for the monthly (3-year-horizon)
# M&A analysis: crsp.msf joined to msenames (US common shares, NYSE/AMEX/
# NASDAQ, shrcd 10/11, exchcd 1-3) with delisting returns from
# crsp.msedelist merged on the delisting month. ret is total return
# (includes dividends); ret_adj compounds the delisting return in the
# delisting month (Shumway -30% imputation for performance delists with
# missing dlret; moot for complete-coverage designs but kept for
# completeness). Credentials: ~/.pgpass (wrds-pgdata).
#
# Writes <PEAD_DinD>/M&A/data/work/crsp_monthly_wrds_1974_2024.csv.gz
# (licensed CRSP content -- lives in the paper's Dropbox workspace, never
# in this repository).
#
# Run from replication/: Rscript ma/wrds_monthly_pull.R
source("config.R")
suppressMessages({library(RPostgres); library(data.table)})

out_file <- dind("M&A", "data", "work", "crsp_monthly_wrds_1974_2024.csv.gz")
pg <- readLines(path.expand("~/.pgpass"))
user <- strsplit(pg[grepl("wrds", pg)][1], ":")[[1]][4]
con <- dbConnect(Postgres(), host = "wrds-pgdata.wharton.upenn.edu",
                 port = 9737, dbname = "wrds", user = user,
                 sslmode = "require")
cat("connected as", user, "\n")
q <- "
select m.permno, m.date, m.ret, d.dlret, d.dlstcd
from crsp.msf m
join crsp.msenames n
  on m.permno = n.permno and m.date between n.namedt and n.nameendt
left join crsp.msedelist d
  on m.permno = d.permno
 and date_trunc('month', m.date) = date_trunc('month', d.dlstdt)
where m.date >= '1974-01-01'
  and n.shrcd in (10, 11)
  and n.exchcd in (1, 2, 3)
"
mo <- as.data.table(dbGetQuery(con, q))
dbDisconnect(con)
cat("rows:", nrow(mo), " permnos:", uniqueN(mo$permno),
    " span:", as.character(range(mo$date)), "\n")

# Shumway: performance delists (500, 520-584) with missing dlret -> -30%
mo[, dlret_use := dlret]
mo[is.na(dlret_use) & !is.na(dlstcd) &
     (dlstcd == 500L | (dlstcd >= 520L & dlstcd <= 584L)),
   dlret_use := -0.30]
mo[, ret_adj := fifelse(is.na(dlret_use), ret,
                        (1 + fcoalesce(ret, 0)) * (1 + dlret_use) - 1)]
mo[, month := format(as.Date(date), "%Y-%m")]
fwrite(mo[, .(permno, month, ret = ret_adj)], out_file)
cat("wrote", out_file, ":", nrow(mo), "rows\n")
