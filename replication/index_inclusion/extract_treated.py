#!/usr/bin/env python3
"""Extract all treated (include==1) firm-day rows from the index-inclusion
event panels into one small csv.gz in the cache dir (outside the git repo).

The original Stata pipeline merged treated firms against raw CRSP daily with
NO shrcd/exchcd filter, so shrcd-ineligible treated firms (REITs etc.) are
absent from permno_date_ret_in_iiwindow.dta. The panels keep treated firms
regardless of share code (only requiring complete 301-day presence); `ret`
retains true missingness while `daret` was zero-filled, so daret is treated
as missing where ret is NA.

Usage: extract_treated.py <cache_dir> <panel1.dta> [<panel2.dta> ...]
"""
import os
import sys

import pandas as pd

COLS = ["permno", "date", "anndate", "event_date", "ret", "daret",
        "include", "index_anndate"]


def main(cache, panels):
    out = os.path.join(cache, "treated_rows.csv.gz")
    if os.path.exists(out):
        print(f"{out} already exists")
        return
    parts = []
    for panel in panels:
        nrows = 0
        with pd.read_stata(panel, columns=COLS, chunksize=5_000_000,
                           convert_dates=True, preserve_dtypes=False) as it:
            for chunk in it:
                nrows += len(chunk)
                sub = chunk[chunk["include"] == 1]
                if len(sub):
                    parts.append(sub.drop(columns=["include"]))
        n = sum(len(p) for p in parts)
        print(f"{os.path.basename(panel)}: {nrows:,} rows scanned, "
              f"{n:,} treated rows so far", flush=True)
    df = pd.concat(parts, ignore_index=True)
    df.sort_values(["index_anndate", "permno", "event_date"], inplace=True)
    tmp = out + ".tmp"
    df.to_csv(tmp, index=False, compression="gzip")
    os.replace(tmp, out)
    print(f"wrote {out}: {len(df):,} rows, "
          f"{df.groupby(['permno', 'anndate']).ngroups} firm-events")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
