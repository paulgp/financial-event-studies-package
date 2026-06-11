#!/usr/bin/env python3
"""Extract per-cohort slices from the index-inclusion event panels.

The panel .dta files (data/work/siblis_anndates/panel_ii_*.dta, up to 27GB,
sorted by permno not cohort) are streamed once with pandas' chunked Stata
reader; rows for requested cohorts are buffered and written as one
cohort_<i>.csv.gz per index_anndate to a cache OUTSIDE the git repo
(licensed CRSP content must never enter the repository).

Usage: extract_cohorts.py <panel.dta> <lo> <hi> <cache_dir>
Cohorts in [lo, hi] already present in cache_dir are skipped; if all are
present the panel scan is skipped entirely.
"""
import os
import sys

import pandas as pd

COLS = ["permno", "event_date", "daret", "include", "index_anndate"]


def main(panel, lo, hi, cache):
    os.makedirs(cache, exist_ok=True)
    want = {i for i in range(lo, hi + 1)
            if not os.path.exists(os.path.join(cache, f"cohort_{i}.csv.gz"))}
    if not want:
        print(f"{os.path.basename(panel)}: cohorts {lo}-{hi} already cached")
        return
    print(f"{os.path.basename(panel)}: scanning for {len(want)} cohorts")
    parts = {}
    nrows = 0
    with pd.read_stata(panel, columns=COLS, chunksize=5_000_000,
                       convert_dates=False, preserve_dtypes=False) as it:
        for chunk in it:
            nrows += len(chunk)
            sub = chunk[chunk["index_anndate"].isin(want)]
            if len(sub):
                for i, g in sub.groupby("index_anndate"):
                    parts.setdefault(int(i), []).append(g)
            if nrows % 50_000_000 < 5_000_000:
                print(f"  ... {nrows:,} rows scanned", flush=True)
    for i, gs in sorted(parts.items()):
        df = pd.concat(gs, ignore_index=True)
        df.sort_values(["permno", "event_date"], inplace=True)
        tmp = os.path.join(cache, f".cohort_{i}.csv.gz.tmp")
        df.to_csv(tmp, index=False, compression="gzip")
        os.replace(tmp, os.path.join(cache, f"cohort_{i}.csv.gz"))
        print(f"  cohort {i}: {len(df):,} rows, "
              f"{df['permno'].nunique():,} units", flush=True)
    missing = sorted(want - set(parts))
    if missing:
        print(f"  WARNING: cohorts not found in panel: {missing}")
    print(f"done ({nrows:,} rows scanned)")


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])
