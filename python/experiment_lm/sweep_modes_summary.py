#!/usr/bin/env python3
"""Aggregate per-chunk sweep_modes_chunk.py TSVs into a per-v summary.

Usage:
  sweep_modes_summary.py SWEEP_DIR

Reads SWEEP_DIR/results/c_*.tsv (each carrying its own header line) and
emits to stdout:
  v   N   sq_TOL  ab_TOL  sq_med_resid  ab_med_resid  sq_max_resid  ab_max_resid
  ...
plus a 'where they disagree' table — CLERS where sq_status != ab_status.

Cheap to run; no external deps beyond stdlib.
"""
import statistics
import sys
from pathlib import Path


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: sweep_modes_summary.py SWEEP_DIR")
    sweep = Path(sys.argv[1])
    results = sweep / "results"
    if not results.is_dir():
        sys.exit(f"no results dir at {results}")

    # rows: list of dicts with V/E/clers and per-mode status,resid
    rows = []
    for tsv in sorted(results.glob("c_*.tsv")):
        with open(tsv) as fh:
            header = fh.readline().strip().split("\t")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) != len(header):
                    continue
                row = dict(zip(header, f))
                try:
                    row["V"] = int(row["V"])
                    row["E"] = int(row["E"])
                    row["sq_resid"] = float(row["sq_resid"])
                    row["ab_resid"] = float(row["ab_resid"])
                except (KeyError, ValueError):
                    continue
                rows.append(row)

    if not rows:
        sys.exit("no rows parsed")

    # group by V
    by_v = {}
    for r in rows:
        by_v.setdefault(r["V"], []).append(r)

    print(f"# total rows: {len(rows)}    distinct V: {len(by_v)}")
    print()
    print(f"{'V':>4}  {'N':>5}  {'sq_TOL':>6}  {'ab_TOL':>6}  "
          f"{'sq_med_r':>10}  {'ab_med_r':>10}  "
          f"{'sq_max_r':>10}  {'ab_max_r':>10}  "
          f"{'sq_med_wall':>11}  {'ab_med_wall':>11}")
    for V in sorted(by_v):
        rs = by_v[V]
        N = len(rs)
        sq_tol = sum(1 for r in rs if r["sq_status"] == "SOLVER_TOL")
        ab_tol = sum(1 for r in rs if r["ab_status"] == "SOLVER_TOL")
        sq_resids = [r["sq_resid"] for r in rs
                     if r["sq_resid"] == r["sq_resid"]]   # not NaN
        ab_resids = [r["ab_resid"] for r in rs
                     if r["ab_resid"] == r["ab_resid"]]
        sq_walls = [float(r["sq_wall"]) for r in rs]
        ab_walls = [float(r["ab_wall"]) for r in rs]
        print(f"{V:>4}  {N:>5}  {sq_tol:>6}  {ab_tol:>6}  "
              f"{statistics.median(sq_resids):>10.3e}  "
              f"{statistics.median(ab_resids):>10.3e}  "
              f"{max(sq_resids):>10.3e}  {max(ab_resids):>10.3e}  "
              f"{statistics.median(sq_walls):>11.3f}  "
              f"{statistics.median(ab_walls):>11.3f}")

    # disagreements
    diffs = [r for r in rows if r["sq_status"] != r["ab_status"]]
    print()
    print(f"# disagreements (sq_status ≠ ab_status): {len(diffs)}")
    if diffs:
        print(f"{'V':>4}  {'sq_status':>30}  {'ab_status':>30}  "
              f"{'sq_resid':>10}  {'ab_resid':>10}  clers")
        for r in diffs[:80]:
            print(f"{r['V']:>4}  {r['sq_status']:>30}  {r['ab_status']:>30}  "
                  f"{r['sq_resid']:>10.3e}  {r['ab_resid']:>10.3e}  "
                  f"{r['clers']}")
        if len(diffs) > 80:
            print(f"... +{len(diffs)-80} more")


if __name__ == "__main__":
    main()
