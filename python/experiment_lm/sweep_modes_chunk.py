#!/usr/bin/env python3
"""Per-chunk LM-mode compare runner: --system square vs --system all-bends.

Reads CLERS strings (one per line) from stdin, runs lm_alpha_march_carry.py
twice per CLERS (square then all-bends, both with --residual quat), and
writes one TSV record per CLERS to stdout:

  clers V E sq_status sq_resid sq_steps sq_retreats sq_wall \\
              ab_status ab_resid ab_steps ab_retreats ab_wall

Status one of: SOLVER_TOL, LAMBDA_SATURATED_POSITIVE_RESIDUAL,
STALLED_POSITIVE_RESIDUAL, REALIZE_FAIL.

No internal parallelism; the caller chunks the corpus and runs many
chunks under GNU parallel. Errors per CLERS are written to stderr and
the row is emitted with REALIZE_FAIL status so downstream summarizers
keep their indexing.

CLI:
  sweep_modes_chunk.py [--target-deg 60] [--init-step-deg 5]
                       [--final-tol 1e-12]
"""
import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]   # python/
THIS = ROOT / "experiment_lm" / "lm_alpha_march_carry.py"

SUMMARY_RE = re.compile(
    r"summary:\s+clers=\S+\s+target_reached=(\S+)\s+steps=(\d+)\s+"
    r"retreats=(\d+)\s+final_resid=(\S+)\s+wall=([\d.]+)s\s+status=(\S+)"
)
HEADER_VE_RE = re.compile(r"V=(\d+)\s+E=(\d+)")


def run_mode(clers: str, system: str, args) -> dict:
    """Run lm_alpha_march_carry.py once and parse the summary."""
    cmd = [sys.executable, str(THIS), clers,
           "--residual", "quat", "--system", system,
           "--target-deg", str(args.target_deg),
           "--init-step-deg", str(args.init_step_deg),
           "--final-tol", str(args.final_tol)]
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    out = {"status": "REALIZE_FAIL", "resid": float("nan"),
           "steps": -1, "retreats": -1, "wall": -1.0,
           "V": -1, "E": -1}
    if p.returncode != 0:
        sys.stderr.write(f"[{system}] {clers}: rc={p.returncode}\n")
        return out
    txt = p.stdout
    mh = HEADER_VE_RE.search(txt)
    if mh:
        out["V"] = int(mh.group(1)); out["E"] = int(mh.group(2))
    ms = SUMMARY_RE.search(txt)
    if ms:
        out["steps"]    = int(ms.group(2))
        out["retreats"] = int(ms.group(3))
        out["resid"]    = float(ms.group(4))
        out["wall"]     = float(ms.group(5))
        out["status"]   = ms.group(6)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-deg",    type=float, default=60.0)
    ap.add_argument("--init-step-deg", type=float, default=5.0)
    ap.add_argument("--final-tol",     type=float, default=1e-12)
    args = ap.parse_args()

    print("clers\tV\tE\t"
          "sq_status\tsq_resid\tsq_steps\tsq_retreats\tsq_wall\t"
          "ab_status\tab_resid\tab_steps\tab_retreats\tab_wall")
    for line in sys.stdin:
        clers = line.strip()
        if not clers or clers.startswith("#"):
            continue
        sq = run_mode(clers, "square",    args)
        ab = run_mode(clers, "all-bends", args)
        V = sq["V"] if sq["V"] > 0 else ab["V"]
        E = sq["E"] if sq["E"] > 0 else ab["E"]
        print(f"{clers}\t{V}\t{E}\t"
              f"{sq['status']}\t{sq['resid']:.6e}\t{sq['steps']}\t"
              f"{sq['retreats']}\t{sq['wall']:.3f}\t"
              f"{ab['status']}\t{ab['resid']:.6e}\t{ab['steps']}\t"
              f"{ab['retreats']}\t{ab['wall']:.3f}", flush=True)


if __name__ == "__main__":
    main()
