#!/usr/bin/env python3
"""Per-chunk full-jump-policy pipeline: LM (all-bends, quat, full jump)
   -> bends -> euclid_realize -> euclid_check.

Reads CLERS strings (one per line) from stdin, runs the full pipeline
once per CLERS, and writes one TSV record per CLERS to stdout. The
caller chunks the corpus and dispatches under GNU parallel.

TSV columns:
  clers V E
  lm_status lm_resid wall_lm
  reached_60 first_jump_to_target retreats_before_first attempts retreats
  accepted_fraction_seq rejected_fraction_seq alpha_seq_deg
  pipe_status dent embed length

pipe_status one of:
  PASS                      LM SOLVER_TOL + realize ok + check PASS
  CHECK_FAIL                LM SOLVER_TOL + realize ok + check failed
                            (i.e. dent<=0 or embed=0 or length too big)
  REALIZE_FAIL              LM SOLVER_TOL but realize/check rejected
  SOLVER_FAIL               LM did not reach SOLVER_TOL

Env:
  REALIZE_BIN   path to euclid_realize          (default: search known)
  CHECK_BIN     path to euclid_check            (default: search known)
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LM   = ROOT / "experiment_lm" / "lm_alpha_march_carry.py"

REALIZE_RC = re.compile(
    r"euclid_check:\s+dent=(\S+)\s+embed=(\d+)\s+defect=\S+\s+"
    r"length=(\S+)\s+sep=\S+\s+->\s+(\S+)"
)
SUMMARY_RE = re.compile(
    r"summary:\s+clers=\S+\s+target_reached=(\S+)\s+steps=\d+\s+"
    r"retreats=(\d+)\s+attempts=(\d+)\s+final_resid=(\S+)\s+wall=([\d.]+)s\s+"
    r"status=(\S+)\s+policy=\S+\s+first_jump_to_target=(\S+)\s+"
    r"retreats_before_first_accept=(\d+)"
)
HEADER_VE_RE = re.compile(r"V=(\d+)\s+E=(\d+)")
ACC_SEQ_RE = re.compile(r"^stats: accepted_fraction_seq=(.+)$", re.M)
REJ_SEQ_RE = re.compile(r"^stats: rejected_fraction_seq=(.+)$", re.M)
ALPHA_SEQ_RE = re.compile(r"^stats: alpha_seq_deg=(.+)$", re.M)


def find_bin(envvar, candidates):
    p = os.environ.get(envvar)
    if p and Path(p).is_file() and os.access(p, os.X_OK):
        return p
    for c in candidates:
        if Path(c).is_file() and os.access(c, os.X_OK):
            return c
    return None


REALIZE_BIN = find_bin("REALIZE_BIN", [
    str(Path.home() / "neo/euclid_realize/bin/euclid_realize"),
    str(Path.home() / "Dropbox/neo/euclid_realize/bin/euclid_realize"),
])
CHECK_BIN = find_bin("CHECK_BIN", [
    str(Path.home() / "neo/euclid_check/bin/euclid_check"),
    str(Path.home() / "Dropbox/neo/euclid_check/bin/euclid_check"),
])


def run_one(clers, args, tmpdir):
    bends = Path(tmpdir) / f"{clers[:30]}.bends"
    obj   = Path(tmpdir) / f"{clers[:30]}.obj"
    out = {
        "V": -1, "E": -1,
        "lm_status": "UNKNOWN", "lm_resid": float("nan"), "wall_lm": -1.0,
        "reached_60": False, "first_jump": False,
        "retreats_before_first": -1, "attempts": -1, "retreats": -1,
        "acc_seq": "-", "rej_seq": "-", "alpha_seq": "-",
        "pipe_status": "SOLVER_FAIL",
        "dent": float("nan"), "embed": -1, "length": float("nan"),
    }
    cmd = [sys.executable, str(LM), clers,
           "--residual", "quat", "--system", "all-bends",
           "--alpha-jump-policy", "full",
           "--target-deg", str(args.target_deg),
           "--quick-min-drop", str(args.quick_min_drop),
           "--max-quick-iter", str(args.max_quick_iter),
           "--final-tol", str(args.final_tol),
           "--bends-out", str(bends)]
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    txt = p.stdout
    mh = HEADER_VE_RE.search(txt)
    if mh:
        out["V"] = int(mh.group(1)); out["E"] = int(mh.group(2))
    ms = SUMMARY_RE.search(txt)
    if ms:
        out["reached_60"]            = (ms.group(1) == "True")
        out["retreats"]              = int(ms.group(2))
        out["attempts"]              = int(ms.group(3))
        out["lm_resid"]              = float(ms.group(4))
        out["wall_lm"]               = float(ms.group(5))
        out["lm_status"]             = ms.group(6)
        out["first_jump"]            = (ms.group(7) == "True")
        out["retreats_before_first"] = int(ms.group(8))
    for rx, key in ((ACC_SEQ_RE, "acc_seq"),
                    (REJ_SEQ_RE, "rej_seq"),
                    (ALPHA_SEQ_RE, "alpha_seq")):
        mm = rx.search(txt)
        if mm:
            out[key] = mm.group(1).strip()

    if out["lm_status"] != "SOLVER_TOL":
        out["pipe_status"] = "SOLVER_FAIL"
        return out
    if not (REALIZE_BIN and CHECK_BIN and bends.is_file()):
        out["pipe_status"] = "REALIZE_FAIL"
        return out
    rcmd = [REALIZE_BIN, "--bends-in", str(bends), "--obj-out", str(obj),
            "--check-bin", CHECK_BIN]
    rp = subprocess.run(rcmd, capture_output=True, text=True, timeout=120)
    rm = REALIZE_RC.search(rp.stdout + "\n" + rp.stderr)
    if rm:
        out["dent"]   = float(rm.group(1))
        out["embed"]  = int(rm.group(2))
        out["length"] = float(rm.group(3))
        verdict = rm.group(4)
        if rp.returncode == 0 and verdict == "PASS":
            out["pipe_status"] = "PASS"
        else:
            out["pipe_status"] = "CHECK_FAIL"
    else:
        out["pipe_status"] = "REALIZE_FAIL"
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-deg",     type=float, default=60.0)
    ap.add_argument("--quick-min-drop", type=float, default=1e-2)
    ap.add_argument("--max-quick-iter", type=int,   default=10)
    ap.add_argument("--final-tol",      type=float, default=1e-12)
    args = ap.parse_args()

    if REALIZE_BIN is None or CHECK_BIN is None:
        sys.stderr.write(f"euclid_realize/check not found; "
                         f"REALIZE_BIN={REALIZE_BIN} CHECK_BIN={CHECK_BIN}\n")
    print("clers\tV\tE\tlm_status\tlm_resid\twall_lm\t"
          "reached_60\tfirst_jump\tretreats_before_first\tattempts\tretreats\t"
          "accepted_fraction_seq\trejected_fraction_seq\talpha_seq_deg\t"
          "pipe_status\tdent\tembed\tlength")
    with tempfile.TemporaryDirectory(prefix="fulljump_") as tmpdir:
        for line in sys.stdin:
            clers = line.strip()
            if not clers or clers.startswith("#"):
                continue
            r = run_one(clers, args, tmpdir)
            print(f"{clers}\t{r['V']}\t{r['E']}\t"
                  f"{r['lm_status']}\t{r['lm_resid']:.6e}\t{r['wall_lm']:.3f}\t"
                  f"{r['reached_60']}\t{r['first_jump']}\t"
                  f"{r['retreats_before_first']}\t{r['attempts']}\t"
                  f"{r['retreats']}\t"
                  f"{r['acc_seq']}\t{r['rej_seq']}\t{r['alpha_seq']}\t"
                  f"{r['pipe_status']}\t{r['dent']:.6e}\t{r['embed']}\t"
                  f"{r['length']:.6e}", flush=True)


if __name__ == "__main__":
    main()
