#!/usr/bin/env python3
"""LM continuation-step sweep on a single hard case.

Tests whether LM can survive larger α-homotopy steps than the cautious
default `step_factor=1/24`. Newton is run once at default as a baseline.

Hard case selected for v=50 first prime — Newton stalls on it under the
default policy (this session, /tmp/lm_exp_large_newton_*.log).
"""
import json
import math
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import horou
import dihedral as dh
from puffup import parse_netcode, edge_str

CLERS_BIN = "/Users/doyle/Dropbox/neo/clers/bin/clers"
PUFFUP = ROOT / "puffup.py"

# Default: v=50 first prime (Newton stalls here under default continuation).
DEFAULT_CLERS = ("CCCACACACACACACACACACACACACACACACACACA"
                 "CACACACACACACACACACACACACACACACACACACACACACACAAE")

def decode(clers: str) -> str:
    return subprocess.run(
        [CLERS_BIN, "decode"], input=clers + "\n",
        capture_output=True, text=True, check=True).stdout.strip()

def build_payload(clers: str, step_factor: float, solver: str,
                  log_path: str, min_step_factor: float = 1e-6) -> dict:
    netcode = decode(clers)
    faces = parse_netcode(netcode)
    u = horou.horou(faces)
    bends = dh.all_dihedrals(faces, u)
    return {
        "netcode": netcode,
        "base_face": list(faces[0]),
        "corner_angle": math.pi / 3,
        "homotopy": {
            "from_corner": 0.0,
            "step_factor": step_factor,
            "min_step_factor": min_step_factor,
        },
        "initial_bends": {edge_str(e): float(b) for e, b in bends.items()},
        "solver": solver,
        "iter_log": log_path,
    }

def run(payload: dict) -> dict:
    t0 = time.monotonic()
    p = subprocess.run(
        ["python3", str(PUFFUP), "--in", "-", "--out", "-"],
        input=json.dumps(payload), capture_output=True, text=True,
    )
    dt = time.monotonic() - t0
    if p.returncode != 0:
        try:
            d = json.loads(p.stdout) if p.stdout else None
        except Exception:
            d = None
        msg = (d or {}).get("solver", {}).get("message", "(no parsed msg)") if d else "(no stdout)"
        return {"_rc": p.returncode, "_msg": msg, "_dt": dt}
    out = json.loads(p.stdout)
    return {"out": out, "_rc": 0, "_dt": dt}

def summarize(label: str, r: dict, log_path: str) -> str:
    if r["_rc"] != 0:
        try:
            log_lines = sum(1 for _ in open(log_path))
        except Exception:
            log_lines = 0
        return (f"  [{label}] FAIL  msg={r.get('_msg','?')}  "
                f"wall={r['_dt']:.2f}s  iter_log_lines={log_lines}")
    sol = r["out"]["solver"]
    log_lines = sum(1 for _ in open(log_path))
    return (f"  [{label}] success={sol['success']}  "
            f"iters={sol['iterations_total']}  "
            f"resid={sol['final_residual_norm']:.3e}  "
            f"steps={len(sol['homotopy_steps'])}  "
            f"wall={r['_dt']:.2f}s  iter_log_lines={log_lines}")

def main():
    clers = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_CLERS
    print(f"CLERS: {clers}  V={len({v for f in parse_netcode(decode(clers)) for v in f})}")
    print()

    # Newton baseline at default step.
    log = "/tmp/lm_step_sweep_newton_default.log"
    r = run(build_payload(clers, 1.0/24.0, "newton", log))
    print(summarize("newton sf=1/24 (baseline)", r, log))
    print()

    # LM at increasing step_factor.
    for sf_label, sf in [("1/24 (default)", 1.0/24.0),
                          ("1/12", 1.0/12.0),
                          ("1/6", 1.0/6.0),
                          ("1/4", 1.0/4.0),
                          ("1/2", 1.0/2.0),
                          ("1 (single shot)", 1.0)]:
        log = f"/tmp/lm_step_sweep_lm_sf_{sf:.4f}.log"
        r = run(build_payload(clers, sf, "lm", log))
        print(summarize(f"lm     sf={sf_label}", r, log))

if __name__ == "__main__":
    main()
