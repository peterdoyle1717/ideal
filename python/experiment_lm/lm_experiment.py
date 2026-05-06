#!/usr/bin/env python3
"""LM-vs-Newton experiment driver for puffup.py.

Takes a CLERS string, computes ideal bends via horou + dihedral, runs
both solvers, writes per-iter logs, and prints a summary.
"""
import json
import math
import subprocess
import sys
from pathlib import Path

ROOT = Path("/Users/doyle/Dropbox/neo/ideal/python")
sys.path.insert(0, str(ROOT))

import horou
import dihedral as dh
from puffup import parse_netcode, edge_str

CLERS_BIN = "/Users/doyle/Dropbox/neo/clers/bin/clers"
PUFFUP = ROOT / "puffup.py"

def decode_clers(clers: str) -> str:
    p = subprocess.run([CLERS_BIN, "decode"], input=clers + "\n",
                       capture_output=True, text=True, check=True)
    return p.stdout.strip()

def build_payload(clers: str) -> dict:
    netcode = decode_clers(clers)
    faces = parse_netcode(netcode)
    u = horou.horou(faces)
    bends = dh.all_dihedrals(faces, u)
    bends_json = {edge_str(e): float(b) for e, b in bends.items()}
    return {
        "netcode": netcode,
        "base_face": list(faces[0]),
        "corner_angle": math.pi / 3,
        "homotopy": {"from_corner": 0.0},
        "initial_bends": bends_json,
    }

def run(payload: dict, solver: str, log_path: str) -> dict:
    payload = dict(payload)
    payload["solver"] = solver
    payload["iter_log"] = log_path
    p = subprocess.run(
        ["python3", str(PUFFUP), "--in", "-", "--out", "-",
         "--solver", solver, "--iter-log", log_path],
        input=json.dumps(payload), capture_output=True, text=True,
    )
    if p.returncode != 0:
        return {"_error": p.stderr[-1500:], "_returncode": p.returncode}
    return json.loads(p.stdout)

def summarize(clers: str, label: str, solver: str, result: dict, log_path: str) -> str:
    if "_error" in result:
        return f"  [{solver}] FAIL rc={result['_returncode']}: {result['_error'][-200:]}"
    sol = result["solver"]
    log_lines = sum(1 for _ in open(log_path))
    return (f"  [{solver}] success={sol['success']} "
            f"iters={sol['iterations_total']} "
            f"resid={sol['final_residual_norm']:.2e} "
            f"spread={sol['vertex_copy_spread']:.2e} "
            f"steps={len(sol['homotopy_steps'])} "
            f"log_lines={log_lines}")

def main():
    if len(sys.argv) < 3:
        print("usage: lm_experiment.py <label> <CLERS> [<CLERS>...]", file=sys.stderr)
        return 2
    label = sys.argv[1]
    cases = sys.argv[2:]
    for clers in cases:
        try:
            payload = build_payload(clers)
        except Exception as e:
            print(f"== {clers} == BUILD FAIL: {e}")
            continue
        print(f"== {label} :: {clers} == "
              f"V={len({v for f in parse_netcode(payload['netcode']) for v in f})} "
              f"E={len(payload['initial_bends'])}")
        for solver in ("newton", "lm"):
            log = f"/tmp/lm_exp_{label}_{solver}_{clers[:12]}.log"
            r = run(payload, solver, log)
            print(summarize(clers, label, solver, r, log))

if __name__ == "__main__":
    sys.exit(main() or 0)
