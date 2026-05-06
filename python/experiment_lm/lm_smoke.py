#!/usr/bin/env python3
"""Smoke test: run puffup.py's newton vs lm on a tetrahedron and report.

Tet is symmetric, so this is *only* a sanity test of the LM machinery,
not a real experiment per Codex's no-tet/octa/icosa caveat. Real
experiment will use asymmetric cases below.
"""
import json
import math
import subprocess
import sys

PUFFUP = "/Users/doyle/Dropbox/neo/ideal/python/puffup.py"

# Tetrahedron netcode: 4 faces, each edge bend = 2π/3 in ideal limit.
TET_PAYLOAD = {
    "netcode": "1,2,3;1,3,4;1,4,2;2,4,3",
    "base_face": [1, 2, 3],
    "corner_angle": math.pi / 3,
    "homotopy": {"from_corner": 0.0},
    "initial_bends": 2.0 * math.pi / 3.0,
}

def run(solver):
    payload = dict(TET_PAYLOAD)
    payload["solver"] = solver
    payload["iter_log"] = f"/tmp/lm_smoke_tet_{solver}.log"
    p = subprocess.run(
        ["python3", PUFFUP, "--in", "-", "--out", "-", "--solver", solver,
         "--iter-log", payload["iter_log"]],
        input=json.dumps(payload), capture_output=True, text=True,
    )
    if p.returncode != 0:
        print(f"[{solver}] EXIT {p.returncode}")
        print("stderr:", p.stderr[-500:])
        return None
    return json.loads(p.stdout)

for s in ("newton", "lm"):
    r = run(s)
    if r is None:
        sys.exit(2)
    sol = r["solver"]
    print(f"[{s}] success={sol['success']} iters={sol['iterations_total']} "
          f"final_residual={sol['final_residual_norm']:.3e} "
          f"vertex_copy_spread={sol['vertex_copy_spread']:.3e} "
          f"steps={len(sol['homotopy_steps'])}")
