#!/usr/bin/env python3
"""LM at α=2.5° from 65° start vs ideal start, on small primes (v=4..12).

For each prime CLERS at given v range:
  - Run LM with start = uniform 65°
  - Run LM with start = horou ideal
  - Report success/iters/residual for each, and whether they reached the
    same point in bend space (max|x_65 - x_ideal| < threshold).
"""
import math
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np

import horou
import dihedral as dh
from puffup import (
    parse_netcode, edge_str, edge_key, Tri,
    holonomy_residual, vertex_turn,
    solver_lm,
)

CLERS_BIN = "/Users/doyle/Dropbox/neo/clers/bin/clers"
PRIMES_DIR = Path("/Users/doyle/Dropbox/neo/data/primes")

def decode(clers):
    return subprocess.run(
        [CLERS_BIN, "decode"], input=clers + "\n",
        capture_output=True, text=True, check=True).stdout.strip()

def setup(clers):
    nc = decode(clers)
    faces = parse_netcode(nc)
    tri = Tri.from_faces(faces)
    u = horou.horou(faces)
    bends_ideal = dh.all_dihedrals(faces, u)
    base_face = faces[0]
    base_edges = [edge_key(base_face[i], base_face[(i+1) % 3]) for i in range(3)]
    base_set = set(base_edges)
    var_edges = [e for e in tri.edges if e not in base_set]
    return tri, base_face, base_edges, var_edges, bends_ideal

def run_lm(tri, base_face, base_edges, var_edges, bends_ideal, alpha,
           start_bends_full, lambda_init=1e-3, tol=1e-12, max_iter=200):
    base_bend = {e: bends_ideal[e] for e in base_edges}
    def F(x):
        return holonomy_residual(tri, base_face, var_edges, x, base_bend, alpha)
    base_face_set = set(base_face)
    def gate(x_trial):
        trial = dict(start_bends_full)
        for e, val in zip(var_edges, x_trial):
            trial[e] = float(val)
        for v in tri.vertices:
            if v in base_face_set:
                continue
            if vertex_turn(tri, v, trial) < 0.0:
                return False
        return True
    x0 = np.array([start_bends_full[e] for e in var_edges])
    log = []
    out = solver_lm(F, x0, tol=tol, max_iter=max_iter,
                    lambda_init=lambda_init, iter_log=log, trial_gate=gate)
    return out

def main():
    v_min = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    v_max = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    alpha_deg = 2.5
    alpha = math.radians(alpha_deg)
    lambda_init = 1e-3
    same_thresh = 1e-6  # bend-vector max-abs difference to count as same point

    print(f"alpha = {alpha_deg}°  start = uniform 65°  vs  start = horou ideal")
    print(f"lambda_init = {lambda_init:.0e}  same threshold max|·| < {same_thresh:.0e}")
    print()
    print(f"  {'v':>3}  {'#':>4}  {'CLERS':>20}  "
          f"{'65°':>22}  {'ideal':>22}  {'same?':>6}  {'max|Δ|':>9}")

    n_total = 0
    n_both_success = 0
    n_same = 0

    for v in range(v_min, v_max + 1):
        path = PRIMES_DIR / f"{v}.txt"
        if not path.exists():
            continue
        clersts = [c.strip() for c in path.read_text().splitlines() if c.strip()]
        for i, clers in enumerate(clersts):
            try:
                tri, base_face, base_edges, var_edges, bends_ideal = setup(clers)
            except Exception as e:
                print(f"  {v:>3}  {i+1:>4}  {clers[:20]:>20}  setup fail: {e}")
                continue
            start_65 = {e: math.radians(65.0) for e in tri.edges}
            o65 = run_lm(tri, base_face, base_edges, var_edges,
                           bends_ideal, alpha, start_65, lambda_init=lambda_init)
            oid = run_lm(tri, base_face, base_edges, var_edges,
                           bends_ideal, alpha, bends_ideal, lambda_init=lambda_init)
            n_total += 1
            r65 = float(np.linalg.norm(o65.residual))
            rid = float(np.linalg.norm(oid.residual))
            tag65 = f"{'OK':>3} {o65.iters:>3}i {r65:>10.2e}" if o65.success else f"FAIL {o65.iters:>3}i {o65.message[:8]:>10}"
            tagid = f"{'OK':>3} {oid.iters:>3}i {rid:>10.2e}" if oid.success else f"FAIL {oid.iters:>3}i {oid.message[:8]:>10}"
            if o65.success and oid.success:
                n_both_success += 1
                d = float(np.max(np.abs(o65.x - oid.x)))
                same = "yes" if d < same_thresh else "no"
                if same == "yes":
                    n_same += 1
                print(f"  {v:>3}  {i+1:>4}  {clers[:20]:>20}  "
                      f"{tag65:>22}  {tagid:>22}  {same:>6}  {d:>9.2e}")
            else:
                print(f"  {v:>3}  {i+1:>4}  {clers[:20]:>20}  "
                      f"{tag65:>22}  {tagid:>22}  {'-':>6}  {'-':>9}")
    print()
    print(f"summary: total={n_total}  both_succeed={n_both_success}  "
          f"same_solution={n_same}  diff={n_both_success - n_same}")

if __name__ == "__main__":
    main()
