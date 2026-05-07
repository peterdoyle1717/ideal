#!/usr/bin/env python3
"""LM at one α (default first-homotopy-target = 2.5°) from various uniform
starts; compare each successful bend solution to the ideal-start solution.
"""
import json
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
           start_bends_full, max_iter=200, tol=1e-12, lambda_init=1e-3):
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
    return out, log

def main():
    if len(sys.argv) < 2:
        print("usage: lm_compare_starts.py <CLERS> [alpha_deg=2.5]", file=sys.stderr)
        sys.exit(2)
    clers = sys.argv[1]
    alpha_deg = float(sys.argv[2]) if len(sys.argv) > 2 else 2.5
    alpha = math.radians(alpha_deg)

    tri, base_face, base_edges, var_edges, bends_ideal = setup(clers)
    print(f"CLERS: {clers}  V={len(tri.vertices)}  E={len(tri.edges)}  "
          f"var_edges={len(var_edges)}  alpha={alpha_deg}°")
    print()

    # Reference: ideal start.
    out_ref, log_ref = run_lm(tri, base_face, base_edges, var_edges,
                                bends_ideal, alpha, bends_ideal)
    if not out_ref.success:
        print(f"REFERENCE (ideal start) FAILED: {out_ref.message}; aborting.")
        sys.exit(1)
    x_ref = out_ref.x
    print(f"reference (ideal start): success in {out_ref.iters} iters; "
          f"|r| = {float(np.linalg.norm(out_ref.residual)):.3e}")
    print()

    print(f"  {'start':>8}  {'success':>8}  {'iters':>6}  {'resid':>11}  "
          f"{'max|x-xref|':>14}  {'mean|x-xref|':>14}")
    starts_deg = list(range(0, 181, 10))
    for sd in starts_deg:
        sb = {e: math.radians(sd) for e in tri.edges}
        out, log = run_lm(tri, base_face, base_edges, var_edges,
                            bends_ideal, alpha, sb)
        if out.success:
            d = out.x - x_ref
            mx = float(np.max(np.abs(d)))
            mn = float(np.mean(np.abs(d)))
            r = float(np.linalg.norm(out.residual))
            print(f"  {sd:>5d}°    {str(out.success):>8}  {out.iters:>6}  "
                  f"{r:>11.3e}  {mx:>14.3e}  {mn:>14.3e}")
        else:
            print(f"  {sd:>5d}°    {str(out.success):>8}  {out.iters:>6}  "
                  f"{'(failed)':>11}  {'-':>14}  {'-':>14}    {out.message}")

if __name__ == "__main__":
    main()
