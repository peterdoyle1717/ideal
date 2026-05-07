#!/usr/bin/env python3
"""LM with fixed uniform start = 65°, sweep lambda_init across α targets.

Larger lambda_init = more conservative (smaller, more gradient-aligned)
first inner step. Smaller lambda_init = more Newton-like first step.
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
    return out, log, x0

def main():
    if len(sys.argv) < 2:
        print("usage: lm_lambda_sweep.py <CLERS>", file=sys.stderr)
        sys.exit(2)
    clers = sys.argv[1]
    tri, base_face, base_edges, var_edges, bends_ideal = setup(clers)
    print(f"CLERS: {clers}  V={len(tri.vertices)}  E={len(tri.edges)}  "
          f"var_edges={len(var_edges)}")
    start_deg = 65.0
    start_bends = {e: math.radians(start_deg) for e in tri.edges}

    alphas_deg = [2.5, 5, 15, 30, 45, 55, 60]
    lambda_inits = [1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0]

    # Reference: ideal start at each α (for comparing convergence point).
    refs = {}
    for ad in alphas_deg:
        out, _, _ = run_lm(tri, base_face, base_edges, var_edges,
                             bends_ideal, math.radians(ad), bends_ideal)
        refs[ad] = (out.x.copy() if out.success else None, out.success)

    print(f"\nstart = uniform {start_deg}° on every edge")
    print()
    print(f"  {'α°':>6}  {'λ_init':>9}  {'success':>8}  {'iters':>6}  "
          f"{'resid':>11}  {'max|x-xref|':>14}  {'ref?':>6}  msg")
    print("  " + "-" * 92)
    for ad in alphas_deg:
        x_ref, ref_ok = refs[ad]
        ref_label = "ok" if ref_ok else "fail"
        for li in lambda_inits:
            out, _, _ = run_lm(tri, base_face, base_edges, var_edges,
                                 bends_ideal, math.radians(ad), start_bends,
                                 lambda_init=li)
            if out.success and x_ref is not None:
                d = out.x - x_ref
                mx = float(np.max(np.abs(d)))
                same = "yes" if mx < 1e-6 else "no"
                resid = float(np.linalg.norm(out.residual))
                print(f"  {ad:>6.2f}  {li:>9.0e}  {str(out.success):>8}  "
                      f"{out.iters:>6}  {resid:>11.3e}  {mx:>14.3e}  "
                      f"{same:>6}  {out.message}")
            else:
                print(f"  {ad:>6.2f}  {li:>9.0e}  {str(out.success):>8}  "
                      f"{out.iters:>6}  {'':>11}  {'':>14}  "
                      f"{'-':>6}  {out.message}  (ref {ref_label})")
        print()

if __name__ == "__main__":
    main()
