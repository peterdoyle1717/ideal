#!/usr/bin/env python3
"""Sweep starting bend value × target α; print compact success grid.

For each starting bend (uniform value or 'ideal' = horou-derived), runs LM
directly at each target α (no homotopy), with dent gating. Prints one row
per start, one column per α, marking cells with iter count or fail symbol.
"""
import json
import math
import subprocess
import sys
import time
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


def decode(clers: str) -> str:
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
    return tri, faces, base_face, base_edges, var_edges, bends_ideal


def run_one(tri, base_face, base_edges, var_edges, bends_ideal, alpha,
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
        print("usage: lm_start_grid.py <CLERS>", file=sys.stderr)
        sys.exit(2)
    clers = sys.argv[1]
    tri, faces, base_face, base_edges, var_edges, bends_ideal = setup(clers)
    print(f"CLERS: {clers}  V={len(tri.vertices)}  E={len(tri.edges)}  "
          f"var_edges={len(var_edges)}")
    alpha_degs = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
    # Uniform-deg starts plus the ideal start as a final row.
    uniform_starts_deg = list(range(0, 181, 10))
    starts = [("uniform", d) for d in uniform_starts_deg] + [("ideal", None)]

    # Header
    print()
    hdr = f"  {'start':>10}  " + " ".join(f"{d:>6.0f}°" for d in alpha_degs) + f"   ok/12"
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))

    for kind, val in starts:
        if kind == "uniform":
            start_bends = {e: math.radians(val) for e in tri.edges}
            label = f"{val:>4d}°    "
        else:
            start_bends = bends_ideal
            label = f"  ideal   "
        row_cells = []
        ok = 0
        for ad in alpha_degs:
            alpha = math.radians(ad)
            try:
                out, log = run_one(tri, base_face, base_edges, var_edges,
                                   bends_ideal, alpha, start_bends)
            except Exception as e:
                row_cells.append(f"  ERR ")
                continue
            if out.success:
                ok += 1
                row_cells.append(f"{out.iters:>6d}")
            else:
                # Encode failure mode briefly:
                #   d = dented saturation (most retries dent)
                #   r = residual saturation (most retries res)
                #   . = saturation neither dominant
                tot_d = sum(rec.get("retries_dent", 0) for rec in log)
                tot_r = sum(rec.get("retries_residual", 0) for rec in log)
                if tot_d > tot_r and tot_d > 0:
                    sym = "    -d"
                elif tot_r > 0:
                    sym = "    -r"
                else:
                    sym = "    --"
                row_cells.append(sym)
        print(f"  {label}  " + " ".join(row_cells) + f"   {ok:>2d}/{len(alpha_degs)}")


if __name__ == "__main__":
    main()
