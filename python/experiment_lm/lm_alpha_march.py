#!/usr/bin/env python3
"""LM α-march: from α=0 ideal bends, head directly to a target α with LM.

No homotopy framework — just one LM call per target α, starting from the
horou-derived ideal bends. Dent protection: trial_gate rejects any LM
trial that introduces a dent (vertex_turn < 0 at a non-base vertex).
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


def build_setup(clers: str):
    netcode = decode(clers)
    faces = parse_netcode(netcode)
    tri = Tri.from_faces(faces)
    u = horou.horou(faces)
    bends = dh.all_dihedrals(faces, u)
    base_face = faces[0]
    base_edges = [edge_key(base_face[i], base_face[(i+1) % 3]) for i in range(3)]
    base_edge_set = set(base_edges)
    var_edges = [e for e in tri.edges if e not in base_edge_set]
    return tri, faces, base_face, base_edges, var_edges, bends


def residual_at(tri, base_face, var_edges, base_edges, bends, alpha):
    base_bend = {e: bends[e] for e in base_edges}
    def F(x):
        return holonomy_residual(tri, base_face, var_edges, x, base_bend, alpha)
    return F


def make_dent_gate(tri, base_face, var_edges, bends_full):
    base_face_set = set(base_face)
    def gate(x_trial):
        trial = dict(bends_full)
        for e, v in zip(var_edges, x_trial):
            trial[e] = float(v)
        for v in tri.vertices:
            if v in base_face_set:
                continue
            if vertex_turn(tri, v, trial) < 0.0:
                return False
        return True
    return gate


def march_one(clers: str, alpha_deg_list, max_iter=200, tol=1e-12,
              lambda_init=1e-3, start_kind: str = "ideal"):
    tri, faces, base_face, base_edges, var_edges, bends = build_setup(clers)
    V = len(tri.vertices)
    E = len(tri.edges)
    print(f"CLERS: {clers}")
    print(f"V={V} E={E} var_edges={len(var_edges)}  start={start_kind}")
    # Build starting bends per request. The trial_gate and base_bend still
    # use horou-derived `bends` for base edges; only var-edge x0 changes.
    if start_kind == "ideal":
        bends_for_x0 = bends
    elif start_kind == "zero":
        bends_for_x0 = {e: 0.0 for e in tri.edges}
    else:
        # numeric degrees, e.g. "30" → π/6 rad on every edge
        deg = float(start_kind)
        rad = math.radians(deg)
        bends_for_x0 = {e: rad for e in tri.edges}
    print()
    print(f"  {'α°':>7}  {'success':>8}  {'iters':>6}  {'resid':>11}  "
          f"{'cond_J':>10}  {'λ_min':>10}  {'λ_max':>10}  "
          f"{'retries_res':>11}  {'retries_dent':>12}  {'wall':>6}  msg")
    for alpha_deg in alpha_deg_list:
        alpha = math.radians(alpha_deg)
        F = residual_at(tri, base_face, var_edges, base_edges, bends, alpha)
        gate = make_dent_gate(tri, base_face, var_edges, bends_for_x0)
        x0 = np.array([bends_for_x0[e] for e in var_edges])
        log = []
        t0 = time.monotonic()
        out = solver_lm(F, x0, tol=tol, max_iter=max_iter,
                        lambda_init=lambda_init, iter_log=log,
                        trial_gate=gate)
        wall = time.monotonic() - t0
        # Aggregate stats from per-iter log.
        if log:
            cond_J = log[-1].get("cond_J", float("nan"))
            lambdas = [r["lambda"] for r in log]
            l_min = min(lambdas)
            l_max = max(lambdas)
            tot_res = sum(r.get("retries_residual", 0) for r in log)
            tot_dent = sum(r.get("retries_dent", 0) for r in log)
        else:
            cond_J = float("nan")
            l_min = l_max = tot_res = tot_dent = 0
        resid = float(np.linalg.norm(out.residual))
        print(f"  {alpha_deg:>7.2f}  {str(out.success):>8}  "
              f"{out.iters:>6}  {resid:>11.3e}  {cond_J:>10.2e}  "
              f"{l_min:>10.2e}  {l_max:>10.2e}  "
              f"{tot_res:>11d}  {tot_dent:>12d}  {wall:>5.2f}s  {out.message}")


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("clers")
    p.add_argument("--start", default="ideal",
                   help="ideal|zero|<deg>; e.g. --start 30 means uniform 30° on every edge")
    p.add_argument("alphas", nargs="*", type=float)
    args = p.parse_args()
    alphas = args.alphas if args.alphas else [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
    march_one(args.clers, alphas, start_kind=args.start)
