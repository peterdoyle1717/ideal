#!/usr/bin/env python3
"""LM at α=1° from a per-edge 'min-degree' start across small primes.

Start rule: bend(i,j) = 2π / min(deg(i), deg(j))
(So deg-3 endpoint → 120°, deg-4 → 90°, deg-5 → 72°, etc.)

Sweep lambda_init; compare to ideal-start solution.
"""
import math
import subprocess
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np

import horou
import dihedral as dh
from puffup import (
    parse_netcode, edge_key, Tri,
    holonomy_residual, vertex_turn,
    solver_lm,
)

CLERS_BIN = "/Users/doyle/Dropbox/neo/clers/bin/clers"
PRIMES_DIR = Path("/Users/doyle/Dropbox/neo/data/primes")


def decode(clers):
    return subprocess.run([CLERS_BIN, "decode"], input=clers + "\n",
                          capture_output=True, text=True, check=True).stdout.strip()


def setup(clers):
    nc = decode(clers)
    faces = parse_netcode(nc)
    tri = Tri.from_faces(faces)
    u = horou.horou(faces)
    bends_ideal = dh.all_dihedrals(faces, u)
    base_face = faces[0]
    base_edges = [edge_key(base_face[i], base_face[(i + 1) % 3]) for i in range(3)]
    var_edges = [e for e in tri.edges if e not in set(base_edges)]
    deg = Counter()
    for f in faces:
        for v in f:
            deg[v] += 1
    # Each face contributes once per vertex per face, so deg[v] = number
    # of faces around v = vertex degree (= number of edges at v).
    return tri, base_face, base_edges, var_edges, bends_ideal, deg


def mindeg_start(tri, deg):
    return {(i, j): 2.0 * math.pi / min(deg[i], deg[j]) for (i, j) in tri.edges}


def run_lm(tri, base_face, base_edges, var_edges, bends_ideal, alpha,
           start_bends_full, lambda_init, tol=1e-12, max_iter=300):
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
    return solver_lm(F, x0, tol=tol, max_iter=max_iter,
                    lambda_init=lambda_init, trial_gate=gate)


def main():
    v_min = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    v_max = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    alpha_deg = 1.0
    alpha = math.radians(alpha_deg)
    lambda_inits = [1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1e3, 1e4]
    same_thresh = 1e-6

    print(f"alpha={alpha_deg}°  start=mindeg(2π/min(deg(i),deg(j)) per edge)  v={v_min}..{v_max}")
    print(f"lambda_init grid: {lambda_inits}")
    print()
    hdr_li = "  ".join(f"{li:>9.0e}" for li in lambda_inits)
    print(f"  {'v':>3}  {'#':>3}  {'CLERS':>22}  {hdr_li}  {'note':>8}")
    print("  " + "-" * (8 + 22 + len(hdr_li) + 12))

    rows = []
    for v in range(v_min, v_max + 1):
        path = PRIMES_DIR / f"{v}.txt"
        if not path.exists():
            continue
        clersts = [c.strip() for c in path.read_text().splitlines() if c.strip()]
        for i, clers in enumerate(clersts):
            try:
                tri, bf, be, ve, bi, deg = setup(clers)
            except Exception as e:
                print(f"  {v:>3}  {i+1:>3}  {clers[:22]:>22}  setup fail: {e}")
                continue
            start = mindeg_start(tri, deg)
            o_ref = run_lm(tri, bf, be, ve, bi, alpha, bi, lambda_init=1e-3)
            x_ref = o_ref.x.copy() if o_ref.success else None
            cell = {}
            for li in lambda_inits:
                o = run_lm(tri, bf, be, ve, bi, alpha, start, lambda_init=li)
                if o.success and x_ref is not None:
                    d = float(np.max(np.abs(o.x - x_ref)))
                    cell[li] = ("=" if d < same_thresh else "≠", o.iters)
                elif o.success:
                    cell[li] = ("?", o.iters)
                else:
                    cell[li] = ("F", o.iters)
            rows.append((v, i + 1, clers, cell, "ref_ok" if x_ref is not None else "ref_fail"))
            cells = "  ".join(f"{cell[li][0]}{cell[li][1]:>7d}" for li in lambda_inits)
            print(f"  {v:>3}  {i+1:>3}  {clers[:22]:>22}  {cells}  {rows[-1][4]:>8}")

    print()
    print(f"  {'lambda_init':>13}  {'=ideal':>7}  {'≠ideal':>7}  {'? noref':>8}  {'fail':>5}")
    for li in lambda_inits:
        cnt = {"=": 0, "≠": 0, "?": 0, "F": 0}
        for (v, i, clers, cell, note) in rows:
            cnt[cell[li][0]] += 1
        print(f"  {li:>13.0e}  {cnt['=']:>7d}  {cnt['≠']:>7d}  {cnt['?']:>8d}  {cnt['F']:>5d}")


if __name__ == "__main__":
    main()
