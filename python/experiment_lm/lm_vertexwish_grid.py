#!/usr/bin/env python3
"""LM at α=1° from vertexwish start across small primes.

The vertexwish start rule projects invavg onto the per-vertex sum-to-1
constraint. Sweeps lambda_init; compares to ideal-start solution.
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

from experiment_lm.vertexwish import vertexwish_start_radians

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
    return tri, base_face, base_edges, var_edges, bends_ideal, deg


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

    print(f"alpha={alpha_deg}°  start=vertexwish  v={v_min}..{v_max}", flush=True)
    print(f"lambda_init grid: {lambda_inits}", flush=True)
    print(flush=True)
    hdr_li = "  ".join(f"{li:>9.0e}" for li in lambda_inits)
    print(f"  {'v':>3}  {'#':>3}  {'CLERS':>22}  {hdr_li}  "
          f"{'min_x':>9}  {'max_x':>9}  {'note':>8}", flush=True)
    print("  " + "-" * (8 + 22 + len(hdr_li) + 25), flush=True)

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
                print(f"  {v:>3}  {i+1:>3}  {clers[:22]:>22}  setup fail: {e}",
                      flush=True)
                continue
            try:
                start = vertexwish_start_radians(tri, deg)
            except Exception as e:
                print(f"  {v:>3}  {i+1:>3}  {clers[:22]:>22}  start fail: {e}",
                      flush=True)
                continue
            sxs = list(start.values())
            min_x_rev = min(sxs) / (2.0 * math.pi)
            max_x_rev = max(sxs) / (2.0 * math.pi)

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
            rows.append((v, i + 1, clers, cell,
                         "ref_ok" if x_ref is not None else "ref_fail"))
            cells = "  ".join(f"{cell[li][0]}{cell[li][1]:>7d}" for li in lambda_inits)
            print(f"  {v:>3}  {i+1:>3}  {clers[:22]:>22}  {cells}  "
                  f"{min_x_rev:>9.3f}  {max_x_rev:>9.3f}  {rows[-1][4]:>8}",
                  flush=True)

    print(flush=True)
    print("  per-lambda totals (over rows where setup AND start succeeded):",
          flush=True)
    print(f"  {'lambda_init':>13}  {'=ideal':>7}  {'≠ideal':>7}  "
          f"{'? noref':>8}  {'fail':>5}", flush=True)
    for li in lambda_inits:
        cnt = {"=": 0, "≠": 0, "?": 0, "F": 0}
        for (v, i, clers, cell, note) in rows:
            cnt[cell[li][0]] += 1
        print(f"  {li:>13.0e}  {cnt['=']:>7d}  {cnt['≠']:>7d}  "
              f"{cnt['?']:>8d}  {cnt['F']:>5d}", flush=True)


if __name__ == "__main__":
    main()
