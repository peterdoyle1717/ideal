#!/usr/bin/env python3
"""LM with vertexwish start at α=1°, lambda_init=1.0, across primes v=4..20,
plus an explicit v=26 hard case at the end.

Per-CLERS: '=' (matches ideal), '≠' (other non-dented solution), or 'F'.
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

V26_HARD = "CCCCACCACCACACCACACACACCACAACCACACAACACCADEABABE"


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
           start_bends_full, lambda_init=1.0, tol=1e-12, max_iter=300):
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
    v_max = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    alpha_deg = 1.0
    alpha = math.radians(alpha_deg)
    lambda_init = 1.0
    same_thresh = 1e-6

    print(f"alpha={alpha_deg}°  start=vertexwish  lambda_init={lambda_init}  "
          f"v={v_min}..{v_max}", flush=True)
    print(flush=True)

    by_v = {}
    for v in range(v_min, v_max + 1):
        path = PRIMES_DIR / f"{v}.txt"
        if not path.exists():
            continue
        clersts = [c.strip() for c in path.read_text().splitlines() if c.strip()]
        # Mutually-exclusive vertexwish classes: =, ≠, ?, F. ref_fail
        # is a separate diagnostic and is NOT included in `total`.
        cnt = {"=": 0, "≠": 0, "?": 0, "F": 0, "ref_fail": 0}
        notable = []
        for i, clers in enumerate(clersts):
            try:
                tri, bf, be, ve, bi, deg = setup(clers)
            except Exception as e:
                cnt["F"] += 1
                continue
            o_ref = run_lm(tri, bf, be, ve, bi, alpha, bi, lambda_init=lambda_init)
            x_ref = o_ref.x.copy() if o_ref.success else None
            if x_ref is None:
                # Diagnostic only — does not affect the four mutually-
                # exclusive vertexwish classes below.
                cnt["ref_fail"] += 1
            try:
                start = vertexwish_start_radians(tri, deg)
            except Exception as e:
                cnt["F"] += 1
                notable.append((clers, "F", float("nan"), 0))
                continue
            o = run_lm(tri, bf, be, ve, bi, alpha, start, lambda_init=lambda_init)
            if o.success and x_ref is not None:
                d = float(np.max(np.abs(o.x - x_ref)))
                if d < same_thresh:
                    cnt["="] += 1
                else:
                    cnt["≠"] += 1
                    notable.append((clers, "≠", d, o.iters))
            elif o.success:
                cnt["?"] += 1
                notable.append((clers, "?", float("nan"), o.iters))
            else:
                cnt["F"] += 1
                notable.append((clers, "F", float("nan"), o.iters))
        # `total` counts vertexwish-attempted cases (=, ≠, ?, F).
        total = cnt["="] + cnt["≠"] + cnt["?"] + cnt["F"]
        by_v[v] = (cnt, total, notable)
        print(f"  v={v:>2}  total={total:>5}  ={cnt['=']:>5}  "
              f"≠={cnt['≠']:>4}  ?={cnt['?']:>3}  F={cnt['F']:>4}  "
              f"ref_fail={cnt['ref_fail']:>3}",
              flush=True)
        for (c, tag, d, iters) in notable[:5]:
            print(f"      {tag} {c}  d={d:.2e}  iters={iters}", flush=True)
        if len(notable) > 5:
            print(f"      ... and {len(notable) - 5} more not =", flush=True)

    print(flush=True)
    g_eq = sum(b[0]["="] for b in by_v.values())
    g_ne = sum(b[0]["≠"] for b in by_v.values())
    g_q = sum(b[0]["?"] for b in by_v.values())
    g_f = sum(b[0]["F"] for b in by_v.values())
    g_rf = sum(b[0]["ref_fail"] for b in by_v.values())
    g_total = sum(b[1] for b in by_v.values())
    print(f"  TOTAL  v={v_min}..{v_max}  total={g_total}  ={g_eq}  ≠={g_ne}  "
          f"?={g_q}  F={g_f}  ref_fail={g_rf}", flush=True)

    print(flush=True)
    print(f"=== v=26 hard case: {V26_HARD} ===", flush=True)
    try:
        tri, bf, be, ve, bi, deg = setup(V26_HARD)
        # Reference: ideal start.
        o_ref = run_lm(tri, bf, be, ve, bi, alpha, bi, lambda_init=lambda_init)
        x_ref = o_ref.x.copy() if o_ref.success else None
        # vertexwish start.
        start = vertexwish_start_radians(tri, deg)
        sxs = list(start.values())
        min_deg = math.degrees(min(sxs))
        max_deg = math.degrees(max(sxs))
        print(f"  V={len(tri.vertices)} E={len(tri.edges)} "
              f"vertexwish bend range = {min_deg:.2f}° .. {max_deg:.2f}°",
              flush=True)
        print(f"  ref (ideal start): success={o_ref.success} iters={o_ref.iters} "
              f"resid={float(np.linalg.norm(o_ref.residual)):.3e}", flush=True)
        # Sweep λ_init for vertexwish on v=26.
        for li in [1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1e3, 1e4]:
            o = run_lm(tri, bf, be, ve, bi, alpha, start, lambda_init=li)
            if o.success and x_ref is not None:
                d = float(np.max(np.abs(o.x - x_ref)))
                tag = "=" if d < same_thresh else "≠"
                print(f"  λ_init={li:>9.0e}  {tag}  iters={o.iters}  "
                      f"resid={float(np.linalg.norm(o.residual)):.2e}  "
                      f"max|Δ|={d:.2e}", flush=True)
            elif o.success:
                # Vertexwish converged; reference unavailable for comparison.
                print(f"  λ_init={li:>9.0e}  ?  iters={o.iters}  "
                      f"resid={float(np.linalg.norm(o.residual)):.2e}  "
                      f"(noref)", flush=True)
            else:
                print(f"  λ_init={li:>9.0e}  F  iters={o.iters}  "
                      f"msg={o.message}", flush=True)
    except Exception as e:
        print(f"  v=26 hard case setup/run fail: {e}", flush=True)


if __name__ == "__main__":
    main()
