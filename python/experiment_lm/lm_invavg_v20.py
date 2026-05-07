#!/usr/bin/env python3
"""LM with invavg start at α=1°, λ_init=1.0, across all primes v=4..20.

Per-CLERS: '=' (matches ideal), '≠' (other non-dented solution), or 'F'.
Per-v summary at end.

The user has clarified that any `≠ideal` solution is not a valid
embedding, so for this experiment we treat `≠` and `F` as both not-good
and only `=` as success.
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
    return tri, base_face, base_edges, var_edges, bends_ideal, deg


def invavg_start(tri, deg):
    return {(i, j): 0.5 * (2.0 * math.pi / deg[i] + 2.0 * math.pi / deg[j])
            for (i, j) in tri.edges}


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
    v_max = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    sample = int(sys.argv[3]) if len(sys.argv) > 3 else 0  # 0 = no cap
    alpha_deg = 1.0
    alpha = math.radians(alpha_deg)
    lambda_init = 1.0
    same_thresh = 1e-6

    sample_note = f"  sample=first {sample} per v" if sample > 0 else ""
    print(f"alpha={alpha_deg}°  start=invavg  lambda_init={lambda_init}  v={v_min}..{v_max}{sample_note}", flush=True)
    print(f"per-v summary at end. = matches ideal, ≠ different non-dented "
          f"solution (not embedded per user note), F lambda_saturated/etc.", flush=True)
    print(flush=True)

    by_v = {}
    for v in range(v_min, v_max + 1):
        path = PRIMES_DIR / f"{v}.txt"
        if not path.exists():
            continue
        clersts = [c.strip() for c in path.read_text().splitlines() if c.strip()]
        if sample > 0 and len(clersts) > sample:
            clersts = clersts[:sample]
        cnt = {"=": 0, "≠": 0, "F": 0, "ref_fail": 0}
        notable = []  # CLERSes that didn't =
        for i, clers in enumerate(clersts):
            try:
                tri, bf, be, ve, bi, deg = setup(clers)
            except Exception as e:
                cnt["F"] += 1
                continue
            o_ref = run_lm(tri, bf, be, ve, bi, alpha, bi, lambda_init=lambda_init)
            x_ref = o_ref.x.copy() if o_ref.success else None
            if x_ref is None:
                cnt["ref_fail"] += 1
                continue
            start = invavg_start(tri, deg)
            o = run_lm(tri, bf, be, ve, bi, alpha, start, lambda_init=lambda_init)
            if o.success:
                d = float(np.max(np.abs(o.x - x_ref)))
                if d < same_thresh:
                    cnt["="] += 1
                else:
                    cnt["≠"] += 1
                    notable.append((clers, "≠", d, o.iters))
            else:
                cnt["F"] += 1
                notable.append((clers, "F", float("nan"), o.iters))
        total = sum(cnt.values())
        by_v[v] = (cnt, total, notable)
        # Print per-v line as we go.
        print(f"  v={v:>2}  total={total:>5}  ={cnt['=']:>5}  "
              f"≠={cnt['≠']:>4}  F={cnt['F']:>4}  ref_fail={cnt['ref_fail']:>3}",
              flush=True)
        # Print up to 5 notable CLERS for this v.
        for (c, tag, d, iters) in notable[:5]:
            print(f"      {tag} {c}  d={d:.2e}  iters={iters}", flush=True)
        if len(notable) > 5:
            print(f"      ... and {len(notable) - 5} more not =", flush=True)

    # Aggregate.
    print(flush=True)
    g_eq = sum(b[0]["="] for b in by_v.values())
    g_ne = sum(b[0]["≠"] for b in by_v.values())
    g_f = sum(b[0]["F"] for b in by_v.values())
    g_rf = sum(b[0]["ref_fail"] for b in by_v.values())
    g_total = sum(b[1] for b in by_v.values())
    print(f"  TOTAL  v={v_min}..{v_max}  total={g_total}  ={g_eq}  ≠={g_ne}  "
          f"F={g_f}  ref_fail={g_rf}", flush=True)


if __name__ == "__main__":
    main()
