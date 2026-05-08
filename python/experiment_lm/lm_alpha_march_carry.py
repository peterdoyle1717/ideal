#!/usr/bin/env python3
"""LM α-march with bend carry between α steps + halve-on-fail retreat.

Differs from sibling lm_alpha_march.py: that script does ONE LM call
per requested α, each starting from horou-derived ideal bends at α=0.
This one carries the converged bends from each successful α step into
the next, which is what 'homotopy with LM as corrector' actually means.

CLI:
  lm_alpha_march_carry.py CLERS \\
      [--target-deg 60] [--init-step-deg 1.0] [--min-step-deg 1e-4] \\
      [--max-quick-iter 10] [--quick-tol 1e-3] [--final-tol 1e-12] \\
      [--lambda-init 1.0] [--bends-out FILE] [--max-attempts 2000]

Output:
  Per-step row table to stdout.
  If --bends-out is set AND target α is reached, writes a
  `puffup-bends 1` file consumable by ~/Dropbox/neo/euclid_realize.
"""
import argparse
import json
import math
import os
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
    parse_netcode, edge_key, Tri,
    holonomy_residual, vertex_turn, solver_lm,
)
from jacobian import analytical_jacobian


def _find_clers_bin():
    """Locate the clers decoder. Order: $CLERS_BIN, then known sibling
    paths on laptop (Dropbox layout) and doob (~/neo/clers). Fail
    clearly rather than at first subprocess call if none is found."""
    env = os.environ.get("CLERS_BIN")
    if env and Path(env).is_file():
        return env
    home = Path.home()
    candidates = [
        home / "Dropbox/neo/clers/bin/clers",
        home / "neo/clers/bin/clers",
        home / "Dropbox/neo/orchestrator/tools/clers/bin/clers",
        home / "neo/orchestrator/tools/clers/bin/clers",
    ]
    for p in candidates:
        if p.is_file():
            return str(p)
    raise FileNotFoundError(
        "clers decoder not found; set $CLERS_BIN or build "
        "{Dropbox/,~/}/neo/clers")


CLERS_BIN = _find_clers_bin()


def decode(clers: str) -> str:
    return subprocess.run(
        [CLERS_BIN, "decode"], input=clers + "\n",
        capture_output=True, text=True, check=True,
    ).stdout.strip()


def canonical_edges_in_build_order(faces, nv):
    """Reproduce homotopy_stage puffup_c.c::build()'s edge enumeration:
    for each face in input order, nbr_add(a,b)/nbr_add(b,a)/nbr_add(b,c)/
    nbr_add(c,b)/nbr_add(a,c)/nbr_add(c,a) — first-seen wins. Then for
    u = 1..NV, walk NBR[u] in insertion order and emit (u, w) iff u < w.

    This is NOT pure (min,max)-lex order in general: for u, the
    neighbors are in face-scan / nbr_add insertion order, which depends
    on which face introduces each edge first. The two orderings happen
    to coincide for many small symmetric CLERSes (e.g. CCAE) but diverge
    on larger asymmetric ones. euclid_realize's parser strictly checks
    row i against the canonical EDGE_A[i]/EDGE_B[i] of its own build, so
    the writer must match this order exactly — codex caught a v=26
    case where pure lex disagrees.
    """
    nbr = {v: [] for v in range(1, nv + 1)}
    seen = set()

    def nbr_add(u, w):
        if (u, w) not in seen:
            nbr[u].append(w)
            seen.add((u, w))

    for (a, b, c) in faces:
        # exact same call order as homotopy_stage build() line ~175
        nbr_add(a, b); nbr_add(b, a)
        nbr_add(b, c); nbr_add(c, b)
        nbr_add(a, c); nbr_add(c, a)

    edges = []
    for u in range(1, nv + 1):
        for w in nbr[u]:
            if u < w:
                edges.append((u, w))
    return edges


def write_puffup_bends(path, netcode, tri, var_edges, base_bend, x_var,
                       alpha_rad):
    """Write bends in the canonical homotopy_stage EDGE_A[i]/EDGE_B[i]
    order so euclid_realize's parser accepts the file."""
    bend_full = dict(base_bend)
    for e, v in zip(var_edges, x_var):
        bend_full[e] = float(v)
    nv = len(tri.vertices)
    edges = canonical_edges_in_build_order(tri.faces, nv)
    if len(edges) != len(tri.edges):
        raise RuntimeError(
            f"canonical-edge enumeration produced {len(edges)} edges, "
            f"tri.edges has {len(tri.edges)}; topology mismatch?")
    with open(path, "w") as fh:
        fh.write("puffup-bends 1\n")
        fh.write(f"NV {nv} NE {len(edges)} "
                 f"alpha_deg {math.degrees(alpha_rad):.15g}\n")
        fh.write(f"faces {netcode}\n")
        fh.write("bends\n")
        for (a, b) in edges:
            # bend dict keyed by edge_key (which is (min,max));
            # (a, b) here already has a < b since we filter u < w above.
            fh.write(f"{a} {b} {bend_full[(a, b)]:.15g}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("clers")
    ap.add_argument("--target-deg",    type=float, default=60.0)
    ap.add_argument("--init-step-deg", type=float, default=1.0)
    ap.add_argument("--min-step-deg",  type=float, default=1e-4)
    ap.add_argument("--max-quick-iter", type=int,  default=10)
    ap.add_argument("--quick-tol",     type=float, default=1e-3)
    ap.add_argument("--final-tol",     type=float, default=1e-12)
    ap.add_argument("--lambda-init",   type=float, default=1.0)
    ap.add_argument("--bends-out",     type=str,   default=None)
    ap.add_argument("--max-attempts",  type=int,   default=2000)
    args = ap.parse_args()

    netcode = decode(args.clers)
    faces = parse_netcode(netcode)
    tri = Tri.from_faces(faces)
    u = horou.horou(faces)
    bends_ideal = dh.all_dihedrals(faces, u)

    base_face = faces[0]
    base_edges = [edge_key(base_face[i], base_face[(i + 1) % 3]) for i in range(3)]
    base_edge_set = set(base_edges)
    var_edges = [e for e in tri.edges if e not in base_edge_set]
    base_bend = {e: bends_ideal[e] for e in base_edges}
    base_face_set = set(base_face)

    # Initial var-edge bends from horou ideal.
    x = np.array([bends_ideal[e] for e in var_edges], dtype=float)

    target = math.radians(args.target_deg)
    step   = math.radians(args.init_step_deg)
    min_step = math.radians(args.min_step_deg)
    # Match homotopy_stage's puffup_c controller: ALPHA_FINAL is target *
    # (1 - 1e-12). Loop exits when alpha_curr crosses that, which it
    # eventually does because the halfway-cap halves remaining-gap each
    # accepted step (geometric → reaches the (1-1e-12) threshold in ~40
    # halvings without further failures).
    alpha_final = target * (1.0 - 1e-12)

    def make_gate(bends_curr_dict):
        def gate(x_trial):
            tr = dict(bends_curr_dict)
            for e, v in zip(var_edges, x_trial):
                tr[e] = float(v)
            for v_ in tri.vertices:
                if v_ in base_face_set:
                    continue
                if vertex_turn(tri, v_, tr) < 0.0:
                    return False
            return True
        return gate

    print(f"clers={args.clers}  V={len(tri.vertices)}  E={len(tri.edges)}  "
          f"var={len(var_edges)}")
    print(f"target={args.target_deg}°  init_step={args.init_step_deg}°  "
          f"min_step={args.min_step_deg}°  λ_init={args.lambda_init}  "
          f"max_quick_iter={args.max_quick_iter}  quick_tol={args.quick_tol}")
    print()
    print(f"  {'#':>4}  {'α°':>9}  {'try°':>9}  {'iters':>5}  "
          f"{'resid':>11}  {'λ_max':>10}  {'verdict':>7}")
    print(f"  {'-':>4}  {'-':>9}  {'-':>9}  {'-':>5}  "
          f"{'-':>11}  {'-':>10}  {'-':>7}")

    alpha_curr = 0.0
    n_step = 0
    n_retreat = 0
    final_resid = float("inf")
    history = []
    target_reached = False

    t_start = time.monotonic()
    for attempt in range(args.max_attempts):
        if alpha_curr >= alpha_final:
            target_reached = True
            break
        cap = (target - alpha_curr) * 0.5
        used = min(step, cap)
        alpha_try = alpha_curr + used

        bends_now = {e: float(b) for e, b in zip(var_edges, x)}
        bends_now.update(base_bend)
        gate = make_gate(bends_now)

        F = lambda xv, a=alpha_try: holonomy_residual(
            tri, base_face, var_edges, xv, base_bend, a)
        Jfn = lambda xv, a=alpha_try: analytical_jacobian(
            tri, base_face, var_edges, xv, base_bend, a)

        log = []
        out = solver_lm(F, x, tol=args.quick_tol,
                        max_iter=args.max_quick_iter,
                        lambda_init=args.lambda_init,
                        iter_log=log, trial_gate=gate,
                        jacobian=Jfn)
        l_max = max((r["lambda"] for r in log), default=args.lambda_init)
        resid = float(np.linalg.norm(out.residual))
        n_step += 1
        # accept iff LM said success AND residual is below quick_tol —
        # LM declares success as soon as |r| <= tol within max_iter,
        # but if it bailed at max_iter it may still be 'success=False'.
        accepted = bool(out.success and resid < args.quick_tol)
        verdict = "ok" if accepted else "halve"
        print(f"  {n_step:>4}  {math.degrees(alpha_curr):>9.4f}  "
              f"{math.degrees(alpha_try):>9.4f}  {out.iters:>5}  "
              f"{resid:>11.3e}  {l_max:>10.2e}  {verdict:>7}")
        history.append({
            "step": n_step,
            "alpha_curr_deg": math.degrees(alpha_curr),
            "alpha_try_deg":  math.degrees(alpha_try),
            "iters":  int(out.iters),
            "resid":  resid,
            "lambda_max": l_max,
            "verdict": verdict,
        })

        if accepted:
            x = out.x
            alpha_curr = alpha_try
            # step does NOT grow on success — match canonical controller.
        else:
            step *= 0.5
            n_retreat += 1
            if step < min_step:
                print(f"  STEP UNDERFLOW (alpha_step < {args.min_step_deg}°) "
                      f"at α={math.degrees(alpha_curr):.6f}°")
                break

    if target_reached:
        # Final tight LM at exactly the target.
        F = lambda xv: holonomy_residual(
            tri, base_face, var_edges, xv, base_bend, target)
        Jfn = lambda xv: analytical_jacobian(
            tri, base_face, var_edges, xv, base_bend, target)
        out_final = solver_lm(F, x, tol=args.final_tol,
                              max_iter=200,
                              lambda_init=args.lambda_init,
                              jacobian=Jfn)
        x = out_final.x
        alpha_curr = target
        final_resid = float(np.linalg.norm(out_final.residual))
        print()
        print(f"FINAL: α={math.degrees(alpha_curr):.10f}°  "
              f"iters={out_final.iters}  resid={final_resid:.6e}  "
              f"{out_final.message}")

        if args.bends_out:
            write_puffup_bends(args.bends_out, netcode, tri,
                               var_edges, base_bend, x, target)
            print(f"wrote bends -> {args.bends_out}")
    else:
        print()
        print(f"DID NOT REACH target. α reached: "
              f"{math.degrees(alpha_curr):.6f}°  retreats: {n_retreat}")

    wall = time.monotonic() - t_start
    print(f"summary: clers={args.clers} target_reached={target_reached} "
          f"steps={n_step} retreats={n_retreat} "
          f"final_resid={final_resid:.6e} wall={wall:.2f}s")


if __name__ == "__main__":
    main()
