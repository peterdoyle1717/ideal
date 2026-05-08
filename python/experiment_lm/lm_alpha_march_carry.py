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
try:
    from jacobian_sparse import analytical_jacobian_sparse
    _SPARSE_J_AVAILABLE = True
except ImportError:
    _SPARSE_J_AVAILABLE = False

try:
    from quat_residual import (
        holonomy_residual_quat, analytical_jacobian_quat_sparse,
        branch_check_post,
    )
    _QUAT_AVAILABLE = True
except ImportError:
    _QUAT_AVAILABLE = False


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
    # Sparse is the default: scales to large NVAR and matches dense to
    # machine precision on smaller cases. Pass --no-sparse to fall back
    # to the dense analytic Jacobian + dense linear solve (only useful
    # for A/B comparison or environments without scipy).
    ap.add_argument("--sparse",
                    action=argparse.BooleanOptionalAction,
                    default=True,
                    help="use scipy.sparse Jacobian + sparse LM solve "
                         "(default; --no-sparse for dense)")
    ap.add_argument("--residual", choices=("matrix", "quat"), default="matrix",
                    help="residual backend: matrix (default; 3 off-diagonals "
                         "of cumulative holomat) or quat (vector part of "
                         "cumulative holonomy quaternion, analytic Jacobian)")
    ap.add_argument("--system", choices=("square", "all-bends"), default="square",
                    help="solver formulation: square (default; var-edges are "
                         "the 3V-9 non-base edges, base bends frozen at horou "
                         "ideal, residuals at V-3 interior verts) or all-bends "
                         "(overdetermined: var-edges = all 3V-6 edges, no base "
                         "freeze, residuals at all V vertices, 3V × (3V-6))")
    ap.add_argument("--branch-vec-tol", type=float, default=1e-6,
                    help="quat post-check: a converged vertex has |vec(Q)| "
                         "below this AND qw should be negative (natural "
                         "-1 sheet for once-around loops); positive qw "
                         "with small vec is flagged as suspect")
    ap.add_argument("--alpha-jump-policy", choices=("half", "full"),
                    default="half",
                    help="α-march cap: half (default; cap = gap/2, then "
                         "halve on failure) or full (experiment; cap = "
                         "gap/1 = full jump to target, then halve on "
                         "failure to gap/2, gap/4, …)")
    ap.add_argument("--quick-min-drop", type=float, default=1e-2,
                    help="full-policy quick-reject: accept iff residual "
                         "≤ quick_tol OR final ≤ quick_min_drop × initial "
                         "(default 1e-2 = 100× shrink)")
    args = ap.parse_args()
    if args.residual == "quat" and not _QUAT_AVAILABLE:
        sys.exit("ERROR: --residual quat requires python/quat_residual.py")
    if args.system == "all-bends" and args.residual == "matrix":
        sys.exit("ERROR: --system all-bends currently requires --residual quat "
                 "(matrix backend not extended to all-bends in this commit)")

    netcode = decode(args.clers)
    faces = parse_netcode(netcode)
    tri = Tri.from_faces(faces)
    u = horou.horou(faces)
    bends_ideal = dh.all_dihedrals(faces, u)

    base_face = faces[0]
    base_edges = [edge_key(base_face[i], base_face[(i + 1) % 3]) for i in range(3)]
    base_edge_set = set(base_edges)
    base_face_set = set(base_face)

    # System layout: square (3V-9 vars, residuals at V-3 interior verts) vs
    # all-bends (3V-6 vars, residuals at all V verts, overdetermined by 6).
    if args.system == "all-bends":
        var_edges = list(tri.edges)               # ALL edges
        base_bend = {}                             # nothing frozen
        residual_vertices = list(tri.vertices)     # all V vertices
        gate_vertices = list(tri.vertices)         # gate over all V vertices
    else:
        var_edges = [e for e in tri.edges if e not in base_edge_set]
        base_bend = {e: bends_ideal[e] for e in base_edges}
        residual_vertices = [v for v in tri.vertices if v not in base_face_set]
        gate_vertices = list(residual_vertices)

    # Initial var-edge bends from horou ideal (covers ALL edges in either mode).
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
            for v_ in gate_vertices:
                if vertex_turn(tri, v_, tr) < 0.0:
                    return False
            return True
        return gate

    V = len(tri.vertices); E = len(tri.edges)
    n_unknowns = len(var_edges)
    n_residuals = 3 * len(residual_vertices)
    overdet = n_residuals - n_unknowns
    print(f"clers={args.clers}  V={V}  E={E}  system={args.system}  "
          f"residual={args.residual}")
    print(f"unknowns={n_unknowns}  residuals={n_residuals}  "
          f"overdet={overdet:+d}")
    if args.residual == "quat":
        # One-time J probe at α=0 with horou ideal bends — shape/nnz are
        # structural (independent of α); the values vary with α but the
        # sparsity pattern does not.
        try:
            J0 = analytical_jacobian_quat_sparse(
                tri, base_face, var_edges, x, base_bend, 0.0,
                vertices=residual_vertices)
            print(f"J shape={J0.shape}  nnz={J0.nnz}")
        except Exception as e:
            print(f"J probe failed: {e!r}")
    print(f"target={args.target_deg}°  init_step={args.init_step_deg}°  "
          f"min_step={args.min_step_deg}°  λ_init={args.lambda_init}  "
          f"max_quick_iter={args.max_quick_iter}  quick_tol={args.quick_tol}  "
          f"jump={args.alpha_jump_policy}")
    print()
    print(f"  {'#':>4}  {'α°':>9}  {'try°':>9}  {'frac':>6}  "
          f"{'iters':>5}  {'r_in':>10}  {'r_out':>10}  "
          f"{'λ_max':>10}  {'verdict':>7}")
    print(f"  {'-':>4}  {'-':>9}  {'-':>9}  {'-':>6}  {'-':>5}  "
          f"{'-':>10}  {'-':>10}  {'-':>10}  {'-':>7}")

    alpha_curr = 0.0
    n_step = 0
    n_retreat = 0
    n_attempts = 0
    final_resid = float("inf")
    history = []
    target_reached = False
    accepted_frac_seq = []
    rejected_frac_seq = []
    alpha_seq = []
    first_attempt_fraction = None
    accepted_first_fraction = None

    cap_fraction = 1.0 if args.alpha_jump_policy == "full" else 0.5
    try_fraction = cap_fraction
    # Legacy half-policy controller uses an absolute angular `step` that
    # halves on failure and gets capped at gap/2 each iteration. The
    # full-jump experiment uses a pure fraction-of-gap controller that
    # halves on failure and resets to cap_fraction (= 1.0) on accept.

    t_start = time.monotonic()
    for attempt in range(args.max_attempts):
        if alpha_curr >= alpha_final:
            target_reached = True
            break
        gap = target - alpha_curr
        if args.alpha_jump_policy == "full":
            used = gap * try_fraction
        else:
            cap = gap * 0.5
            used = min(step, cap)
            try_fraction = used / gap   # for stats
        alpha_try = alpha_curr + used

        bends_now = {e: float(b) for e, b in zip(var_edges, x)}
        bends_now.update(base_bend)
        gate = make_gate(bends_now)

        if args.residual == "quat":
            F = lambda xv, a=alpha_try: holonomy_residual_quat(
                tri, base_face, var_edges, xv, base_bend, a,
                vertices=residual_vertices)
            Jfn = lambda xv, a=alpha_try: analytical_jacobian_quat_sparse(
                tri, base_face, var_edges, xv, base_bend, a,
                vertices=residual_vertices)
        else:
            F = lambda xv, a=alpha_try: holonomy_residual(
                tri, base_face, var_edges, xv, base_bend, a)
            if args.sparse and _SPARSE_J_AVAILABLE:
                Jfn = lambda xv, a=alpha_try: analytical_jacobian_sparse(
                    tri, base_face, var_edges, xv, base_bend, a)
            else:
                Jfn = lambda xv, a=alpha_try: analytical_jacobian(
                    tri, base_face, var_edges, xv, base_bend, a)

        # Stale-bend residual at the new α — the LM's input residual.
        # Used for the full-policy "must shrink by quick_min_drop" test.
        try:
            initial_resid = float(np.linalg.norm(F(x)))
        except Exception:
            initial_resid = float("nan")

        log = []
        out = solver_lm(F, x, tol=args.quick_tol,
                        max_iter=args.max_quick_iter,
                        lambda_init=args.lambda_init,
                        iter_log=log, trial_gate=gate,
                        jacobian=Jfn)
        l_max = max((r["lambda"] for r in log), default=args.lambda_init)
        resid = float(np.linalg.norm(out.residual))
        n_step += 1
        n_attempts += 1
        if first_attempt_fraction is None:
            first_attempt_fraction = try_fraction

        # Acceptance:
        # half:  legacy — out.success AND resid < quick_tol.
        # full:  resid ≤ quick_tol OR resid ≤ quick_min_drop × initial_resid
        #        (fast-path for big jumps that drop a lot but not yet to
        #        the quick_tol floor; LM at the next α step takes them
        #        the rest of the way). No out.success requirement: a big
        #        jump that drops residual by ≥100× already counts as
        #        forward progress regardless of LM's tol-based "success".
        if args.alpha_jump_policy == "full":
            denom = max(initial_resid, 1e-30) if initial_resid == initial_resid else 1.0
            accepted_tol  = bool(resid <= args.quick_tol)
            accepted_drop = bool(resid <= args.quick_min_drop * denom)
            accepted = accepted_tol or accepted_drop
        else:
            accepted = bool(out.success and resid < args.quick_tol)
        verdict = "ok" if accepted else "retreat"
        print(f"  {n_step:>4}  {math.degrees(alpha_curr):>9.4f}  "
              f"{math.degrees(alpha_try):>9.4f}  {try_fraction:>6.3f}  "
              f"{out.iters:>5}  {initial_resid:>10.3e}  {resid:>10.3e}  "
              f"{l_max:>10.2e}  {verdict:>7}")
        history.append({
            "step": n_step,
            "alpha_curr_deg": math.degrees(alpha_curr),
            "alpha_try_deg":  math.degrees(alpha_try),
            "fraction": try_fraction,
            "iters":  int(out.iters),
            "initial_resid": initial_resid,
            "resid":  resid,
            "lambda_max": l_max,
            "verdict": verdict,
        })

        if accepted:
            x = out.x
            alpha_curr = alpha_try
            accepted_frac_seq.append(try_fraction)
            alpha_seq.append(math.degrees(alpha_curr))
            if accepted_first_fraction is None:
                accepted_first_fraction = try_fraction
            if args.alpha_jump_policy == "full":
                try_fraction = cap_fraction          # reset
            # else: legacy step does NOT grow on success.
        else:
            rejected_frac_seq.append(try_fraction)
            n_retreat += 1
            if args.alpha_jump_policy == "full":
                try_fraction *= 0.5
                if try_fraction * gap < min_step:
                    print(f"  FRACTION UNDERFLOW "
                          f"(fraction*gap < {args.min_step_deg}°) "
                          f"at α={math.degrees(alpha_curr):.6f}°")
                    break
            else:
                step *= 0.5
                if step < min_step:
                    print(f"  STEP UNDERFLOW "
                          f"(alpha_step < {args.min_step_deg}°) "
                          f"at α={math.degrees(alpha_curr):.6f}°")
                    break

    # Catch the edge case where the for-loop exhausts max_attempts on
    # the same iteration that reached alpha_final: the top-of-loop guard
    # never fired, but alpha_curr is at target. The final tight LM still
    # belongs.
    if alpha_curr >= alpha_final:
        target_reached = True

    if target_reached:
        # Final tight LM at exactly the target.
        if args.residual == "quat":
            F = lambda xv: holonomy_residual_quat(
                tri, base_face, var_edges, xv, base_bend, target,
                vertices=residual_vertices)
            Jfn = lambda xv: analytical_jacobian_quat_sparse(
                tri, base_face, var_edges, xv, base_bend, target,
                vertices=residual_vertices)
        else:
            F = lambda xv: holonomy_residual(
                tri, base_face, var_edges, xv, base_bend, target)
            if args.sparse and _SPARSE_J_AVAILABLE:
                Jfn = lambda xv: analytical_jacobian_sparse(
                    tri, base_face, var_edges, xv, base_bend, target)
            else:
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

        if args.residual == "quat":
            # Post-convergence sanity: each vertex should have |vec(Q)| ~0
            # AND qw on the −1 sheet (natural for once-around loops).
            # Suspect = converged but qw > 0.
            report = branch_check_post(tri, base_face, var_edges, x,
                                        base_bend, target,
                                        vec_tol=args.branch_vec_tol,
                                        vertices=residual_vertices)
            n_susp = sum(1 for r in report if r["suspect"])
            print()
            print(f"branch post-check (vec_tol={args.branch_vec_tol:.1e}):")
            print(f"  {'v':>4}  {'qw':>14}  {'|vec(Q)|':>12}  suspect")
            for r in report:
                flag = "SUSPECT" if r["suspect"] else ""
                print(f"  {r['vertex']:>4}  {r['qw']:>14.6e}  "
                      f"{r['vec_norm']:>12.3e}  {flag}")
            print(f"  total suspect vertices: {n_susp}")

        if args.bends_out:
            write_puffup_bends(args.bends_out, netcode, tri,
                               var_edges, base_bend, x, target)
            print(f"wrote bends -> {args.bends_out}")
    else:
        print()
        print(f"DID NOT REACH target. α reached: "
              f"{math.degrees(alpha_curr):.6f}°  retreats: {n_retreat}")

    wall = time.monotonic() - t_start
    if not target_reached:
        status = "STALLED_POSITIVE_RESIDUAL"
    elif final_resid <= args.final_tol:
        status = "SOLVER_TOL"
    else:
        msg = (out_final.message or "").lower()
        if "lambda" in msg:
            status = "LAMBDA_SATURATED_POSITIVE_RESIDUAL"
        else:
            status = "STALLED_POSITIVE_RESIDUAL"

    # Did the very first try go straight to alpha=60 and succeed?
    first_jump_to_target = False
    if (first_attempt_fraction is not None and
            abs(first_attempt_fraction - 1.0) < 1e-12 and
            accepted_first_fraction is not None and
            abs(accepted_first_fraction - 1.0) < 1e-12 and
            target_reached):
        first_jump_to_target = True

    # Retreats before the first accept (= rejects until first accept).
    retreats_before_first_accept = 0
    for h in history:
        if h["verdict"] == "ok":
            break
        retreats_before_first_accept += 1

    accepted_seq_str = ",".join(f"{f:.6g}" for f in accepted_frac_seq) or "-"
    rejected_seq_str = ",".join(f"{f:.6g}" for f in rejected_frac_seq) or "-"
    alpha_seq_str = ",".join(f"{a:.6g}" for a in alpha_seq) or "-"
    max_attempt_fraction = max([f for f in (accepted_frac_seq +
                                             rejected_frac_seq)] or [0.0])

    print(f"summary: clers={args.clers} target_reached={target_reached} "
          f"steps={n_step} retreats={n_retreat} attempts={n_attempts} "
          f"final_resid={final_resid:.6e} wall={wall:.2f}s "
          f"status={status} policy={args.alpha_jump_policy} "
          f"first_jump_to_target={first_jump_to_target} "
          f"retreats_before_first_accept={retreats_before_first_accept} "
          f"first_attempt_fraction={first_attempt_fraction!s} "
          f"accepted_first_fraction={accepted_first_fraction!s} "
          f"max_attempt_fraction={max_attempt_fraction:.6g}")
    print(f"stats: accepted_fraction_seq={accepted_seq_str}")
    print(f"stats: rejected_fraction_seq={rejected_seq_str}")
    print(f"stats: alpha_seq_deg={alpha_seq_str}")


if __name__ == "__main__":
    main()
