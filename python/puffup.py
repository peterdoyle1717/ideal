#!/usr/bin/env python3
"""
puffup.py — neoplatonic dihedral solver via homotopy continuation in
the corner angle.

Direct port of the working plat1000.wl approach:
  - Variables: bend angles (radians) on each edge.
  - Geometric primitive: 3x3 rotations matz(yaw) @ matx(-roll).
  - Vertex holonomy: walk the face cycle around a vertex, compose
    movemats; the final matrix should be identity. Residual = three
    off-diagonal entries (M[0,1], M[0,2], M[1,2]).
  - Square system: 3v-9 non-base bends as variables vs 3v-9 holonomy
    residuals from non-base vertices.
  - Homotopy: ramp corner angle from `from_corner` (default 0, the ideal
    horoball limit) to the target. Adaptive step: halve on turn-check
    failure.

Input/output: JSON. See README at top of file or run with --help.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import defaultdict, deque
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

Vertex = int
Face = Tuple[Vertex, Vertex, Vertex]
Edge = Tuple[Vertex, Vertex]


# ── parsing ──────────────────────────────────────────────────────────────

def edge_key(a: Vertex, b: Vertex) -> Edge:
    if a == b:
        raise ValueError("loop edge")
    return (a, b) if a < b else (b, a)


def edge_str(e: Edge) -> str:
    return f"{e[0]}-{e[1]}"


def parse_edge_str(s: str) -> Edge:
    a, b = s.split("-")
    return edge_key(int(a), int(b))


def parse_netcode(text: str) -> List[Face]:
    faces: List[Face] = []
    for raw in text.replace("\n", ";").split(";"):
        raw = raw.strip()
        if not raw:
            continue
        parts = [p.strip() for p in raw.split(",") if p.strip()]
        if len(parts) != 3:
            raise ValueError(f"bad face {raw!r}; expected a,b,c")
        tri = tuple(int(p) for p in parts)
        if len(set(tri)) != 3:
            raise ValueError(f"degenerate face {tri}")
        faces.append(tri)  # type: ignore[arg-type]
    if not faces:
        raise ValueError("empty netcode")
    return faces


# ── triangulation indices ────────────────────────────────────────────────

@dataclass
class Tri:
    faces: List[Face]
    vertices: List[Vertex]
    edges: List[Edge]
    directed_edge_face: Dict[Tuple[Vertex, Vertex], int]
    vertex_flower: Dict[Vertex, List[int]]

    @staticmethod
    def from_faces(faces: List[Face]) -> "Tri":
        verts = sorted({v for f in faces for v in f})
        directed: Dict[Tuple[Vertex, Vertex], int] = {}
        edge_count: Dict[Edge, int] = defaultdict(int)
        for fi, f in enumerate(faces):
            for i in range(3):
                a, b = f[i], f[(i + 1) % 3]
                if (a, b) in directed:
                    raise ValueError(f"directed edge {(a, b)} appears twice; orientation inconsistent")
                directed[(a, b)] = fi
                edge_count[edge_key(a, b)] += 1
        bad = {e: c for e, c in edge_count.items() if c != 2}
        if bad:
            raise ValueError(f"not a closed surface; bad edge incidences: {bad}")
        edges = sorted(edge_count)
        flowers = {v: _flower(v, faces, directed) for v in verts}
        return Tri(faces, verts, edges, directed, flowers)


def _flower(v: Vertex, faces: List[Face], directed: Dict[Tuple[Vertex, Vertex], int]) -> List[int]:
    """Cyclic list of face indices around vertex v, oriented to match plat1000.wl.

    From a face (v, b, c) we step to the face containing the directed edge (v, c).
    """
    start = next(fi for fi, f in enumerate(faces) if v in f)
    cycle = [start]
    seen = {start}
    while True:
        f = faces[cycle[-1]]
        i = f.index(v)
        c = f[(i + 2) % 3]  # third vertex (after v, b)
        nxt = directed.get((v, c))
        if nxt is None:
            raise ValueError(f"can't walk around vertex {v}: missing directed edge ({v}, {c})")
        if nxt == start:
            return cycle
        if nxt in seen:
            raise ValueError(f"flower at vertex {v} closes early")
        cycle.append(nxt)
        seen.add(nxt)


def face_with_v_first(face: Face, v: Vertex) -> Face:
    i = face.index(v)
    return (face[i], face[(i + 1) % 3], face[(i + 2) % 3])


# ── 3D rotations ─────────────────────────────────────────────────────────

def matz(a: float) -> np.ndarray:
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def matx(a: float) -> np.ndarray:
    c, s = math.cos(a), math.sin(a)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])


def movemat(yaw: float, roll: float) -> np.ndarray:
    """matz(yaw) @ matx(-roll), per plat1000.wl."""
    return matz(yaw) @ matx(-roll)


# ── holonomy residuals ───────────────────────────────────────────────────

def vertex_holomat(tri: Tri, v: Vertex, alpha: float, bend: Dict[Edge, float]) -> np.ndarray:
    M = np.eye(3)
    for fi in tri.vertex_flower[v]:
        f = face_with_v_first(tri.faces[fi], v)  # (v, b, c)
        c = f[2]
        roll = bend[edge_key(v, c)]
        M = M @ movemat(alpha, roll)
    return M


def take3(M: np.ndarray) -> np.ndarray:
    """Three off-diagonals that vanish iff M = I (linearization of SO(3))."""
    return np.array([M[0, 1], M[0, 2], M[1, 2]])


def vertex_residual(tri: Tri, v: Vertex, alpha: float, bend: Dict[Edge, float]) -> np.ndarray:
    return take3(vertex_holomat(tri, v, alpha, bend))


def holonomy_residual(
    tri: Tri,
    base_face: Face,
    var_edges: Sequence[Edge],
    var_values: np.ndarray,
    fixed_bend: Dict[Edge, float],
    alpha: float,
) -> np.ndarray:
    bend = dict(fixed_bend)
    for e, val in zip(var_edges, var_values):
        bend[e] = float(val)
    interior_vs = [v for v in tri.vertices if v not in base_face]
    parts = [vertex_residual(tri, v, alpha, bend) for v in interior_vs]
    return np.concatenate(parts) if parts else np.zeros(0)


def vertex_turn(tri: Tri, v: Vertex, bend: Dict[Edge, float]) -> float:
    seen: set = set()
    total = 0.0
    for fi in tri.vertex_flower[v]:
        f = face_with_v_first(tri.faces[fi], v)
        e = edge_key(v, f[2])
        if e not in seen:
            total += bend[e]
            seen.add(e)
    return total


# ── Newton solver ────────────────────────────────────────────────────────

def fd_jacobian(fun, x: np.ndarray, h: float = 1e-7) -> np.ndarray:
    n = len(x)
    f0 = fun(x)
    J = np.zeros((len(f0), n))
    for j in range(n):
        step = h * max(1.0, abs(float(x[j])))
        xp = x.copy(); xp[j] += step
        xm = x.copy(); xm[j] -= step
        J[:, j] = (fun(xp) - fun(xm)) / (2.0 * step)
    return J


@dataclass
class NewtonOut:
    x: np.ndarray
    residual: np.ndarray
    success: bool
    iters: int
    message: str


def newton(fun, x0: np.ndarray, tol: float = 1e-12, max_iter: int = 60) -> NewtonOut:
    x = np.array(x0, float)
    r = fun(x)
    norm = float(np.linalg.norm(r))
    for it in range(max_iter):
        if norm <= tol:
            return NewtonOut(x, r, True, it, "tol")
        J = fd_jacobian(fun, x)
        try:
            dx, *_ = np.linalg.lstsq(J, -r, rcond=None)
        except np.linalg.LinAlgError as exc:
            return NewtonOut(x, r, False, it, f"lstsq fail: {exc}")
        lam = 1.0
        improved = False
        for _ in range(40):
            r_try = fun(x + lam * dx)
            n_try = float(np.linalg.norm(r_try))
            if n_try < norm:
                x = x + lam * dx
                r = r_try
                norm = n_try
                improved = True
                break
            lam *= 0.5
        if not improved:
            return NewtonOut(x, r, False, it, "line search fail")
    return NewtonOut(x, r, norm <= tol, max_iter, "max iter")


# ── dent index (port of undented/src/dent_check.c) ──────────────────────

def vertex_neighbors_cyclic(tri: Tri, v: Vertex) -> List[Vertex]:
    """Cyclic list of v's neighbors, in flower order."""
    out: List[Vertex] = []
    seen: set = set()
    for fi in tri.vertex_flower[v]:
        f = face_with_v_first(tri.faces[fi], v)  # (v, b, c)
        if f[1] not in seen:
            out.append(f[1]); seen.add(f[1])
        if f[2] not in seen:
            out.append(f[2]); seen.add(f[2])
    return out


def link_turning(v: Vertex, neighbors: Sequence[Vertex],
                 coords: Dict[Vertex, np.ndarray]) -> float:
    """Sum of signed spherical exterior angles of the link polygon at v.
    Positive = undented at v. Direct port of dent_check.c::dent_index inner loop."""
    pv = coords[v]
    dirs = []
    for nb in neighbors:
        d = coords[nb] - pv
        n = float(np.linalg.norm(d))
        if n < 1e-15:
            n = 1e-15
        dirs.append(d / n)
    k = len(dirs)
    total = 0.0
    for i in range(k):
        A = dirs[(i - 1) % k]
        B = dirs[i]
        C = dirs[(i + 1) % k]
        cross = np.cross(A, C)
        num = float(B @ cross)
        den = float((A @ B) * (B @ C) - (A @ C))
        total += math.atan2(num, den)
    return total


def dent_index(tri: Tri, coords: Dict[Vertex, np.ndarray]) -> float:
    """Minimum link turning across all vertices. Positive = undented.

    The cyclic neighbor order comes from the face orientation in the input
    netcode. If you get a sign you didn't expect, the input orientation is
    wrong — fix it upstream, do NOT silently flip here."""
    min_turn = float("inf")
    for v in tri.vertices:
        nbrs = vertex_neighbors_cyclic(tri, v)
        if len(nbrs) < 3:
            continue
        t = link_turning(v, nbrs, coords)
        if t < min_turn:
            min_turn = t
    return min_turn


# ── homotopy driver ──────────────────────────────────────────────────────

@dataclass
class HomotopyOut:
    bend_full: Dict[Edge, float]
    final_residual: float
    success: bool
    steps: List[dict]
    iters_total: int
    message: str


def solve_homotopy(
    tri: Tri,
    base_face: Face,
    target_alpha: float,
    initial_bends: Dict[Edge, float],
    from_alpha: float = 0.0,
    step_factor: float = 1.0 / 24.0,
    min_step_factor: float = 1e-6,
    turn_eps: float = -1e-9,
    use_dent_check: bool = True,
    dent_eps: float = 1e-9,
    grow_factor: float = 1.5,            # multiply step on K consecutive successes
    grow_after: int = 3,                  # K
    max_step_factor: float = 0.5,         # cap on grown step
) -> HomotopyOut:
    var_edges = [e for e in tri.edges
                 if e != edge_key(base_face[0], base_face[1])
                 and e != edge_key(base_face[1], base_face[2])
                 and e != edge_key(base_face[2], base_face[0])]
    base_edges = [edge_key(base_face[0], base_face[1]),
                  edge_key(base_face[1], base_face[2]),
                  edge_key(base_face[2], base_face[0])]
    bend = dict(initial_bends)

    span = target_alpha - from_alpha
    if abs(span) < 1e-15:
        # Single solve, no homotopy needed
        return _final_solve(tri, base_face, target_alpha, bend, var_edges, base_edges)

    t = 0.0
    step = step_factor
    successes_in_a_row = 0
    iters_total = 0
    steps: List[dict] = []
    steps.append({"t": 0.0, "alpha": from_alpha,
                  "residual": float(np.linalg.norm(holonomy_residual(
                      tri, base_face, var_edges,
                      np.array([bend[e] for e in var_edges]),
                      {e: bend[e] for e in base_edges}, from_alpha)))})

    while t < 1.0:
        t1 = min(1.0, t + step)
        alpha1 = from_alpha + t1 * span
        x_init = np.array([bend[e] for e in var_edges])

        def F(x: np.ndarray) -> np.ndarray:
            return holonomy_residual(tri, base_face, var_edges, x,
                                     {e: bend[e] for e in base_edges}, alpha1)

        out = newton(F, x_init)
        iters_total += out.iters

        # Provisional bend dict with the trial solution
        trial_bend = dict(bend)
        for e, v in zip(var_edges, out.x):
            trial_bend[e] = float(v)

        # Loose turn check on non-base vertices (cheap proxy for "not
        # totally folded up"). The geometric dent check we have works
        # only at alpha = pi/3 because reconstruction assumes Euclidean
        # unit triangles; at intermediate alpha that's a fictitious
        # shape. So during homotopy we use the proxy; the caller can
        # dent_index() the final coords separately when alpha = pi/3.
        # Pancake polyhedra correctly converge to bends ∈ {0, ±π} with
        # turn=0 at flat interior vertices — keep this loose so we
        # don't reject them.
        turns_ok = True
        for v in tri.vertices:
            if v in base_face:
                continue
            t_v = vertex_turn(tri, v, trial_bend)
            if t_v < turn_eps:
                turns_ok = False
                break

        if not (out.success and turns_ok):
            step *= 0.5
            successes_in_a_row = 0
            steps.append({"t": t, "alpha_attempted": alpha1,
                          "step_halved_to": step,
                          "newton_success": out.success,
                          "turns_ok": turns_ok,
                          "newton_message": out.message,
                          "residual": float(np.linalg.norm(out.residual))})
            if step < min_step_factor:
                return HomotopyOut(bend, float(np.linalg.norm(out.residual)),
                                   False, steps, iters_total,
                                   f"step shrunk below min_step_factor at t={t:.6f}")
            continue

        bend = trial_bend
        t = t1
        successes_in_a_row += 1
        steps.append({"t": t, "alpha": alpha1,
                      "residual": float(np.linalg.norm(out.residual)),
                      "newton_iters": out.iters,
                      "step": step})
        # Grow step on a streak of clean successes
        if successes_in_a_row >= grow_after and step < max_step_factor:
            step = min(max_step_factor, step * grow_factor)
            successes_in_a_row = 0

    # Now final-solve at target_alpha (already there if t reached 1.0)
    return _final_solve(tri, base_face, target_alpha, bend, var_edges, base_edges, steps, iters_total)


def _final_solve(
    tri: Tri, base_face: Face, target_alpha: float,
    bend: Dict[Edge, float], var_edges: List[Edge], base_edges: List[Edge],
    steps: Optional[List[dict]] = None, iters_total: int = 0,
) -> HomotopyOut:
    if steps is None:
        steps = []
    # Square solve at target alpha
    def F(x: np.ndarray) -> np.ndarray:
        return holonomy_residual(tri, base_face, var_edges, x,
                                 {e: bend[e] for e in base_edges}, target_alpha)
    out = newton(F, np.array([bend[e] for e in var_edges]))
    iters_total += out.iters
    final_bend = dict(bend)
    for e, v in zip(var_edges, out.x):
        final_bend[e] = float(v)

    # Solve base-edge bends from base-vertex holonomy residuals
    final_bend = _complete_base_bends(tri, base_face, target_alpha, final_bend, base_edges)

    final_res = float(np.linalg.norm(out.residual))
    steps.append({"t": 1.0, "alpha": target_alpha, "residual": final_res,
                  "newton_iters": out.iters, "phase": "final_solve"})
    return HomotopyOut(final_bend, final_res, out.success, steps,
                       iters_total, "ok" if out.success else f"final newton: {out.message}")


def _complete_base_bends(
    tri: Tri, base_face: Face, alpha: float,
    bend: Dict[Edge, float], base_edges: List[Edge],
) -> Dict[Edge, float]:
    """Closed-form imputation of the 3 base-edge bends via PYP factorization.

    For each base vertex v_i, the puffup product is
        M_v = Y(α) P(b_first) · K · Y(α) P(b_last) = I
    where K is the product of v's d−2 non-base movemats. Then
        X = K · Y(α),   X^T = P(b_last) · Y(α) · P(b_first)
    so factoring X^T as P(A) Y(α) P(B) gives A = b_last (the bend on the
    edge from base_face[i] to base_face[(i+1)%3]). Each base bend is
    recovered exactly once at the appropriate base vertex.
    """
    Yalpha = matz(alpha)
    new_bend = dict(bend)
    recovered_A = [None] * 3
    for bi, v in enumerate(base_face):
        flower = tri.vertex_flower[v]
        k = len(flower)
        if k < 3:
            raise ValueError(f"vertex {v} flower length {k} < 3")
        # K = product of non-base movemats (skip first and last factors,
        # which involve the two base bends incident to v)
        K = np.eye(3)
        for t in range(1, k - 1):
            f = face_with_v_first(tri.faces[flower[t]], v)
            e = edge_key(v, f[2])
            K = K @ movemat(alpha, bend[e])
        X = K @ Yalpha
        # M = X^T factors as P(A) Y(α) P(B); recover A.
        # X^T[2,0] = X[0,2], X^T[1,0] = X[0,1]
        recovered_A[bi] = math.atan2(-X[0, 2], X[0, 1])
    for bi in range(3):
        e = edge_key(base_face[bi], base_face[(bi + 1) % 3])
        new_bend[e] = recovered_A[bi]
    return new_bend


# ── reconstruction (vertex coordinates) ──────────────────────────────────

def reconstruct(
    tri: Tri, base_face: Face, bend: Dict[Edge, float],
) -> Tuple[Dict[Vertex, np.ndarray], float]:
    """Place base_face's three vertices at the canonical gauge
        base_face[0] → (0, 0,  1/2)
        base_face[1] → (0, 0, -1/2)
        base_face[2] → (√3/2, 0, 0)
    (CCW base normal points to −ŷ; polyhedron interior in +ŷ.)

    BFS the dual graph. For each new face whose two shared-edge endpoints
    a, b are already placed and whose previous neighbor face has third
    vertex p, place the new third vertex c by:
        m = (a+b)/2,  ê = b−a (unit),  û = (p−m)/|p−m|,  v̂ = û × ê
        c = m + (√3/2)·(−û·cos θ + v̂·sin θ)
    where θ is the bend on edge (a,b). v̂ points into +interior side (RH
    rule with ê = b−a), so positive θ folds toward interior.

    Returns (vertex positions, 0.0).  Each vertex is placed exactly once,
    so spread is meaningful only as a constant 0; left in for API compat.
    """
    h = math.sqrt(3.0) / 2.0
    coords: Dict[Vertex, np.ndarray] = {}
    coords[base_face[0]] = np.array([0.0, 0.0,  0.5])
    coords[base_face[1]] = np.array([0.0, 0.0, -0.5])
    coords[base_face[2]] = np.array([h,   0.0,  0.0])

    base_idx = next(fi for fi, f in enumerate(tri.faces) if set(f) == set(base_face))
    placed = {base_idx}
    queue = deque([base_idx])
    while queue:
        fi = queue.popleft()
        f = tri.faces[fi]
        for i in range(3):
            a, b = f[i], f[(i + 1) % 3]
            other = tri.directed_edge_face.get((b, a))
            if other is None or other in placed:
                continue
            other_f = tri.faces[other]
            c = next(v for v in other_f if v != a and v != b)
            p = next(v for v in f       if v != a and v != b)
            theta = bend[edge_key(a, b)]
            pa, pb, pp = coords[a], coords[b], coords[p]
            m = 0.5 * (pa + pb)
            e_hat = pb - pa
            p_perp = pp - m
            u_hat = p_perp / np.linalg.norm(p_perp)
            v_hat = np.cross(u_hat, e_hat)
            v_hat /= np.linalg.norm(v_hat)
            c_perp = -u_hat * math.cos(theta) + v_hat * math.sin(theta)
            coords[c] = m + h * c_perp
            placed.add(other)
            queue.append(other)
    return coords, 0.0


def _transition(
    pf: Dict[Vertex, np.ndarray], pg: Dict[Vertex, np.ndarray],
    e_uw: Tuple[Vertex, Vertex], theta: float,
) -> Tuple[np.ndarray, np.ndarray]:
    u, w = e_uw
    af = pf[w] - pf[u]
    ag = pg[w] - pg[u]
    # Rotate g's frame so its (u->w) aligns with f's (u->w)
    phi = math.atan2(af[1], af[0]) - math.atan2(ag[1], ag[0])
    cp, sp = math.cos(phi), math.sin(phi)
    A0 = np.array([[cp, -sp, 0.0], [sp, cp, 0.0], [0.0, 0.0, 1.0]])
    t0 = pf[u] - A0 @ pg[u]
    # Now fold by -theta about the shared edge axis
    axis = af / np.linalg.norm(af)
    R_fold = _rotation_about_axis(axis, -theta)
    A = R_fold @ A0
    t = R_fold @ (t0 - pf[u]) + pf[u]
    return A, t


def _rotation_about_axis(n: np.ndarray, theta: float) -> np.ndarray:
    nx, ny, nz = n
    c, s = math.cos(theta), math.sin(theta)
    C = 1.0 - c
    return np.array([
        [c + nx * nx * C, nx * ny * C - nz * s, nx * nz * C + ny * s],
        [ny * nx * C + nz * s, c + ny * ny * C, ny * nz * C - nx * s],
        [nz * nx * C - ny * s, nz * ny * C + nx * s, c + nz * nz * C],
    ])


# ── input / output glue ──────────────────────────────────────────────────

def normalize_initial_bends(spec, edges: List[Edge]) -> Dict[Edge, float]:
    if isinstance(spec, (int, float)):
        return {e: float(spec) for e in edges}
    if not isinstance(spec, dict):
        raise ValueError("initial_bends must be a number or {edge: angle} dict")
    out: Dict[Edge, float] = {}
    for e in edges:
        key = edge_str(e)
        if key in spec:
            out[e] = float(spec[key])
        elif f"{e[1]}-{e[0]}" in spec:
            out[e] = float(spec[f"{e[1]}-{e[0]}"])
        else:
            raise ValueError(f"initial_bends missing edge {key}")
    return out


def run_one(payload: dict) -> dict:
    netcode = payload["netcode"]
    faces = parse_netcode(netcode)
    tri = Tri.from_faces(faces)
    base_face = tuple(payload.get("base_face") or list(faces[0]))
    if base_face not in faces and tuple(reversed(base_face)) not in faces:
        # try cyclic permutations
        for f in faces:
            if set(f) == set(base_face):
                base_face = f
                break
    if set(base_face) not in [set(f) for f in faces]:
        raise ValueError(f"base_face {base_face} is not a face of the netcode")

    alpha = float(payload.get("corner_angle", math.pi / 3))
    init = normalize_initial_bends(payload["initial_bends"], tri.edges)

    homo = payload.get("homotopy")
    if homo is None:
        out = solve_homotopy(tri, base_face, alpha, init, from_alpha=alpha)
    else:
        out = solve_homotopy(
            tri, base_face, alpha, init,
            from_alpha=float(homo.get("from_corner", 0.0)),
            step_factor=float(homo.get("step_factor", 1.0 / 24.0)),
            min_step_factor=float(homo.get("min_step_factor", 1e-3)),
        )

    avg_pos, spread = reconstruct(tri, base_face, out.bend_full)
    turns = {str(v): vertex_turn(tri, v, out.bend_full) for v in tri.vertices}

    return {
        "input": payload,
        "solver": {
            "method": "matz_matx_vertex_holonomy_homotopy",
            "success": out.success,
            "homotopy": homo is not None,
            "homotopy_steps": out.steps,
            "iterations_total": out.iters_total,
            "final_residual_norm": out.final_residual,
            "vertex_copy_spread": spread,
            "message": out.message,
        },
        "bends": {edge_str(e): out.bend_full[e] for e in tri.edges},
        "vertex_turns": turns,
        "vertices": {str(v): avg_pos[v].tolist() for v in tri.vertices},
    }


# ── CLI ──────────────────────────────────────────────────────────────────

def main(argv: List[str]) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--in", dest="inp", default="-",
                   help="input JSON file ('-' for stdin)")
    p.add_argument("--out", dest="outp", default="-",
                   help="output JSON file ('-' for stdout)")
    p.add_argument("--indent", type=int, default=2,
                   help="JSON indent (default 2; 0 for compact)")
    args = p.parse_args(argv)

    text = sys.stdin.read() if args.inp == "-" else open(args.inp).read()
    payload = json.loads(text)
    result = run_one(payload)
    out_text = json.dumps(result, indent=args.indent if args.indent > 0 else None)
    if args.outp == "-":
        sys.stdout.write(out_text + "\n")
    else:
        with open(args.outp, "w") as fh:
            fh.write(out_text + "\n")
    return 0 if result["solver"]["success"] else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
