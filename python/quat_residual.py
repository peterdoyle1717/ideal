"""Quaternion holonomy residual + analytic sparse Jacobian.

Mirror of puffup.holonomy_residual + jacobian_sparse.analytical_jacobian_sparse,
using quaternion products instead of 3×3 matrices. Convention identical to
the matrix backend: at each flower step, movemat(α, β) = matz(α)·matx(−β).
The corresponding quaternion factor is q_step(α, β) = q_z(α) · q_x(−β):

    q_z(α)   = (cos(α/2),  0,           0, sin(α/2))
    q_x(−β)  = (cos(β/2), -sin(β/2),    0,        0)

Hamilton product:
    (aw, ax, ay, az) · (bw, bx, by, bz) =
      ( aw·bw − ax·bx − ay·by − az·bz,
        aw·bx + ax·bw + ay·bz − az·by,
        aw·by − ax·bz + ay·bw + az·bx,
        aw·bz + ax·by − ay·bx + az·bw )

So:
    q_step(α, β) = ( ca·cb,
                    -ca·sb,
                    -sa·sb,
                     sa·cb )    where ca=cos(α/2), sa=sin(α/2),
                                       cb=cos(β/2), sb=sin(β/2)

Convention spot-check: building Q via q_step + q_mul and converting
back through R(q) reproduces vertex_holomat to within 2.2e-16 on the
v=4 tet at α=30° with horou ideal bends.

Residual at vertex v: vector part (qx, qy, qz) of the cumulative
holonomy quaternion Q_v walking v's flower. At a true solution
Q_v = ±identity ⇒ vector part = 0.

Branch convention (a once-around-vertex loop accumulates a 2π
rotation in the universal cover, so Q_v naturally lives near
q = (−1, 0, 0, 0), not (+1, 0, 0, 0); measured at α=0 with horou
ideal bends on V=4,6,8 nets, every interior vertex's holonomy is
exactly Q ≈ (−1, 0, 0, 0)):

  - q and −q represent the same rotation; we do NOT silently flip.
  - The residual is just the vector part (qx, qy, qz). LM minimizes it.
  - We do NOT gate per-iteration trials on qw — that's expected to
    drift in transient steps. The residual itself is the only
    convergence criterion.
  - As a post-convergence sanity check, `branch_check_post(...)`
    inspects the final per-vertex (qw, |vec(q)|): if any vertex has
    |vec(q)| small (converged) but qw is *positive* (wrong sheet,
    near +identity instead of −identity), it is a sign that
    something is wrong. The check returns a list of suspect
    vertices; callers decide what to do.
"""
from __future__ import annotations

import math
from typing import Dict, List, Sequence, Tuple

import numpy as np
import scipy.sparse as sp

from puffup import (
    Tri, edge_key, face_with_v_first,
)

Vertex = int
Edge = Tuple[int, int]
Face = Tuple[int, int, int]


def q_mul(a, b):
    aw, ax, ay, az = a
    bw, bx, by, bz = b
    return (
        aw*bw - ax*bx - ay*by - az*bz,
        aw*bx + ax*bw + ay*bz - az*by,
        aw*by - ax*bz + ay*bw + az*bx,
        aw*bz + ax*by - ay*bx + az*bw,
    )


def q_step(alpha: float, beta: float):
    """q_z(α) · q_x(−β). Matches matz(α) · matx(−β)."""
    ca = math.cos(alpha * 0.5); sa = math.sin(alpha * 0.5)
    cb = math.cos(beta  * 0.5); sb = math.sin(beta  * 0.5)
    return (ca*cb, -ca*sb, -sa*sb, sa*cb)


def q_step_dbeta(alpha: float, beta: float):
    """∂/∂β q_step(α, β)."""
    ca = math.cos(alpha * 0.5); sa = math.sin(alpha * 0.5)
    cb = math.cos(beta  * 0.5); sb = math.sin(beta  * 0.5)
    half = 0.5
    return (-half*ca*sb,
            -half*ca*cb,
            -half*sa*cb,
            -half*sa*sb)


def _flower_steps(tri: Tri, v: Vertex, alpha: float, bend: Dict[Edge, float]):
    """Per-flower-step quaternion factors and derivatives for vertex v."""
    qs = []
    dqs = []
    edges: List[Edge] = []
    for fi in tri.vertex_flower[v]:
        f = face_with_v_first(tri.faces[fi], v)
        e = edge_key(v, f[2])
        b = bend[e]
        qs.append(q_step(alpha, b))
        dqs.append(q_step_dbeta(alpha, b))
        edges.append(e)
    return len(qs), qs, dqs, edges


def vertex_holonomy_quat(tri: Tri, v: Vertex, alpha: float,
                          bend: Dict[Edge, float]):
    """Cumulative holonomy quaternion at vertex v. Pure compute, no
    branch check — call branch_check_post on a converged solution."""
    Q = (1.0, 0.0, 0.0, 0.0)
    for fi in tri.vertex_flower[v]:
        f = face_with_v_first(tri.faces[fi], v)
        b = bend[edge_key(v, f[2])]
        Q = q_mul(Q, q_step(alpha, b))
    return Q


def holonomy_residual_quat(
    tri: Tri,
    base_face: Face,
    var_edges: Sequence[Edge],
    var_values: np.ndarray,
    fixed_bend: Dict[Edge, float],
    alpha: float,
    vertices: "Sequence[Vertex] | None" = None,
) -> np.ndarray:
    """Stacked vector parts (qx, qy, qz) of per-vertex holonomy quats.

    Default (vertices=None): interior vertices only (V−3 of them).
    Pass vertices=tri.vertices for the all-bends overdetermined system
    where every vertex's holonomy contributes 3 residual rows."""
    bend = dict(fixed_bend)
    for e, val in zip(var_edges, var_values):
        bend[e] = float(val)
    if vertices is None:
        vertices = [v for v in tri.vertices if v not in base_face]
    out = np.empty(3 * len(vertices), dtype=float)
    for i, v in enumerate(vertices):
        Q = vertex_holonomy_quat(tri, v, alpha, bend)
        out[3*i + 0] = Q[1]
        out[3*i + 1] = Q[2]
        out[3*i + 2] = Q[3]
    return out


def analytical_jacobian_quat_sparse(
    tri: Tri,
    base_face: Face,
    var_edges: Sequence[Edge],
    var_values: np.ndarray,
    fixed_bend: Dict[Edge, float],
    alpha: float,
    vertices: "Sequence[Vertex] | None" = None,
) -> sp.csc_matrix:
    """Sparse CSC analytic Jacobian of holonomy_residual_quat.

    Default (vertices=None): rows correspond to interior-vertex
    holonomies (square system).
    All-bends mode: pass vertices=tri.vertices so rows correspond to
    every vertex (3V × E with E = len(var_edges) = 3V − 6 typically).

    Prefix/suffix product trick:
      Q = q_0 · … · q_{k-1};   P[t] = ∏_{s<t} q_s,   S[t] = ∏_{s≥t} q_s
      ∂Q/∂β_t = P[t] · (dq_t/dβ_t) · S[t+1]
    The vector part of ∂Q/∂β_t gives the 3 Jacobian entries for that
    (vertex, var-edge) at rows 3·vi + (0,1,2)."""
    bend = dict(fixed_bend)
    for e, val in zip(var_edges, var_values):
        bend[e] = float(val)
    if vertices is None:
        vertices = [v for v in tri.vertices if v not in base_face]
    var_idx = {e: i for i, e in enumerate(var_edges)}
    n_vars = len(var_edges)
    n_rows = 3 * len(vertices)

    rows: List[int] = []
    cols: List[int] = []
    vals: List[float] = []

    for vi, v in enumerate(vertices):
        k, qs, dqs, edges = _flower_steps(tri, v, alpha, bend)
        # Prefix products P[0..k]: P[t] = q_0 · … · q_{t-1}, P[0] = identity.
        P = [(1.0, 0.0, 0.0, 0.0)]
        for t in range(k):
            P.append(q_mul(P[t], qs[t]))
        # Suffix products S[0..k]: S[t] = q_t · … · q_{k-1}, S[k] = identity.
        S = [(1.0, 0.0, 0.0, 0.0)] * (k + 1)
        for t in range(k - 1, -1, -1):
            S[t] = q_mul(qs[t], S[t + 1])

        row_x = 3 * vi + 0
        row_y = 3 * vi + 1
        row_z = 3 * vi + 2
        for t in range(k):
            e = edges[t]
            col = var_idx.get(e)
            if col is None:
                continue
            dq = dqs[t]
            tmp = q_mul(P[t], dq)
            dQ = q_mul(tmp, S[t + 1])
            rows.append(row_x); cols.append(col); vals.append(dQ[1])
            rows.append(row_y); cols.append(col); vals.append(dQ[2])
            rows.append(row_z); cols.append(col); vals.append(dQ[3])

    if not rows:
        return sp.csc_matrix((n_rows, n_vars))
    J = sp.coo_matrix((vals, (rows, cols)),
                      shape=(n_rows, n_vars)).tocsc()
    J.sum_duplicates()
    return J


def branch_check_post(
    tri: Tri,
    base_face: Face,
    var_edges: Sequence[Edge],
    var_values: np.ndarray,
    fixed_bend: Dict[Edge, float],
    alpha: float,
    vec_tol: float = 1e-6,
    vertices: "Sequence[Vertex] | None" = None,
):
    """Post-convergence sanity. For each chosen vertex, compute Q
    and report (vec_norm, qw). Default = interior vertices (square
    system). All-bends mode: pass vertices=tri.vertices.

    Flag any vertex where the residual has converged (|vec(Q)| <
    vec_tol) but qw is positive — that is the wrong sheet relative
    to the natural −1 lift for a once-around loop."""
    bend = dict(fixed_bend)
    for e, val in zip(var_edges, var_values):
        bend[e] = float(val)
    if vertices is None:
        vertices = [v for v in tri.vertices if v not in base_face]
    out = []
    for v in vertices:
        Q = vertex_holonomy_quat(tri, v, alpha, bend)
        vec_norm = math.sqrt(Q[1]**2 + Q[2]**2 + Q[3]**2)
        suspect = (vec_norm < vec_tol) and (Q[0] > 0.0)
        out.append({
            "vertex":   v,
            "qw":       Q[0],
            "qx":       Q[1],
            "qy":       Q[2],
            "qz":       Q[3],
            "vec_norm": vec_norm,
            "suspect":  suspect,
        })
    return out


def fd_jacobian_quat(tri, base_face, var_edges, var_values, fixed_bend,
                      alpha, h: float = 1e-7) -> np.ndarray:
    """Reference FD Jacobian for the analytic-vs-FD spot check only.
    NOT used on the solver path. Dense ndarray, central differences."""
    var_values = np.asarray(var_values, dtype=float)
    n = len(var_values)
    f0 = holonomy_residual_quat(tri, base_face, var_edges, var_values,
                                 fixed_bend, alpha)
    J = np.zeros((len(f0), n))
    for j in range(n):
        x_plus = var_values.copy(); x_plus[j] += h
        x_minus = var_values.copy(); x_minus[j] -= h
        f_plus  = holonomy_residual_quat(tri, base_face, var_edges, x_plus,
                                           fixed_bend, alpha)
        f_minus = holonomy_residual_quat(tri, base_face, var_edges, x_minus,
                                           fixed_bend, alpha)
        J[:, j] = (f_plus - f_minus) / (2 * h)
    return J
