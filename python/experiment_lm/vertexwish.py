"""Vertex-wish start rule for the LM α-march experiment.

Solve, in revolution units (1 rev = 2π rad):

  minimize   sum_e  [(x_e - 1/d_i)^2 + (x_e - 1/d_j)^2]    e = {i, j}
  subject to sum_{e incident to v}  x_e = 1   for every vertex v

The unconstrained minimizer per edge is invavg_e = (1/d_i + 1/d_j)/2
(the existing `invavg` start rule). This QP projects invavg onto the
per-vertex sum-to-1 constraint via a Lagrange-multiplier (KKT) solve.

In standard form (minimize (1/2) x^T H x - g^T x s.t. B x = c):

  H        = 2 I_E
  g_e      = 1/d_i + 1/d_j
  c        = 1_V
  B (V x E): B_{i,e} = 1 if e is incident to i

KKT system:

  [ H   B^T ] [ x ]   [ g ]
  [ B    0  ] [ λ ] = [ c ]

Eliminating x = (g - B^T λ) / 2 gives the V x V system

  (B B^T) λ = B g - 2 c

and then x = (g - B^T λ) / 2.

(B B^T) is the signless graph Laplacian D + A. It is positive definite
for any connected non-bipartite graph; a closed triangulated sphere
contains triangles (odd cycles) and is connected, so the system is
non-singular (Codex audit, 2026-05-07).

The radians version multiplies by 2π.
"""
from __future__ import annotations

from typing import Dict, List, Tuple

import math
import numpy as np


def vertexwish_start_revolutions(
    tri,
    deg: Dict[int, int],
) -> Dict[Tuple[int, int], float]:
    edges = list(tri.edges)
    verts = list(tri.vertices)
    E = len(edges)
    V = len(verts)
    vert_idx = {v: k for k, v in enumerate(verts)}

    # Incidence B (V x E).
    B = np.zeros((V, E))
    for k, (i, j) in enumerate(edges):
        B[vert_idx[i], k] = 1.0
        B[vert_idx[j], k] = 1.0

    # g_e = 1/d_i + 1/d_j; c = 1_V.
    g = np.zeros(E)
    for k, (i, j) in enumerate(edges):
        g[k] = 1.0 / deg[i] + 1.0 / deg[j]
    c = np.ones(V)

    # (B B^T) λ = B g - 2 c
    BBt = B @ B.T
    rhs = B @ g - 2.0 * c
    try:
        lam = np.linalg.solve(BBt, rhs)
    except np.linalg.LinAlgError:
        lam, *_ = np.linalg.lstsq(BBt, rhs, rcond=None)

    x = (g - B.T @ lam) / 2.0

    # Sanity invariant: B x ≈ 1_V (within numerical tolerance).
    res = float(np.max(np.abs(B @ x - c)))
    if res > 1e-9:
        raise RuntimeError(
            f"vertexwish: per-vertex sum-to-1 violated by {res:.3e}"
        )

    return {edges[k]: float(x[k]) for k in range(E)}


def vertexwish_start_radians(tri, deg) -> Dict[Tuple[int, int], float]:
    rev = vertexwish_start_revolutions(tri, deg)
    tau = 2.0 * math.pi
    return {e: tau * v for e, v in rev.items()}
