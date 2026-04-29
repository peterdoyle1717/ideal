#!/usr/bin/env python3
"""Analytical (closed-form) Jacobian of holonomy_residual.

Each non-base vertex v contributes a 3-row block. Walking v's flower of
movemats M_1 ... M_k, the partial of the cumulative product w.r.t. the
bend on flower-edge i is

    ∂(M_1 ... M_k)/∂β_i = P[i] · (∂M_i/∂β_i) · S[i+1]
    ∂M_i/∂β_i           = matz(α) · (−J_x) · matx(−β_i)

where P / S are prefix / suffix products and J_x = [[0,0,0],[0,0,-1],[0,1,0]]
(generator of x-rotations: dR_x(θ)/dθ = J_x · R_x(θ), so dR_x(−β)/dβ =
−J_x · R_x(−β)).

Sparsity: each column has nonzeros only at the (≤2) non-base vertices
that the edge touches. We don't materialize a sparse matrix here — we
build a dense n×n J that is structurally sparse — but a sparse build is
trivial if we want it later."""

import math
import numpy as np
from typing import Dict, List, Sequence, Tuple

# import shared types/utilities from puffup
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from puffup import (Tri, edge_key, face_with_v_first, matz, matx, movemat)

Vertex = int
Edge = Tuple[Vertex, Vertex]

_JX = np.array([[0., 0., 0.],
                [0., 0., -1.],
                [0., 1., 0.]])


def analytical_jacobian(tri: Tri, base_face, var_edges: Sequence[Edge],
                         var_values: np.ndarray, fixed_bend: Dict[Edge, float],
                         alpha: float) -> np.ndarray:
    """Closed-form Jacobian of holonomy_residual."""
    bend = dict(fixed_bend)
    for e, val in zip(var_edges, var_values):
        bend[e] = float(val)
    interior_vs = [v for v in tri.vertices if v not in base_face]
    var_idx = {e: i for i, e in enumerate(var_edges)}
    n_vars = len(var_edges)
    n_rows = 3 * len(interior_vs)
    J = np.zeros((n_rows, n_vars))

    Mz = matz(alpha)
    neg_Mz_Jx = -Mz @ _JX  # constant within this α

    for vi, v in enumerate(interior_vs):
        flower = tri.vertex_flower[v]
        Ms: List[np.ndarray] = []
        edges: List[Edge] = []
        for fi in flower:
            f = face_with_v_first(tri.faces[fi], v)
            c = f[2]
            e = edge_key(v, c)
            edges.append(e)
            Ms.append(movemat(alpha, bend[e]))
        k = len(Ms)

        # prefix products: P[i] = M_0 @ ... @ M_{i-1};  P[0] = I
        P = [np.eye(3)]
        for M in Ms:
            P.append(P[-1] @ M)
        # suffix products: S[i] = M_i @ ... @ M_{k-1};  S[k] = I
        S = [None] * (k + 1)
        S[k] = np.eye(3)
        for i in range(k - 1, -1, -1):
            S[i] = Ms[i] @ S[i + 1]

        row = 3 * vi
        for i, e in enumerate(edges):
            if e not in var_idx:
                continue
            # dM_i / dβ_i = -Mz · Jx · matx(-β_i)
            dM = neg_Mz_Jx @ matx(-bend[e])
            contrib = P[i] @ dM @ S[i + 1]
            col = var_idx[e]
            J[row + 0, col] += contrib[0, 1]
            J[row + 1, col] += contrib[0, 2]
            J[row + 2, col] += contrib[1, 2]
    return J
