#!/usr/bin/env python3
"""Sparse build of the analytical holonomy Jacobian.

Same math as jacobian.analytical_jacobian; the only difference is that
the result is assembled in COO format and returned as a scipy CSC matrix.
The Jacobian is structurally sparse: each column (one bend variable on
edge e) has nonzeros only on the (≤2) interior-vertex rows that the edge
touches in its flower walk. NNZ is O(6V), out of (3V)² dense entries.

For V ≈ 2000 this is ~12K nonzeros vs 36M dense — and scipy's sparse
LU (UMFPACK / SuperLU) factorizes in tens of ms vs seconds for LAPACK
dgesv on a dense ~6000×6000.
"""

from __future__ import annotations

import math
import os
import sys
from typing import Dict, List, Sequence, Tuple

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from puffup import Tri, edge_key, face_with_v_first, matz, matx, movemat

Vertex = int
Edge = Tuple[Vertex, Vertex]

_JX = np.array([[0., 0., 0.],
                [0., 0., -1.],
                [0., 1., 0.]])


def analytical_jacobian_sparse(tri: Tri, base_face, var_edges: Sequence[Edge],
                               var_values: np.ndarray, fixed_bend: Dict[Edge, float],
                               alpha: float) -> sp.csc_matrix:
    bend = dict(fixed_bend)
    for e, val in zip(var_edges, var_values):
        bend[e] = float(val)
    interior_vs = [v for v in tri.vertices if v not in base_face]
    var_idx = {e: i for i, e in enumerate(var_edges)}
    n_vars = len(var_edges)
    n_rows = 3 * len(interior_vs)

    Mz = matz(alpha)
    neg_Mz_Jx = -Mz @ _JX

    # Accumulate (row, col, val) triples. For each interior vertex's flower
    # we get k edges × 3 rows = ~18 entries per vertex; total NNZ ≈ 18·N_INT.
    rows: List[int] = []
    cols: List[int] = []
    vals: List[float] = []

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

        P = [np.eye(3)]
        for M in Ms:
            P.append(P[-1] @ M)
        S: List[np.ndarray] = [None] * (k + 1)  # type: ignore[list-item]
        S[k] = np.eye(3)
        for i in range(k - 1, -1, -1):
            S[i] = Ms[i] @ S[i + 1]

        row0 = 3 * vi
        for i, e in enumerate(edges):
            if e not in var_idx:
                continue
            dM = neg_Mz_Jx @ matx(-bend[e])
            contrib = P[i] @ dM @ S[i + 1]
            col = var_idx[e]
            for r_off, mij in ((0, contrib[0, 1]), (1, contrib[0, 2]), (2, contrib[1, 2])):
                if mij != 0.0:
                    rows.append(row0 + r_off); cols.append(col); vals.append(mij)

    if not vals:
        return sp.csc_matrix((n_rows, n_vars))
    return sp.coo_matrix((vals, (rows, cols)), shape=(n_rows, n_vars)).tocsc()
