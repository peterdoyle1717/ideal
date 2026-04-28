#!/usr/bin/env python3
"""
Dihedral (bend) angle on each edge of a horoball-packed triangulation.

Given a horou solution u (vertex -> weight), compute the bend angle
along each edge. Convention: bend = π − interior_dihedral, so flat = 0
and a convex fold gives bend > 0. Bends sum to 2π around every vertex
(holonomy identity in the α=0 ideal limit).

The "vertex at infinity" v∞ is the vertex with u[v∞] = ∞; its neighbors
have u = 1 (these are the boundary of the finite triangulation).

Three cases for an edge {i, j}:
  1. v∞ ∈ {i, j}        (boundary edge)        → π − Σ finite-face petals at the other endpoint
  2. both adjacent faces finite                → π − petal_a − petal_d
  3. one adjacent face finite, the other not   → π − petal at the third vertex of the finite face
"""

import math
from typing import Dict, List, Sequence, Tuple

Vertex = int
Face = Tuple[Vertex, Vertex, Vertex]
Edge = Tuple[Vertex, Vertex]


def _inf_vertex(u: Dict[Vertex, float]) -> Vertex:
    for v, x in u.items():
        if math.isinf(x):
            return v
    raise ValueError("u dict has no infinite vertex")


def _other(face: Face, edge: Edge) -> Vertex:
    a, b = edge
    return next(v for v in face if v != a and v != b)


def _petal(v: Vertex, face: Face, u: Dict[Vertex, float]) -> float:
    """Angle at v in face, in the horou metric (edge length = u[i]·u[j])."""
    p, q = (w for w in face if w != v)
    x, y, z = 1.0/u[v], 1.0/u[p], 1.0/u[q]
    cos_t = (y*y + z*z - x*x) / (2.0 * y * z)
    return math.acos(max(-1.0, min(1.0, cos_t)))


def all_edges(faces: Sequence[Face]) -> List[Edge]:
    es = set()
    for f in faces:
        for i in range(3):
            a, b = f[i], f[(i+1) % 3]
            es.add((a, b) if a < b else (b, a))
    return sorted(es)


def dihedral(edge: Edge, u: Dict[Vertex, float],
             faces: Sequence[Face], v_inf: Vertex = None) -> float:
    if v_inf is None:
        v_inf = _inf_vertex(u)

    a, b = edge
    fs = [f for f in faces if a in f and b in f]
    if len(fs) != 2:
        raise ValueError(f"edge {edge} adjoins {len(fs)} faces, need 2")

    if v_inf in edge:
        v = a if b == v_inf else b
        return math.pi - sum(_petal(v, f, u) for f in faces
                             if v in f and v_inf not in f)

    finite = [f for f in fs if v_inf not in f]
    if len(finite) == 2:
        return (math.pi
                - _petal(_other(finite[0], edge), finite[0], u)
                - _petal(_other(finite[1], edge), finite[1], u))
    # exactly one finite, one through v_inf
    f0 = finite[0]
    return math.pi - _petal(_other(f0, edge), f0, u)


def all_dihedrals(faces: Sequence[Face],
                  u: Dict[Vertex, float]) -> Dict[Edge, float]:
    v_inf = _inf_vertex(u)
    return {e: dihedral(e, u, faces, v_inf) for e in all_edges(faces)}
