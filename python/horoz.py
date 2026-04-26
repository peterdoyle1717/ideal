#!/usr/bin/env python3
"""
Ideal horoball positions: BFS placement of vertices in the upper half-plane.

Given a triangulated sphere with horoball weights u[v] (from horou),
compute flat positions (x[v], y[v]) for each vertex.

The edge length between adjacent vertices i, j is u[i]*u[j].
Vertex v0 (point at infinity) has no finite position.
Vertices v1, v2 (first two boundary neighbors) are pinned at (0,0) and (1,0).
Remaining vertices are placed by BFS over directed edges using the thirdpoint formula.

Interface
---------
    pos = horoz(poly, u)
    # Returns dict: vertex -> (x, y).  v0 is absent (point at infinity).
    # v1=(0,0), v2=(1,0) by convention.
"""

import math
from collections import deque


def _build_em(poly):
    """Build directed-edge map: (a,b) -> c for each oriented face (a,b,c)."""
    em = {}
    for a, b, c in poly:
        em[(a, b)] = c
        em[(b, c)] = a
        em[(c, a)] = b
    return em


def _thirdpoint(xa, ya, xb, yb, dA, dB):
    """
    Place vertex c given directed edge a→b (face (a,b,c) CCW).

    dA = u[a]*u[c]  (edge length a-c in Euclidean horoball metric)
    dB = u[b]*u[c]  (edge length b-c)

    Returns (xc, yc).
    """
    dx, dy = xb - xa, yb - ya
    L = math.sqrt(dx*dx + dy*dy)
    d0, d1 = dA / L, dB / L
    xn = (1.0 + d0*d0 - d1*d1) / 2.0
    yn2 = d0*d0 - xn*xn
    yn = math.sqrt(max(0.0, yn2))
    return xa + dx*xn + dy*yn, ya + dy*xn - dx*yn


def horoz(poly, u):
    """
    Compute flat vertex positions from horoball weights.

    Parameters
    ----------
    poly : list of (a, b, c) oriented triangles (integer 1-indexed vertices).
           poly[0][0] is taken as v0 (point at infinity).
    u    : dict vertex -> weight, as returned by horou().
           v0 maps to float('inf'); boundary vertices map to 1.0.

    Returns
    -------
    dict: vertex -> (x, y).  v0 is absent.
    """
    v0 = poly[0][0]
    em  = _build_em(poly)

    # Identify first two finite vertices (neighbors of v0 in the first face).
    # By convention: v1 = poly[0][1], v2 = poly[0][2].
    v1, v2 = poly[0][1], poly[0][2]

    pos = {v1: (0.0, 0.0), v2: (1.0, 0.0)}
    placed = {v0, v1, v2}

    # BFS over directed edges.  Seed: edge v2→v1 to enter the first interior face.
    q = deque()
    q.append((v2, v1))

    while q:
        a, b = q.popleft()
        c = em.get((a, b))
        if c is None or c in placed:
            continue
        xa, ya = pos[a]
        xb, yb = pos[b]
        dA = u[a] * u[c]
        dB = u[b] * u[c]
        pos[c] = _thirdpoint(xa, ya, xb, yb, dA, dB)
        placed.add(c)
        q.append((c, b))
        q.append((a, c))

    # Cleanup: place any vertex not reached by BFS (rare topology edge cases).
    changed = True
    while changed:
        changed = False
        for a, b, c in poly:
            if a == v0 or b == v0 or c == v0:
                continue
            if a in placed and b in placed and c not in placed:
                dA = u[a] * u[c]; dB = u[b] * u[c]
                pos[c] = _thirdpoint(pos[a][0], pos[a][1], pos[b][0], pos[b][1], dA, dB)
                placed.add(c); changed = True
            elif b in placed and c in placed and a not in placed:
                dB = u[b] * u[a]; dC = u[c] * u[a]
                pos[a] = _thirdpoint(pos[b][0], pos[b][1], pos[c][0], pos[c][1], dB, dC)
                placed.add(a); changed = True
            elif a in placed and c in placed and b not in placed:
                dC = u[c] * u[b]; dA = u[a] * u[b]
                pos[b] = _thirdpoint(pos[c][0], pos[c][1], pos[a][0], pos[a][1], dC, dA)
                placed.add(b); changed = True

    return pos


# ── CLI ───────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    import sys, os
    sys.path.insert(0, os.path.dirname(__file__))
    from horou import horou

    for line in sys.stdin:
        s = line.strip()
        if not s:
            continue
        poly = [tuple(int(x) for x in f.split(',')) for f in s.split(';')]
        u = horou(poly)
        pos = horoz(poly, u)
        v0 = poly[0][0]
        parts = ' '.join(f'{v}:({pos[v][0]:.6g},{pos[v][1]:.6g})' for v in sorted(pos))
        print(parts)
