#!/usr/bin/env python3
"""gen_thick_pancakes.py — CLERS of thick pancakes with hip-roof side caps.

Take a planar lattice polygon (same region as gen_pancakes.py). Cross it
with a unit interval to get a prism: top cap + bottom cap + rectangular
side faces, one per polygon boundary edge. Replace each rectangle with a
hip roof (ridge of ridge-vertices parallel to the polygon edge,
triangulated outward). Per unit of polygon boundary, the hip roof adds
one ridge vertex.

Vertex count: V = 2*|polygon pts| + perimeter  (perimeter counted in
lattice units).
Face count:   F = 2*T_poly + 4*perimeter       (T_poly = triangles in
polygon triangulation; 2 triangles per unit-rectangle become 4 for the
pyramidal/ridge cap).
"""
from __future__ import annotations
import sys, os, argparse

from . import _CLERS_SRC  # noqa: F401  (canonical clers on path)
from .gen_pancakes import region_points, classify
from clers import official_unoriented_name


def enumerate_polygons(vmax):
    """Yield (k, l, r, q) for every lattice polygon (3-6 sides) whose
    thick pancake fits in v <= vmax. Includes degenerate (pentagon, quad,
    triangle) cases where an edge length collapses to zero."""
    for k in range(1, vmax + 1):
        for l in range(k, vmax + 1):
            for r in range(0, k + l + 1):
                for q in range(r + 1, k + l + 1):
                    # Quick V upper bound: thick V = 2|pts| + perimeter.
                    # |pts| <= (k+1)(l+1), perimeter <= k+l+k+l = 2(k+l).
                    if 2 * (k + 1) * (l + 1) + 2 * (k + l) > 4 * vmax:
                        continue
                    yield (k, l, r, q)


def _boundary_cycle(pts_set):
    """Walk the polygon boundary CCW. Returns list of (a,b) lattice points
    in order, visiting every point at most once. Uses the hex-lattice
    adjacency convention (6 dirs) with a consistent rotation."""
    # Neighbor directions in CCW order (hex lattice)
    dirs_ccw = [(1, 0), (1, -1), (0, -1), (-1, 0), (-1, 1), (0, 1)]

    def is_boundary(p):
        return any((p[0] + d[0], p[1] + d[1]) not in pts_set for d in dirs_ccw)

    bdry = [p for p in pts_set if is_boundary(p)]
    if not bdry:
        return []
    start = min(bdry)
    cycle = [start]
    cur = start
    prev_dir = 3  # came from "-x"; start outgoing search from that direction
    while True:
        found = False
        for i in range(6):
            d_idx = (prev_dir + i) % 6
            d = dirs_ccw[d_idx]
            nxt = (cur[0] + d[0], cur[1] + d[1])
            if nxt in pts_set and is_boundary(nxt):
                # verify nxt is actually adjacent on the boundary (shared
                # edge): the edge (cur, nxt) must have outside on one side.
                # Cheap heuristic: accept the first CCW neighbor that isn't
                # the one we came from.
                if len(cycle) > 1 and nxt == cycle[-2]:
                    continue
                cycle.append(nxt)
                # next outgoing search starts one step past the reversed direction
                prev_dir = (d_idx + 3 + 1) % 6
                cur = nxt
                found = True
                break
        if not found or cur == start:
            break
    if cycle and cycle[-1] == start:
        cycle.pop()
    return cycle


def build_thick_pancake(k, l, r, q):
    """Return (V, faces) or None if polygon is degenerate. Faces 1-indexed."""
    pts = region_points(k, l, r, q)
    if len(pts) < 3:
        return None
    pts_set = set(pts)
    _, inter = classify(pts)

    # Index: every polygon point has a top copy (T) and a bottom copy (B).
    idx = {}
    n = 0
    for p in sorted(pts):
        n += 1; idx[('T', p)] = n
        n += 1; idx[('B', p)] = n

    # Boundary cycle (CCW)
    cycle = _boundary_cycle(pts_set)
    if len(cycle) < 3:
        return None

    # One apex per unit of polygon boundary. The boundary walk moves
    # one lattice step at a time; one apex is created per step.
    apex_idx = []
    for _ in cycle:
        n += 1
        apex_idx.append(n)

    faces = []
    # Top cap: triangulate polygon (same pattern as pancake top sheet).
    for a in range(-1, k + 1):
        for b in range(-1, l + 1):
            if (a, b) in pts_set and (a + 1, b) in pts_set and (a, b + 1) in pts_set:
                faces.append((idx[('T', (a, b))], idx[('T', (a + 1, b))], idx[('T', (a, b + 1))]))
            if (a + 1, b) in pts_set and (a + 1, b + 1) in pts_set and (a, b + 1) in pts_set:
                faces.append((idx[('T', (a + 1, b))], idx[('T', (a + 1, b + 1))], idx[('T', (a, b + 1))]))
    # Bottom cap: same triangulation, opposite orientation.
    for a in range(-1, k + 1):
        for b in range(-1, l + 1):
            if (a, b) in pts_set and (a + 1, b) in pts_set and (a, b + 1) in pts_set:
                faces.append((idx[('B', (a, b))], idx[('B', (a, b + 1))], idx[('B', (a + 1, b))]))
            if (a + 1, b) in pts_set and (a + 1, b + 1) in pts_set and (a, b + 1) in pts_set:
                faces.append((idx[('B', (a + 1, b))], idx[('B', (a, b + 1))], idx[('B', (a + 1, b + 1))]))

    # Side hip roofs: each unit rectangle of the prism side gets an apex.
    # Adjacent apexes are connected by a ridge edge only if the polygon
    # boundary continues STRAIGHT between them (same direction in/out at
    # the shared mid-edge vertex) — so the ridge runs along the polygon
    # edge. At polygon CORNERS (direction change), the adjacent pyramids
    # are independent (4-face pyramid each).
    n_steps = len(cycle)

    def step_dir(i):
        a = cycle[i]
        b = cycle[(i + 1) % n_steps]
        return (b[0] - a[0], b[1] - a[1])

    # is_corner[i] = True if the polygon turns at cycle[i] (between the
    # incoming edge (i-1)->(i) and outgoing edge (i)->(i+1)).
    is_corner = [step_dir(i - 1) != step_dir(i) for i in range(n_steps)]

    for i in range(n_steps):
        cur = cycle[i]
        nxt = cycle[(i + 1) % n_steps]
        a_cur = apex_idx[i]
        a_nxt = apex_idx[(i + 1) % n_steps]
        Tcur, Bcur = idx[('T', cur)], idx[('B', cur)]
        Tnxt, Bnxt = idx[('T', nxt)], idx[('B', nxt)]
        # Two slope triangles of this rectangle:
        faces.append((Tcur, Tnxt, a_cur))
        faces.append((Bnxt, Bcur, a_cur))
        if is_corner[(i + 1) % n_steps]:
            # Corner at nxt: close the pyramid with the two hip triangles
            # that meet only this apex (no ridge to the next apex).
            faces.append((Tnxt, Bnxt, a_cur))
            # and the Bcur-Tcur side triangle is already "closed" by the
            # previous corner's hip — if the previous vertex (cur) is ALSO
            # a corner, this apex is a full pyramid; otherwise this side
            # is closed by the ridge-triangle from the previous step.
        else:
            # Straight: ridge segment to next apex, sharing mid-vertex nxt.
            faces.append((Tnxt, a_nxt, a_cur))
            faces.append((Bnxt, a_cur, a_nxt))
        if is_corner[i]:
            # Close the start of this pyramid (the Bcur-Tcur side) since
            # there's no incoming ridge at cur.
            faces.append((Bcur, Tcur, a_cur))

    V = n
    return V, faces


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vmax', type=int, default=30)
    args = ap.parse_args()
    seen = set()
    for k, l, r, q in enumerate_polygons(args.vmax):
        res = build_thick_pancake(k, l, r, q)
        if res is None:
            continue
        V, faces = res
        if V > args.vmax:
            continue
        if len(faces) != 2 * V - 4:
            continue
        code = official_unoriented_name(faces)
        if code in seen:
            continue
        seen.add(code)
        print(f"v={V} (k,l,r,q)=({k},{l},{r},{q})  {code}")


if __name__ == '__main__':
    main()
