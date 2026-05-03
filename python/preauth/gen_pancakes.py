#!/usr/bin/env python3
"""gen_pancakes.py — generate CLERS codes of all prime pancakes up to v vertices.

A pancake is the double cover of a 2D lattice region. The region is
parameterized by (k, l, r, q):

    0 <= a <= k,  0 <= b <= l,  r <= a+b <= q

giving a polygon of up to six sides (a=0, a=k, b=0, b=l, a+b=r, a+b=q).

Double-cover construction: boundary vertices get one 3D copy; interior
vertices get two copies (top and bottom). Each 2D triangle becomes two
3D triangles (front + reversed back).

Constraints for a proper hex-footprint pancake (all 120° corners, no
pointy tips, at least one interior vertex), canonical form:

    k <= l,
    k >= 2,
    1 <= r <= k - 1,
    l + 1 <= q <= k + l - 1,
    r + q <= k + l        (reflection canonical)

Direct nested iteration over this envelope.
"""
from __future__ import annotations
import sys, os, argparse

# Pulls canonical clers in via preauth/__init__.py side effect.
from . import _CLERS_SRC  # noqa: F401
from clers import official_unoriented_name


def region_points(k, l, r, q):
    return [(a, b) for a in range(k+1) for b in range(l+1) if r <= a+b <= q]


def classify(pts):
    s = set(pts)
    nbr_dirs = [(1,0), (-1,0), (0,1), (0,-1), (1,-1), (-1,1)]
    bdry, inter = set(), set()
    for p in pts:
        if all((p[0]+d[0], p[1]+d[1]) in s for d in nbr_dirs):
            inter.add(p)
        else:
            bdry.add(p)
    return bdry, inter


def pancake_facelist(k, l, r, q):
    """Return (V, faces) for the double-cover pancake."""
    pts = region_points(k, l, r, q)
    _, inter = classify(pts)
    idx = {}
    n = 0
    for p in sorted(pts):
        n += 1; idx[('T', p)] = n
        if p in inter:
            n += 1; idx[('B', p)] = n
        else:
            idx[('B', p)] = idx[('T', p)]
    faces = []
    pts_set = set(pts)
    for a in range(-1, k+1):
        for b in range(-1, l+1):
            if (a,b) in pts_set and (a+1,b) in pts_set and (a,b+1) in pts_set:
                faces.append((idx[('T',(a,b))], idx[('T',(a+1,b))], idx[('T',(a,b+1))]))
                faces.append((idx[('B',(a,b))], idx[('B',(a,b+1))], idx[('B',(a+1,b))]))
            if (a+1,b) in pts_set and (a+1,b+1) in pts_set and (a,b+1) in pts_set:
                faces.append((idx[('T',(a+1,b))], idx[('T',(a+1,b+1))], idx[('T',(a,b+1))]))
                faces.append((idx[('B',(a+1,b))], idx[('B',(a,b+1))], idx[('B',(a+1,b+1))]))
    V = n
    if len(faces) != 2 * V - 4:
        return None
    return V, faces


def pancake_v(k, l, r, q):
    pts = region_points(k, l, r, q)
    _, inter = classify(pts)
    return len(pts) + len(inter)


def pancake_v_formula(k, l, r, q):
    """v = 2kl + 2 − r² − (k+l−q)². Valid for proper-hex params."""
    s = k + l - q
    return 2 * k * l + 2 - r * r - s * s


def enumerate_pancakes(vmax):
    """Yield (k, l, r, q, V) for every canonical proper-hex pancake with V <= vmax.

    v(k,l,r,q) = 2kl + 2 − r² − (k+l−q)². Max r, max (k+l-q) both = k−1,
    so min v at fixed (k,l) is 2kl + 2 − 2(k−1)². Break when that's > vmax.
    """
    for k in range(2, vmax + 1):
        for l in range(k, vmax + 1):
            if 2 * k * l + 2 - 2 * (k - 1) ** 2 > vmax:
                break
            for r in range(1, k):
                for q in range(l + 1, k + l):
                    if r + q > k + l:
                        continue
                    v = pancake_v_formula(k, l, r, q)
                    if 4 <= v <= vmax:
                        yield (k, l, r, q, v)


def edge_lengths(k, l, r, q):
    """The 6 hex-footprint edge lengths in CCW order."""
    return (k - r, q - k, k + l - q, q - l, l - r, r)


def canonical_sextuple(sx):
    """Canonical (min under rotation + reflection) of a 6-tuple."""
    rots = [sx[i:] + sx[:i] for i in range(6)]
    rots += [tuple(reversed(rr)) for rr in rots]
    return min(tuple(r) for r in rots)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vmax', type=int, default=50)
    args = ap.parse_args()

    seen_sigs = {}  # canonical edge sextuple -> (k, l, r, q)
    for k, l, r, q, v in enumerate_pancakes(args.vmax):
        sig = canonical_sextuple(edge_lengths(k, l, r, q))
        if sig not in seen_sigs:
            seen_sigs[sig] = (k, l, r, q, v)

    codes = []
    for sig, (k, l, r, q, v) in seen_sigs.items():
        result = pancake_facelist(k, l, r, q)
        if result is None:
            continue
        V, faces = result
        codes.append((v, official_unoriented_name(faces)))

    codes.sort()
    for _, code in codes:
        print(code)


if __name__ == '__main__':
    main()
