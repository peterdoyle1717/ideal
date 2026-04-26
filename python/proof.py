#!/usr/bin/env python3
"""
proof.py — sub/supersolution existence proof for ideal horoball packings.

Given a triangulated sphere, proves that an ideal horoball packing exists
by bracketing the solution between a sub-solution (angle sums > 2π) and a
super-solution (angle sums < 2π) and verifying five conditions that together
guarantee a true solution lies strictly between them.

Input:  face lists from stdin, one per line: "a,b,c;d,e,f;..."
        (use clers_decode.py to convert CLERS strings)

Output: one line per net:
          PROVED   slack=...  eps=1/N
          UNPROVED ...
          FAILED   (solver did not converge)

Usage:
    echo "CCAE" | python3 ../clers/python/clers_decode.py | python3 proof.py
    python3 ../clers/python/clers_decode.py < ../clers/fuller/prime/20.txt | python3 proof.py
    python3 proof.py --verbose   (show all 5 check values for each net)
"""

import sys
import os
import math
import argparse

sys.path.insert(0, os.path.dirname(__file__))
from horou import horou

TAU = 2.0 * math.pi


# ── triangulation helpers ─────────────────────────────────────────────────────

def _build_em(poly):
    em = {}
    for a, b, c in poly:
        em[(a, b)] = c
        em[(b, c)] = a
        em[(c, a)] = b
    return em


def _cyclic_nbrs(v, em, any_nbr):
    """Return neighbors of v in cyclic order using the directed-edge map."""
    start = any_nbr
    ring = [start]
    cur = em.get((v, start))
    while cur is not None and cur != start:
        ring.append(cur)
        cur = em.get((v, cur))
    return ring


def _build_topology(poly):
    """Return (v0, boundary, interior, em, rings, finite_faces)."""
    v0 = poly[0][0]
    em = _build_em(poly)

    # Cyclic neighbors of v0 = boundary vertices.
    v0_nbrs_any = next(b for (a, b) in em if a == v0)
    bndry = set(_cyclic_nbrs(v0, em, v0_nbrs_any))

    all_verts = sorted({v for f in poly for v in f} - {v0})
    interior = [v for v in all_verts if v not in bndry]

    # Precompute cyclic rings for interior vertices.
    rings = {}
    for v in interior:
        any_nbr = next(b for (a, b) in em if a == v)
        rings[v] = _cyclic_nbrs(v, em, any_nbr)

    finite_faces = [(a, b, c) for a, b, c in poly if v0 not in (a, b, c)]

    return v0, bndry, interior, em, rings, finite_faces


# ── angle geometry ────────────────────────────────────────────────────────────

def _petal(ui, uj, uk):
    """Angle at vertex i in triangle with edge lengths ui*uj, ui*uk, uj*uk."""
    a, b, c = ui*uj, ui*uk, uj*uk
    cos_t = (a*a + b*b - c*c) / (2.0*a*b)
    return math.acos(max(-1.0, min(1.0, cos_t)))


def _flower(v, u, ring):
    """Sum of petal angles around vertex v."""
    k = len(ring)
    return sum(_petal(u[v], u[ring[j]], u[ring[(j+1)%k]]) for j in range(k))


# ── five proof checks ─────────────────────────────────────────────────────────
#
# umin = sub-solution: horou(-eps), so flower(umin[v]) = 2π+eps > 2π everywhere.
# umax = super-solution: horou(+eps), so flower(umax[v]) = 2π-eps < 2π everywhere.
#
# Together they bracket a true solution u* with umin < u* < umax, proved by:
#   1. monocheck   — umax[v] > umin[v] for all interior vertices
#   2. excesscheck — umin is sub-solution, umax is super-solution (verified numerically)
#   3. tricheck    — triangle inequality holds for all u in [umin, umax]
#   4. convcheck   — Delaunay/butterfly condition for all interior edges
#   5. bndrycheck  — boundary horoball angles < π
#
# Each check returns a slack > 0 if it passes. The overall proof slack is the
# minimum across all checks.

def _check_mono(interior, umin, umax):
    """1. umax[v] > umin[v] for all interior v."""
    return min(umax[v] - umin[v] for v in interior)


def _check_excess(interior, rings, umin, umax):
    """2. umin is sub-solution (flower > 2π) and umax is super-solution (flower < 2π)."""
    slack = math.inf
    for v in interior:
        fn_min = _flower(v, umin, rings[v])
        fn_max = _flower(v, umax, rings[v])
        slack = min(slack, fn_min - TAU, TAU - fn_max)
    return slack


def _check_triangle(finite_faces, umin, umax):
    """3. Triangle inequality for any u in [umin, umax].
    Worst case for edge (a,b,c): smallest center (umin[a]), largest sides (umax[b], umax[c]).
    Slack: umin_a/umax_c + umin_a/umax_b - 1 (all three rotations)."""
    slack = math.inf
    for a, b, c in finite_faces:
        uma, umb, umc = umin[a], umin[b], umin[c]
        xma, xmb, xmc = umax[a], umax[b], umax[c]
        slack = min(slack,
                    uma*(xmb+xmc)/(xmb*xmc) - 1.0,
                    umb*(xma+xmc)/(xma*xmc) - 1.0,
                    umc*(xma+xmb)/(xma*xmb) - 1.0)
    return slack


def _check_convex(em, umin, umax, v0):
    """4. Butterfly/Delaunay condition for every interior edge.
    For edge {v,w} shared by faces (v,w,b) and (w,v,d):
      2*umin[b]²*umin[d]²*(umin[v]²+umin[w]²) / (umax[v]²*umax[w]²*(umax[b]²+umax[d]²)) > 1."""
    slack = math.inf
    seen = set()
    for (v, w) in em:
        if v == v0 or w == v0 or (w, v) in seen:
            continue
        b = em.get((v, w))
        d = em.get((w, v))
        if b is None or d is None or b == v0 or d == v0:
            continue
        seen.add((v, w))
        ua2, uc2 = umin[v]**2, umin[w]**2
        ub2, ud2 = umin[b]**2, umin[d]**2
        Va2, Vc2 = umax[v]**2, umax[w]**2
        Vb2, Vd2 = umax[b]**2, umax[d]**2
        num = 2.0 * ub2 * ud2 * (ua2 + uc2)
        den = Va2 * Vc2 * (Vb2 + Vd2)
        slack = min(slack, num/den - 1.0)
    return slack


def _check_boundary(poly, em, umin, umax, v0, bndry):
    """5. Boundary horoball angle < π.
    For each boundary vertex b, sum petals over its open arc of finite neighbors.
    Worst case: umin[b] for center, umax for neighbors."""
    slack = math.inf
    v0_nbr_any = next(nb for (a, nb) in em if a == v0)
    bndry_ring = _cyclic_nbrs(v0, em, v0_nbr_any)
    for b in bndry_ring:
        b_nbr_any = next(nb for (a, nb) in em if a == b)
        b_ring = _cyclic_nbrs(b, em, b_nbr_any)
        # Arc = all neighbors of b except v0.
        try:
            start = (b_ring.index(v0) + 1) % len(b_ring)
        except ValueError:
            continue
        arc = [b_ring[(start + j) % len(b_ring)] for j in range(len(b_ring) - 1)]
        slide = sum(_petal(umin[b], umax[arc[j]], umax[arc[j+1]])
                    for j in range(len(arc) - 1))
        slack = min(slack, (math.pi - slide) / TAU)
    return slack


# ── prover ────────────────────────────────────────────────────────────────────

EPS_START    = 1/500
EPS_MIN      = 1/500000   # stop halving below this


def prove(poly, eps=None, verbose=False):
    """
    Attempt to prove existence of an ideal horoball packing for the given net.

    Parameters
    ----------
    poly    : list of (a,b,c) oriented triangles.
    eps     : fixed defect value, or None for adaptive halving from 1/500.
    verbose : if True, print each check's slack value.

    Returns
    -------
    (proved, slack, eps_used)
    proved  : True if all 5 checks pass.
    slack   : minimum slack (positive = proved; negative = which check failed).
    eps_used: the eps value at which the proof succeeded (or last tried).
    """
    v0, bndry, interior, em, rings, finite_faces = _build_topology(poly)

    if not interior:
        # Trivial: no interior vertices, angle sums are automatic.
        return True, math.inf, 0.0

    def try_eps(e):
        # umin: defect = -e → target = 2π+e → angle sum > 2π (sub-solution)
        # umax: defect = +e → target = 2π-e → angle sum < 2π (super-solution)
        umin_d = horou(poly, defect=-e)
        umax_d = horou(poly, defect=+e)
        if umin_d is None or umax_d is None:
            return None, None

        checks = [
            ('mono',    _check_mono(interior, umin_d, umax_d)),
            ('excess',  _check_excess(interior, rings, umin_d, umax_d)),
            ('triangle',_check_triangle(finite_faces, umin_d, umax_d)),
            ('convex',  _check_convex(em, umin_d, umax_d, v0)),
            ('boundary',_check_boundary(poly, em, umin_d, umax_d, v0, bndry)),
        ]
        s = min(v for _, v in checks)
        return s, checks

    if eps is not None:
        slack, checks = try_eps(eps)
        if verbose and checks:
            for name, val in checks:
                print(f'    {name}: {val:.6g}')
        if slack is None:
            return False, float('nan'), eps
        return slack > 0, slack, eps

    # Adaptive: start at 1/500, halve until proved or give up.
    e = EPS_START
    while e >= EPS_MIN:
        slack, checks = try_eps(e)
        if slack is None:
            e *= 0.5
            continue
        if slack > 0:
            if verbose and checks:
                for name, val in checks:
                    print(f'    {name}: {val:.6g}')
            return True, slack, e
        e *= 0.5

    return False, slack if slack is not None else float('nan'), e * 2


# ── CLI ───────────────────────────────────────────────────────────────────────

def _eps_label(e):
    """Format eps as 1/N if possible."""
    if e == 0:
        return '0'
    n = round(1.0 / e)
    if abs(1.0/n - e) < 1e-12:
        return f'1/{n}'
    return f'{e:.6g}'


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--eps', type=float, default=None,
                    help='fixed eps (default: adaptive from 1/500)')
    ap.add_argument('--verbose', '-v', action='store_true',
                    help='print all 5 check values per net')
    args = ap.parse_args()

    n_proved = n_failed = n_unproved = 0

    for line in sys.stdin:
        s = line.strip()
        if not s:
            continue
        try:
            poly = [tuple(int(x) for x in f.split(',')) for f in s.split(';')]
        except ValueError:
            print(f'PARSE ERROR: {s[:60]}', flush=True)
            continue

        proved, slack, eps_used = prove(poly, eps=args.eps, verbose=args.verbose)

        if math.isnan(slack):
            print(f'FAILED   (solver did not converge)', flush=True)
            n_failed += 1
        elif proved:
            print(f'PROVED   slack={slack:.4f}  eps={_eps_label(eps_used)}', flush=True)
            n_proved += 1
        else:
            print(f'UNPROVED slack={slack:.4f}  eps={_eps_label(eps_used)}', flush=True)
            n_unproved += 1

    total = n_proved + n_failed + n_unproved
    print(f'--- {total} nets: {n_proved} proved, {n_unproved} unproved, {n_failed} solver failures',
          file=sys.stderr)


if __name__ == '__main__':
    main()
