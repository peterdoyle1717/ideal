#!/usr/bin/env python3
"""
Ideal horoball packing: Newton solver for vertex weights u[v].

Given a triangulated sphere (list of oriented triangles), find weights u[v] such
that the angle sum around every interior vertex equals 2π (or 2π − defect).

Edge length between adjacent vertices i, j is u[i]*u[j]  (the product metric).
Vertex v0 = poly[0][0] is the "point at infinity"; its neighbors are pinned to 1.
Interior vertices (not adjacent to v0) are the free variables.

Newton's method starts from u = 1 everywhere.  A back-tracking line search
guards against triangle-inequality violations (which would make arccos blow up).

Interface
---------
    u = horou(poly)               # dict: vertex -> weight,  v0 -> float('inf')
    u = horou(poly, defect=d)     # with uniform angular defect d at interior verts
"""

import sys
import numpy as np

TAU = 2.0 * np.pi


# ── triangulation helpers ─────────────────────────────────────────────────────

def _cyclic_nbrs(v, faces):
    """Return the neighbors of v in cyclic order (either CW or CCW consistently)."""
    nxt = {}
    for fi, fj, fk in faces:
        if   fi == v: nxt[fk] = fj
        elif fj == v: nxt[fi] = fk
        elif fk == v: nxt[fj] = fi
    if not nxt:
        return []
    start = next(iter(nxt))
    ring = [start]
    cur = nxt[start]
    while cur != start:
        ring.append(cur)
        cur = nxt[cur]
    return ring


# ── angle / flower ────────────────────────────────────────────────────────────

def _petal(ui, uj, uk):
    """Angle at the ui-vertex of a triangle with edge lengths ui*uj, ui*uk, uj*uk."""
    a, b, c = ui * uj, ui * uk, uj * uk
    cos_t = (a*a + b*b - c*c) / (2.0 * a * b)
    return np.arccos(np.clip(cos_t, -1.0, 1.0))


def _petal_grad(ui, uj, uk):
    """Gradients of petal(ui,uj,uk) w.r.t. each argument.

    d(petal)/d(ui) = -(uj*uk / ui^3) / sin(theta)
    d(petal)/d(uj) = -(1/(2uk) - uk/(2uj^2) - uk/(2ui^2)) / sin(theta)
    d(petal)/d(uk) = -(1/(2uj) - uj/(2uk^2) - uj/(2ui^2)) / sin(theta)

    Returns (dui, duj, duk).  Returns (0,0,0) at degenerate triangles.
    """
    a, b, c = ui * uj, ui * uk, uj * uk
    cos_t = (a*a + b*b - c*c) / (2.0 * a * b)
    sin_t = np.sqrt(max(0.0, 1.0 - np.clip(cos_t, -1.0, 1.0)**2))
    if sin_t < 1e-15:
        return 0.0, 0.0, 0.0
    s = -1.0 / sin_t
    dui = s * (uj * uk / ui**3)
    duj = s * (1.0/(2*uk) - uk/(2*uj**2) - uk/(2*ui**2))
    duk = s * (1.0/(2*uj) - uj/(2*uk**2) - uj/(2*ui**2))
    return dui, duj, duk


def _flower(ui, ring_u):
    """Sum of angles around the ui-vertex; ring_u is the cyclic neighbor u-values."""
    k = len(ring_u)
    return sum(_petal(ui, ring_u[j], ring_u[(j + 1) % k]) for j in range(k))


# ── Newton solver ─────────────────────────────────────────────────────────────

def horou(poly, defect=0.0, tol=1e-12, max_iter=200):
    """
    Solve for horoball weights on a triangulated sphere.

    Parameters
    ----------
    poly     : list of (a, b, c) oriented triangles (integer vertex indices).
               poly[0][0] is taken as the point at infinity (v0).
    defect   : uniform angular defect target at interior vertices (default 0).
               Equation: flower(v) = 2π − defect.
    tol      : convergence tolerance on max |residual|.
    max_iter : Newton iteration cap.

    Returns
    -------
    dict  vertex -> u-value.  v0 maps to float('inf'), boundary vertices to 1.0.
    """
    v0 = poly[0][0]
    all_verts = sorted({v for f in poly for v in f})

    bndry    = set(_cyclic_nbrs(v0, poly))
    interior = [v for v in all_verts if v != v0 and v not in bndry]
    n        = len(interior)
    vidx     = {v: i for i, v in enumerate(interior)}

    bndry_vals = {v: 1.0 for v in bndry}

    if n == 0:
        return {v0: float('inf'), **bndry_vals}

    # Precompute cyclic neighbor rings for every interior vertex.
    nbr_lists = {v: _cyclic_nbrs(v, poly) for v in interior}

    # Faces that don't touch v0 — needed for the triangle-inequality guard.
    finite_faces = [(a, b, c) for a, b, c in poly if v0 not in (a, b, c)]

    target = TAU - defect

    def u_val(v, x):
        """u-value for vertex v given current interior solution x."""
        i = vidx.get(v)
        return x[i] if i is not None else 1.0   # boundary fixed at 1

    def residual(x):
        F = np.empty(n)
        for i, v in enumerate(interior):
            ring_u = [u_val(nb, x) for nb in nbr_lists[v]]
            F[i] = _flower(x[i], ring_u) - target
        return F

    def jacobian(x):
        J = np.zeros((n, n))
        for i, v in enumerate(interior):
            ns = nbr_lists[v]; k = len(ns)
            ui = x[i]
            for j in range(k):
                vj = ns[j]; vk = ns[(j + 1) % k]
                uj = u_val(vj, x); uk = u_val(vk, x)
                dui, duj, duk = _petal_grad(ui, uj, uk)
                J[i, i] += dui
                if vj in vidx: J[i, vidx[vj]] += duj
                if vk in vidx: J[i, vidx[vk]] += duk
        return J

    def min_tri_slack(x):
        """Minimum triangle-inequality slack over all finite faces."""
        s = float('inf')
        for a, b, c in finite_faces:
            ua, ub, uc = u_val(a, x), u_val(b, x), u_val(c, x)
            p, q, r = ua*ub, ua*uc, ub*uc
            s = min(s, p+q-r, p+r-q, q+r-p)
        return s

    x = np.ones(n)

    for _ in range(max_iter):
        F = residual(x)
        res = np.max(np.abs(F))
        if res < tol:
            break

        J = jacobian(x)
        try:
            dx = np.linalg.solve(J, -F)
        except np.linalg.LinAlgError:
            dx = np.linalg.lstsq(J, -F, rcond=None)[0]

        # Back-tracking line search.
        # Requirements: all u > 0, no triangle-inequality violation, residual decreases.
        step = 1.0
        for _ in range(60):
            x_new = x + step * dx
            if (x_new.min() > 0
                    and min_tri_slack(x_new) > 0
                    and np.max(np.abs(residual(x_new))) < res):
                break
            step *= 0.5
        else:
            break   # can't make progress; return best so far

        x = x_new

    result = {v0: float('inf'), **bndry_vals}
    for i, v in enumerate(interior):
        result[v] = x[i]
    return result


# ── CLI ───────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    for line in sys.stdin:
        s = line.strip()
        if not s:
            continue
        poly = [tuple(int(x) for x in f.split(',')) for f in s.split(';')]
        u    = horou(poly)
        v0   = poly[0][0]
        vals = sorted((v, u[v]) for v in u if v != v0)
        parts = ' '.join(f'{v}:{w:.8g}' for v, w in vals)
        print(parts)
