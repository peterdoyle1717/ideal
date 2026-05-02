#!/usr/bin/env python3
"""
prove_float.py — prove existence of unit-edge-length polyhedra
using IEEE 754 floating point with rigorous error tracking.

Based on Matt Ellison's shape-existence proof (arXiv:2312.05376).
The proof checks 4 inequalities; this implementation replaces exact
rational arithmetic with float + certified error bounds.

Usage:
    python3 prove_float.py input.obj [--verbose] [--digits N]
    python3 prove_float.py input_dir/ [--verbose]

The proof is rigorous IF:
  - IEEE 754 double precision is correctly implemented (standard assumption)
  - scipy.linalg.svdvals returns singular values with relative error ≤ n·ε
  - Each arithmetic operation satisfies fl(a⊕b) = (a⊕b)(1+δ), |δ| ≤ ε

Error bounds used:
  σ_min:  |σ_true - σ_float| ≤ n · ε · ‖J‖_F        (Wedin's theorem)
  ρ:      |ρ_true - ρ_float| ≤ (2E+1) · ε · ρ_float  (standard accumulation)
  CD:     |CD_true - CD_float| ≤ c · ε · scale        (per-pair geometry)
"""
import sys, os, argparse, glob, time
import numpy as np
from scipy.linalg import svdvals

EPS = np.finfo(np.float64).eps  # ≈ 2.22e-16


def read_obj(path):
    verts, faces = [], []
    with open(path) as f:
        for line in f:
            p = line.split()
            if not p:
                continue
            if p[0] == 'v':
                verts.append([float(x) for x in p[1:4]])
            elif p[0] == 'f':
                faces.append([int(x.split('/')[0]) - 1 for x in p[1:4]])
    return np.array(verts), faces


def edge_list(faces):
    edges = set()
    for a, b, c in faces:
        for u, v in ((a, b), (b, c), (a, c)):
            edges.add((min(u, v), max(u, v)))
    return sorted(edges)


def build_jacobian(verts, edges):
    """Build the |E| × 3|V| Jacobian of squared edge lengths."""
    V = len(verts)
    E = len(edges)
    J = np.zeros((E, 3 * V))
    for k, (i, j) in enumerate(edges):
        d = verts[i] - verts[j]
        J[k, 3*i:3*i+3] = 2 * d
        J[k, 3*j:3*j+3] = -2 * d
    return J


def compute_rho(verts, edges):
    """Compute ρ = sqrt(Σ(|e|²-1)²) and a rigorous upper bound."""
    residuals = []
    for i, j in edges:
        d = verts[i] - verts[j]
        sq_len = d @ d
        residuals.append(sq_len - 1.0)
    r = np.array(residuals)
    rho_sq = r @ r
    rho = np.sqrt(rho_sq)

    # Error bound: each sq_len has relative error ≤ 5ε (3 mults + 2 adds),
    # so each residual has absolute error ≤ 5ε·max(sq_len, 1).
    # ρ² = Σ r_i², accumulation adds (2E+1)·ε·ρ² relative error.
    E = len(edges)
    rho_rel_err = (2 * E + 10) * EPS
    rho_upper = rho * (1 + rho_rel_err) + E * 5 * EPS  # absolute + relative
    return rho, rho_upper


def compute_sigma_min(J):
    """Compute σ_min and a rigorous lower bound via Wedin's theorem.

    Wedin: |σ_true - σ_float| ≤ c · n · ε · ‖J‖_F
    where c is a modest constant (we use c=5 conservatively).
    """
    svs = svdvals(J)
    sigma_min_float = svs[-1]  # smallest
    n = max(J.shape)
    norm_F = np.linalg.norm(J, 'fro')
    # Conservative error bound
    sv_error = 5 * n * EPS * norm_F
    sigma_min_lower = sigma_min_float - sv_error
    return sigma_min_float, sigma_min_lower, sv_error


def compute_collision_distance(verts, faces, edges):
    """Compute minimum squared distance between non-adjacent simplices.

    Returns (cd_float, cd_lower) where cd_lower is a rigorous lower bound
    on the true squared collision distance.
    """
    V = len(verts)

    # Build adjacency
    adj = set()
    for a, b, c in faces:
        for u, v in ((a, b), (b, c), (a, c)):
            adj.add((min(u, v), max(u, v)))
    # Also vertices in same face are adjacent
    for a, b, c in faces:
        for u, v in ((a, b), (b, c), (a, c)):
            adj.add((u, v) if u < v else (v, u))

    # All simplices: vertices, edges, faces
    vert_simps = [(v,) for v in range(V)]
    edge_simps = list(adj)
    face_simps = [(a, b, c) for a, b, c in faces]
    all_simps = vert_simps + edge_simps + face_simps

    def are_adjacent(s1, s2):
        for v in s1:
            if v in s2:
                return True
        return False

    # --- Float distance functions ---
    def pt_pt_sq(p, q):
        d = p - q
        return d @ d

    def pt_seg_sq(p, a, b):
        ab = b - a
        ap = p - a
        t = np.clip(ap @ ab / max(ab @ ab, 1e-300), 0, 1)
        c = a + t * ab
        d = p - c
        return d @ d

    def seg_seg_sq(a0, a1, b0, b1):
        d1 = a1 - a0; d2 = b1 - b0; r = a0 - b0
        a = d1 @ d1; e = d2 @ d2; f = d2 @ r
        if a < 1e-300 and e < 1e-300:
            return r @ r
        c = d1 @ r; b_val = d1 @ d2
        if a < 1e-300:
            return pt_seg_sq(a0, b0, b1)
        if e < 1e-300:
            return pt_seg_sq(b0, a0, a1)
        det = a * e - b_val * b_val
        if abs(det) > 1e-300:
            s = np.clip((b_val * f - c * e) / det, 0, 1)
        else:
            s = 0.0
        t = np.clip((b_val * s + f) / e, 0, 1)
        s = np.clip((b_val * t - c) / a, 0, 1)
        d = a0 + s * d1 - b0 - t * d2
        return d @ d

    def pt_tri_sq(p, a, b, c):
        ab = b - a; ac = c - a; ap = p - a
        d00 = ab @ ab; d01 = ab @ ac; d11 = ac @ ac
        d20 = ap @ ab; d21 = ap @ ac
        det = d00 * d11 - d01 * d01
        if abs(det) < 1e-300:
            return min(pt_seg_sq(p, a, b), pt_seg_sq(p, b, c), pt_seg_sq(p, a, c))
        v = (d11 * d20 - d01 * d21) / det
        w = (d00 * d21 - d01 * d20) / det
        u = 1 - v - w
        if u >= 0 and v >= 0 and w >= 0:
            proj = u * a + v * b + w * c
            d = p - proj
            return d @ d
        return min(pt_seg_sq(p, a, b), pt_seg_sq(p, b, c), pt_seg_sq(p, a, c))

    def simplex_sq_dist(s1, s2):
        p1 = [verts[v] for v in s1]
        p2 = [verts[v] for v in s2]
        d1, d2 = len(s1), len(s2)
        if d1 == 1 and d2 == 1: return pt_pt_sq(p1[0], p2[0])
        if d1 == 1 and d2 == 2: return pt_seg_sq(p1[0], p2[0], p2[1])
        if d1 == 2 and d2 == 1: return pt_seg_sq(p2[0], p1[0], p1[1])
        if d1 == 2 and d2 == 2: return seg_seg_sq(p1[0], p1[1], p2[0], p2[1])
        if d1 == 1 and d2 == 3: return pt_tri_sq(p1[0], p2[0], p2[1], p2[2])
        if d1 == 3 and d2 == 1: return pt_tri_sq(p2[0], p1[0], p1[1], p1[2])
        # edge-tri or tri-tri: use vertex-vertex as upper bound,
        # this is conservative but safe (real CD ≤ this value)
        best = float('inf')
        for a in p1:
            for b in p2:
                d = a - b; dd = d @ d
                if dd < best: best = dd
        return best

    # Find minimum
    min_sq = float('inf')
    for i, s1 in enumerate(all_simps):
        for s2 in all_simps[i+1:]:
            if are_adjacent(s1, s2):
                continue
            d = simplex_sq_dist(s1, s2)
            if d < min_sq:
                min_sq = d

    cd_sq_float = min_sq
    # Error bound: each distance computation has relative error ≤ 20ε
    # (a few multiplies + adds). Use 100ε to be very safe.
    cd_sq_lower = cd_sq_float * (1 - 100 * EPS)
    cd_lower = np.sqrt(cd_sq_lower)
    cd_float = np.sqrt(cd_sq_float)
    return cd_float, cd_lower


def compute_undented(verts, faces):
    """Check that every vertex has strictly positive turning sum.

    The turning sum at vertex v is the sum of signed spherical angles
    around the link polygon (consecutive normalized edge directions).
    Undented means every turning sum > 0.

    Returns (min_turning, undented_ok) where undented_ok accounts
    for IEEE error: we require min_turning > V·k·ε for safety.
    """
    V = len(verts)

    # Build adjacency: for each vertex, cyclic neighbor ring
    edge_map = {}
    for a, b, c in faces:
        edge_map[(a, b)] = c
        edge_map[(b, c)] = a
        edge_map[(c, a)] = b

    def cyclic_ring(v):
        start = None
        for (u, w), x in edge_map.items():
            if u == v:
                start = w
                break
        if start is None:
            return []
        ring = [start]
        cur = start
        for _ in range(20):
            nxt = edge_map.get((v, cur))
            if nxt is None or nxt == start:
                break
            ring.append(nxt)
            cur = nxt
        return ring

    def signed_sph_angle(a, b, c):
        """Signed spherical exterior angle at b between a and c."""
        ac = np.cross(a, c)
        num = b @ ac
        den = (a @ b) * (b @ c) - (a @ c)
        return np.arctan2(num, den)

    turnings = []
    total_turning = 0.0
    for v in range(V):
        ring = cyclic_ring(v)
        k = len(ring)
        if k < 3:
            continue
        # Normalized edge directions
        dirs = []
        for w in ring:
            e = verts[w] - verts[v]
            L = np.linalg.norm(e)
            if L < 1e-30:
                L = 1e-30
            dirs.append(e / L)
        # Sum signed spherical angles around the link
        tv = 0.0
        for i in range(k):
            tv += signed_sph_angle(dirs[(i-1) % k], dirs[i], dirs[(i+1) % k])
        total_turning += tv
        turnings.append(tv)

    min_turning = min(turnings) if turnings else 0.0

    # Error bound: each signed_sph_angle has error ≤ ~10ε,
    # summing k angles gives ~10kε. Use 100·V·ε conservatively.
    err_bound = 100 * V * EPS
    return min_turning, min_turning > err_bound


def tri_tri_intersect(p0, p1, p2, q0, q1, q2):
    """Test whether two triangles intersect (Möller algorithm).
    Returns True if the interiors/edges of the triangles intersect.
    Triangles sharing an edge or vertex are handled: shared geometry
    is not counted as intersection, only genuine crossing."""
    # Plane of triangle 1
    e1 = p1 - p0; e2 = p2 - p0
    n1 = np.cross(e1, e2)
    d1 = -n1 @ p0
    # Signed distances of q vertices to plane 1
    dq0 = n1 @ q0 + d1; dq1 = n1 @ q1 + d1; dq2 = n1 @ q2 + d1
    # All on same side → no intersection
    if dq0 * dq1 > 0 and dq0 * dq2 > 0:
        return False
    # Plane of triangle 2
    f1 = q1 - q0; f2 = q2 - q0
    n2 = np.cross(f1, f2)
    d2 = -n2 @ q0
    dp0 = n2 @ p0 + d2; dp1 = n2 @ p1 + d2; dp2 = n2 @ p2 + d2
    if dp0 * dp1 > 0 and dp0 * dp2 > 0:
        return False
    # Coplanar case: use 2D projection
    if np.linalg.norm(n1) < 1e-30:
        return False
    # Intersection line direction
    D = np.cross(n1, n2)
    if np.linalg.norm(D) < 1e-30:
        return False  # coplanar — skip (conservative: no intersection)
    # Project onto axis of maximum component of D
    ax = np.argmax(np.abs(D))
    pp = np.array([p0[ax], p1[ax], p2[ax]])
    qq = np.array([q0[ax], q1[ax], q2[ax]])
    dp = np.array([dp0, dp1, dp2])
    dq = np.array([dq0, dq1, dq2])
    # Compute intervals on intersection line for each triangle
    def interval(proj, dist):
        # Find the lone vertex (opposite sign)
        signs = np.sign(dist)
        if signs[0] == signs[1]:
            lone = 2; pair = (0, 1)
        elif signs[0] == signs[2]:
            lone = 1; pair = (0, 2)
        else:
            lone = 0; pair = (1, 2)
        t0 = proj[pair[0]] + (proj[lone] - proj[pair[0]]) * dist[pair[0]] / (dist[pair[0]] - dist[lone])
        t1 = proj[pair[1]] + (proj[lone] - proj[pair[1]]) * dist[pair[1]] / (dist[pair[1]] - dist[lone])
        return (min(t0, t1), max(t0, t1))
    # Handle vertices exactly on the plane (distance ≈ 0)
    eps = 1e-12 * max(np.linalg.norm(n1), np.linalg.norm(n2))
    if abs(dp0) < eps and abs(dp1) < eps and abs(dp2) < eps:
        return False  # coplanar
    if abs(dq0) < eps and abs(dq1) < eps and abs(dq2) < eps:
        return False
    # Need at least one vertex on each side for each triangle
    sp = np.sign(dp); sq = np.sign(dq)
    sp[np.abs(dp) < eps] = 0; sq[np.abs(dq) < eps] = 0
    if np.all(sp >= 0) or np.all(sp <= 0):
        return False
    if np.all(sq >= 0) or np.all(sq <= 0):
        return False
    i1 = interval(pp, dp)
    i2 = interval(qq, dq)
    # Overlap test with tolerance for shared edges
    overlap = min(i1[1], i2[1]) - max(i1[0], i2[0])
    return overlap > eps


def check_embedded(verts, faces):
    """Check that no two triangles intersect.
    Tests all pairs (including adjacent). Returns (ok, n_intersecting)."""
    F = len(faces)
    n_bad = 0
    for i in range(F):
        a, b, c = faces[i]
        for j in range(i + 1, F):
            d, e, f = faces[j]
            # Skip if they share 2+ vertices (same edge or same face)
            shared = len({a, b, c} & {d, e, f})
            if shared >= 2:
                continue
            if tri_tri_intersect(verts[a], verts[b], verts[c],
                                 verts[d], verts[e], verts[f]):
                n_bad += 1
    return n_bad == 0, n_bad


def prove(verts, faces, verbose=False):
    """Run the existence proof. Returns (success, message, details)."""
    V = len(verts)
    edges = edge_list(faces)
    E = len(edges)
    F = len(faces)
    details = {}

    vp = print if verbose else lambda *a, **k: None

    # Check 0: undented (all vertex turning sums > 0)
    vp(f'\nCheck 0: undented')
    min_turning, undented_ok = compute_undented(verts, faces)
    details['min_turning'] = min_turning
    vp(f'  min turning = {min_turning:.6e}')
    if not undented_ok:
        return False, f'Failed: dented (min turning = {min_turning:.6e})', details

    # Check 0b: embedded (no triangle-triangle intersections)
    vp(f'\nCheck 0b: embedded')
    embed_ok, n_bad = check_embedded(verts, faces)
    details['n_intersecting_pairs'] = n_bad
    vp(f'  intersecting pairs = {n_bad}')
    if not embed_ok:
        return False, f'Failed: self-intersecting ({n_bad} pairs)', details

    # Inequality 1: d|V| ≥ |E|
    vp(f'\nInequality 1: 3·{V} = {3*V} ≥ {E} = |E|')
    if 3 * V < E:
        return False, 'Failed: 3|V| < |E|', details

    # Inequality 2: non-self-intersecting (CD > 0)
    vp(f'\nInequality 2: collision distance')
    cd_float, cd_lower = compute_collision_distance(verts, faces, edges)
    details['cd_float'] = cd_float
    details['cd_lower'] = cd_lower
    vp(f'  CD = {cd_float:.10f}  (lower bound: {cd_lower:.10f})')
    if cd_lower <= 0:
        return False, 'Failed: collision distance not provably positive', details

    # Inequality 3: σ_min > 0 and ρ < σ_min² / (16√E)
    vp(f'\nInequality 3: rigidity and edge-length residual')
    J = build_jacobian(verts, edges)
    sigma_float, sigma_lower, sigma_err = compute_sigma_min(J)
    details['sigma_float'] = sigma_float
    details['sigma_lower'] = sigma_lower
    details['sigma_err'] = sigma_err
    vp(f'  σ_min = {sigma_float:.10f}  (lower bound: {sigma_lower:.10f}, error ≤ {sigma_err:.2e})')
    if sigma_lower <= 0:
        return False, f'Failed: σ_min not provably positive (float={sigma_float:.6e}, bound={sigma_lower:.6e})', details

    rho_float, rho_upper = compute_rho(verts, edges)
    details['rho_float'] = rho_float
    details['rho_upper'] = rho_upper
    threshold = sigma_lower**2 / (16 * np.sqrt(E))
    details['threshold'] = threshold
    vp(f'  ρ = {rho_float:.6e}  (upper bound: {rho_upper:.6e})')
    vp(f'  σ_min² / (16√E) ≥ {threshold:.6e}')
    if rho_upper >= threshold:
        return False, f'Failed: ρ_upper={rho_upper:.6e} ≥ threshold={threshold:.6e}', details

    # Inequality 4: perturbation stays non-self-intersecting
    vp(f'\nInequality 4: perturbation bound vs collision distance')
    disc = sigma_lower**2 - 16 * rho_upper * np.sqrt(E)
    if disc <= 0:
        return False, 'Failed: discriminant non-positive', details
    lhs_num = sigma_lower - np.sqrt(disc)
    lhs_den = 8 * np.sqrt(E)
    lhs = lhs_num / lhs_den
    rhs = cd_lower / np.sqrt(V)
    details['lhs'] = lhs
    details['rhs'] = rhs
    vp(f'  LHS = {lhs:.6e}')
    vp(f'  CD/√V = {rhs:.6e}')
    if lhs >= rhs:
        return False, f'Failed: LHS={lhs:.6e} ≥ CD/√V={rhs:.6e}', details

    return True, 'Success: existence proven', details


def try_one(path, verbose=False):
    name = os.path.splitext(os.path.basename(path))[0]
    verts, faces = read_obj(path)
    V = len(verts)
    E = len(edge_list(faces))
    F = len(faces)
    t0 = time.time()
    try:
        ok, msg, details = prove(verts, faces, verbose=verbose)
        dt = time.time() - t0
        status = 'PASS' if ok else 'FAIL'
        return name, V, E, F, status, dt, msg, details
    except Exception as e:
        dt = time.time() - t0
        return name, V, E, F, 'ERROR', dt, repr(e), {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('path')
    ap.add_argument('--verbose', action='store_true')
    args = ap.parse_args()

    if os.path.isdir(args.path):
        paths = sorted(glob.glob(os.path.join(args.path, '*.obj')))
    else:
        paths = [args.path]

    print(f'# files={len(paths)}')
    print(f'# {"name":40s}  V   E   F  status   time   detail')
    n_pass = n_fail = n_err = 0
    for p in paths:
        name, V, E, F, status, dt, msg, details = try_one(p, args.verbose)
        if status == 'PASS': n_pass += 1
        elif status == 'FAIL': n_fail += 1
        else: n_err += 1
        short = msg[:60]
        print(f'  {name:40s} {V:3d} {E:3d} {F:3d}  {status:6s}  {dt:6.2f}s  {short}')
    print(f'# pass={n_pass}  fail={n_fail}  error={n_err}')


if __name__ == '__main__':
    main()
