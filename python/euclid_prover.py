#!/usr/bin/env python3
"""
ieee_shape_prover_witness.py

Standalone Ellison-style verifier for a triangulated OBJ surface.

This version uses an NP-style witness/certificate design for the singular-value
bound.  Ordinary floating-point code proposes a factorization of

    A = J J^T - s^2 I.

The trusted verification code proves that A is positive definite by checking that
A is close enough to the proposed factorization C C^T.  NumPy/LAPACK are used
only to guess witnesses.  The accepted inequalities are checked by small
IEEE-float interval routines in this file.
"""

from __future__ import annotations

import argparse
import itertools
import math
from dataclasses import dataclass
from typing import Iterable, List, Sequence, Tuple

try:
    import numpy as np
except Exception as exc:
    raise SystemExit("This script needs numpy for witness generation.") from exc

NEG_INF = -math.inf
POS_INF = math.inf


def down(x: float) -> float:
    return math.nextafter(float(x), NEG_INF)


def up(x: float) -> float:
    return math.nextafter(float(x), POS_INF)


@dataclass(frozen=True)
class I:
    lo: float
    hi: float

    def __post_init__(self) -> None:
        if math.isnan(self.lo) or math.isnan(self.hi):
            raise ArithmeticError("NaN interval endpoint")
        if self.lo > self.hi:
            raise ArithmeticError(f"invalid interval [{self.lo}, {self.hi}]")

    @staticmethod
    def point(x: float) -> "I":
        x = float(x)
        if math.isnan(x):
            raise ArithmeticError("NaN point interval")
        return I(x, x)

    def __add__(self, other: "I") -> "I":
        return I(down(self.lo + other.lo), up(self.hi + other.hi))

    def __sub__(self, other: "I") -> "I":
        return I(down(self.lo - other.hi), up(self.hi - other.lo))

    def __neg__(self) -> "I":
        return I(-self.hi, -self.lo)

    def __mul__(self, other: "I") -> "I":
        vals = [
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        ]
        return I(down(min(vals)), up(max(vals)))

    def __truediv__(self, other: "I") -> "I":
        if other.lo <= 0.0 <= other.hi:
            raise ArithmeticError(f"division by interval containing zero: {other}")
        vals = [
            self.lo / other.lo,
            self.lo / other.hi,
            self.hi / other.lo,
            self.hi / other.hi,
        ]
        return I(down(min(vals)), up(max(vals)))

    def square(self) -> "I":
        if self.lo <= 0.0 <= self.hi:
            return I(0.0, up(max(self.lo * self.lo, self.hi * self.hi)))
        vals = [self.lo * self.lo, self.hi * self.hi]
        return I(down(min(vals)), up(max(vals)))

    def sqrt(self) -> "I":
        if self.lo < 0.0:
            raise ArithmeticError(f"sqrt of interval with negative lower end: {self}")
        return I(down(math.sqrt(self.lo)), up(math.sqrt(self.hi)))


def izero() -> I:
    return I.point(0.0)


def isum(items: Iterable[I]) -> I:
    out = izero()
    for x in items:
        out = out + x
    return out


def iabs_upper(x: I) -> float:
    return up(max(abs(x.lo), abs(x.hi)))


def frobenius_upper_from_intervals(entries: Iterable[I]) -> float:
    total = izero()
    for x in entries:
        a = I.point(iabs_upper(x))
        total = total + a.square()
    return total.sqrt().hi


def parse_obj(path: str) -> Tuple[List[Tuple[float, float, float]], List[Tuple[int, int, int]]]:
    vertices: List[Tuple[float, float, float]] = []
    faces: List[Tuple[int, int, int]] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if parts[0] == "v":
                if len(parts) < 4:
                    raise ValueError(f"bad vertex line: {line.rstrip()}")
                vertices.append((float(parts[1]), float(parts[2]), float(parts[3])))
            elif parts[0] == "f":
                if len(parts) != 4:
                    raise ValueError("only triangular OBJ faces are supported")
                face = []
                for p in parts[1:]:
                    token = p.split("/")[0]
                    idx = int(token)
                    if idx == 0:
                        raise ValueError("OBJ index 0 is invalid")
                    if idx < 0:
                        idx = len(vertices) + idx + 1
                    face.append(idx - 1)
                faces.append(tuple(face))  # type: ignore[arg-type]
    if not vertices or not faces:
        raise ValueError("OBJ must contain vertices and triangular faces")
    for f in faces:
        if len(set(f)) != 3:
            raise ValueError(f"degenerate face with repeated vertex: {f}")
        if min(f) < 0 or max(f) >= len(vertices):
            raise ValueError(f"face index out of range: {f}")
    return vertices, faces


def edges_from_faces(faces: Sequence[Tuple[int, int, int]]) -> List[Tuple[int, int]]:
    edges = set()
    for a, b, c in faces:
        edges.add(tuple(sorted((a, b))))
        edges.add(tuple(sorted((b, c))))
        edges.add(tuple(sorted((c, a))))
    return sorted(edges)


def build_jacobian(vertices: Sequence[Tuple[float, float, float]], edges: Sequence[Tuple[int, int]]) -> np.ndarray:
    n = len(vertices)
    j = np.zeros((len(edges), 3 * n), dtype=float)
    v = np.array(vertices, dtype=float)
    for r, (a, b) in enumerate(edges):
        d = 2.0 * (v[a] - v[b])
        j[r, 3 * a:3 * a + 3] = d
        j[r, 3 * b:3 * b + 3] = -d
    return j


def interval_length_squared(vertices: Sequence[Tuple[float, float, float]], a: int, b: int) -> I:
    terms = []
    for k in range(3):
        d = I.point(vertices[a][k]) - I.point(vertices[b][k])
        terms.append(d.square())
    return isum(terms)


def certify_rho_upper(vertices: Sequence[Tuple[float, float, float]], edges: Sequence[Tuple[int, int]]) -> float:
    one = I.point(1.0)
    terms = []
    for a, b in edges:
        err = interval_length_squared(vertices, a, b) - one
        terms.append(err.square())
    total = isum(terms)
    return total.sqrt().hi


def interval_jjt_minus_s2i(j: np.ndarray, s: float) -> List[List[I]]:
    m = j.shape[0]
    out: List[List[I]] = [[izero() for _ in range(m)] for _ in range(m)]
    s2 = I.point(s).square()
    for a in range(m):
        rowa = j[a]
        for b in range(a + 1):
            rowb = j[b]
            terms = []
            for k in range(j.shape[1]):
                if rowa[k] != 0.0 and rowb[k] != 0.0:
                    terms.append(I.point(rowa[k]) * I.point(rowb[k]))
            val = isum(terms)
            if a == b:
                val = val - s2
            out[a][b] = val
            out[b][a] = val
    return out


def float_a_matrix(j: np.ndarray, s: float) -> np.ndarray:
    a = j @ j.T
    a = (a + a.T) * 0.5
    idx = np.diag_indices_from(a)
    a[idx] -= s * s
    return a


def float_ldl_nopivot(a: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    n = a.shape[0]
    l = np.eye(n, dtype=float)
    d = np.zeros(n, dtype=float)
    for k in range(n):
        val = float(a[k, k])
        for p in range(k):
            val -= l[k, p] * l[k, p] * d[p]
        if not (val > 0.0 and math.isfinite(val)):
            raise np.linalg.LinAlgError(f"nonpositive LDL pivot at {k}: {val}")
        d[k] = val
        for i in range(k + 1, n):
            num = float(a[i, k])
            for p in range(k):
                num -= l[i, p] * l[k, p] * d[p]
            l[i, k] = num / val
    return l, d


def witness_factor(a_float: np.ndarray, method: str) -> Tuple[np.ndarray, str]:
    if method == "cholesky":
        c = np.linalg.cholesky(a_float)
        return c, "Cholesky witness C with A approx C C^T"
    if method == "ldl":
        l, d = float_ldl_nopivot(a_float)
        c = l.copy()
        for k in range(len(d)):
            c[:, k] *= math.sqrt(d[k])
        return c, "LDL witness converted to C = L sqrt(D)"
    raise ValueError(f"unknown factorization method: {method}")


def interval_product_entry(x: np.ndarray, y: np.ndarray, i: int, j: int) -> I:
    return isum(I.point(x[i, k]) * I.point(y[k, j]) for k in range(x.shape[1]))


def interval_cct_entry(c: np.ndarray, i: int, j: int) -> I:
    return isum(I.point(c[i, k]) * I.point(c[j, k]) for k in range(c.shape[1]))


def certify_factor_witness(a_interval: List[List[I]], c: np.ndarray) -> Tuple[bool, str, float, float, float]:
    n = c.shape[0]
    try:
        b = np.linalg.inv(c)
    except np.linalg.LinAlgError:
        return False, "witness factor is numerically singular", POS_INF, 0.0, POS_INF

    identity_errors: List[I] = []
    for i in range(n):
        for j in range(n):
            bc = interval_product_entry(b, c, i, j)
            target = I.point(1.0 if i == j else 0.0)
            identity_errors.append(target - bc)
    delta = frobenius_upper_from_intervals(identity_errors)
    if not (delta < 1.0):
        return False, f"could not certify inverse witness: ||I - B C||_F <= {delta:.17g}", delta, 0.0, POS_INF

    bnorm = frobenius_upper_from_intervals(I.point(float(x)) for x in b.ravel())
    factor_sigma_lower = down((1.0 - delta) / bnorm)
    if not (factor_sigma_lower > 0.0):
        return False, "factor singular-value lower bound is not positive", delta, factor_sigma_lower, POS_INF

    residuals: List[I] = []
    for i in range(n):
        for j in range(n):
            residuals.append(a_interval[i][j] - interval_cct_entry(c, i, j))
    residual = frobenius_upper_from_intervals(residuals)
    margin = down(factor_sigma_lower * factor_sigma_lower - residual)
    if margin > 0.0:
        msg = (
            f"witness verified: ||I-BC||_F <= {delta:.17g}, "
            f"sigma_min(C) >= {factor_sigma_lower:.17g}, "
            f"||A-CCT||_F <= {residual:.17g}, margin >= {margin:.17g}"
        )
        return True, msg, delta, factor_sigma_lower, residual
    msg = (
        f"witness margin failed: sigma_min(C)^2 lower={down(factor_sigma_lower * factor_sigma_lower):.17g}, "
        f"residual upper={residual:.17g}, delta={delta:.17g}"
    )
    return False, msg, delta, factor_sigma_lower, residual


def interval_ldl_positive_definite(a: List[List[I]]) -> Tuple[bool, float, str]:
    n = len(a)
    l: List[List[I]] = [[izero() for _ in range(n)] for _ in range(n)]
    d: List[I] = [izero() for _ in range(n)]
    min_pivot = POS_INF
    for k in range(n):
        correction = izero()
        for p in range(k):
            correction = correction + l[k][p].square() * d[p]
        pivot = a[k][k] - correction
        if not (pivot.lo > 0.0):
            return False, min_pivot, f"interval LDL pivot {k} not positive: [{pivot.lo}, {pivot.hi}]"
        d[k] = pivot
        min_pivot = min(min_pivot, pivot.lo)
        for i in range(k + 1, n):
            correction = izero()
            for p in range(k):
                correction = correction + l[i][p] * l[k][p] * d[p]
            l[i][k] = (a[i][k] - correction) / pivot
    return True, min_pivot, "all interval LDL pivots positive"


def certify_sigma_lower(
    j: np.ndarray,
    guess: float,
    method: str,
    shrink: float,
    max_tries: int,
) -> Tuple[float, str]:
    if not math.isfinite(guess) or guess <= 0.0:
        raise ValueError("singular-value guess must be positive")
    s = down(guess * shrink)
    failures = []
    for _ in range(max_tries):
        a_interval = interval_jjt_minus_s2i(j, s)
        if method == "interval-ldl":
            ok, min_pivot, msg = interval_ldl_positive_definite(a_interval)
            if ok:
                return s, f"interval LDL certificate: min pivot lower >= {min_pivot:.17g}"
            failures.append(f"s={s:.17g}: {msg}")
        else:
            try:
                a_float = float_a_matrix(j, s)
                c, witness_msg = witness_factor(a_float, method)
                ok, msg, _, _, _ = certify_factor_witness(a_interval, c)
                if ok:
                    return s, f"{witness_msg}; {msg}"
                failures.append(f"s={s:.17g}: {msg}")
            except (np.linalg.LinAlgError, ArithmeticError, ValueError) as exc:
                failures.append(f"s={s:.17g}: witness failed: {exc}")
        s = down(s * shrink)
    raise RuntimeError("could not certify sigma lower bound; last failures:\n" + "\n".join(failures[-8:]))


def all_simplices(nv: int, edges: Sequence[Tuple[int, int]], faces: Sequence[Tuple[int, int, int]]) -> List[Tuple[int, ...]]:
    simplices: List[Tuple[int, ...]] = [(i,) for i in range(nv)]
    simplices.extend(tuple(e) for e in edges)
    simplices.extend(tuple(f) for f in faces)
    return simplices


def active_set_closest_points(pverts: np.ndarray, qverts: np.ndarray) -> Tuple[float, np.ndarray, np.ndarray]:
    na = len(pverts)
    nb = len(qverts)
    best_d = POS_INF
    best_p = pverts[0].copy()
    best_q = qverts[0].copy()
    for mask_a in range(1, 1 << na):
        ia = [i for i in range(na) if (mask_a >> i) & 1]
        a = pverts[ia]
        for mask_b in range(1, 1 << nb):
            ib = [i for i in range(nb) if (mask_b >> i) & 1]
            b = qverts[ib]
            ca = len(a)
            cb = len(b)
            gram = np.block([[a @ a.T, -a @ b.T], [-b @ a.T, b @ b.T]]) * 2.0
            mat = np.zeros((ca + cb + 2, ca + cb + 2), dtype=float)
            rhs = np.zeros(ca + cb + 2, dtype=float)
            mat[:ca + cb, :ca + cb] = gram
            mat[:ca, ca + cb] = 1.0
            mat[ca + cb, :ca] = 1.0
            mat[ca:ca + cb, ca + cb + 1] = 1.0
            mat[ca + cb + 1, ca:ca + cb] = 1.0
            rhs[ca + cb] = 1.0
            rhs[ca + cb + 1] = 1.0
            try:
                sol = np.linalg.solve(mat, rhs)
            except np.linalg.LinAlgError:
                sol = np.linalg.lstsq(mat, rhs, rcond=None)[0]
            lam = sol[:ca]
            mu = sol[ca:ca + cb]
            if lam.min(initial=0.0) >= -1e-9 and mu.min(initial=0.0) >= -1e-9:
                lam = np.maximum(lam, 0.0)
                mu = np.maximum(mu, 0.0)
                if lam.sum() == 0.0 or mu.sum() == 0.0:
                    continue
                lam = lam / lam.sum()
                mu = mu / mu.sum()
                p = lam @ a
                q = mu @ b
                d = float(np.linalg.norm(p - q))
                if d < best_d:
                    best_d = d
                    best_p = p
                    best_q = q
    return best_d, best_p, best_q


def interval_dot(n: Sequence[float], v: Sequence[float]) -> I:
    return isum(I.point(n[k]) * I.point(v[k]) for k in range(3))


def norm_upper(n: Sequence[float]) -> float:
    return isum(I.point(n[k]).square() for k in range(3)).sqrt().hi


def separation_lower_for_normal(
    vertices: Sequence[Tuple[float, float, float]],
    simplex_a: Tuple[int, ...],
    simplex_b: Tuple[int, ...],
    normal: Sequence[float],
) -> float:
    nup = norm_upper(normal)
    if not (nup > 0.0 and math.isfinite(nup)):
        return NEG_INF
    dots_a = [interval_dot(normal, vertices[i]) for i in simplex_a]
    dots_b = [interval_dot(normal, vertices[i]) for i in simplex_b]
    min_a_lo = min(x.lo for x in dots_a)
    max_b_hi = max(x.hi for x in dots_b)
    gap1 = down(min_a_lo - max_b_hi)
    min_b_lo = min(x.lo for x in dots_b)
    max_a_hi = max(x.hi for x in dots_a)
    gap2 = down(min_b_lo - max_a_hi)
    gap = max(gap1, gap2)
    if gap <= 0.0:
        return NEG_INF
    return down(gap / nup)


def certify_collision_distance_lower(
    vertices: Sequence[Tuple[float, float, float]],
    edges: Sequence[Tuple[int, int]],
    faces: Sequence[Tuple[int, int, int]],
) -> Tuple[float, Tuple[int, ...], Tuple[int, ...], str]:
    verts_np = np.array(vertices, dtype=float)
    simplices = all_simplices(len(vertices), edges, faces)
    best_lower = POS_INF
    best_pair: Tuple[Tuple[int, ...], Tuple[int, ...]] = ((), ())
    failures = []
    checked = 0
    for a, b in itertools.combinations(simplices, 2):
        if set(a).intersection(b):
            continue
        checked += 1
        av = verts_np[list(a)]
        bv = verts_np[list(b)]
        _, p, q = active_set_closest_points(av, bv)
        candidates = []
        diff = p - q
        if np.linalg.norm(diff) > 0.0:
            candidates.append(tuple(float(x) for x in diff))
        cdiff = av.mean(axis=0) - bv.mean(axis=0)
        if np.linalg.norm(cdiff) > 0.0:
            candidates.append(tuple(float(x) for x in cdiff))
        for va in av:
            for vb in bv:
                d = va - vb
                if np.linalg.norm(d) > 0.0:
                    candidates.append(tuple(float(x) for x in d))
        pair_lower = NEG_INF
        for normal in candidates:
            lb = separation_lower_for_normal(vertices, a, b, normal)
            if lb > pair_lower:
                pair_lower = lb
        if not (pair_lower > 0.0):
            failures.append((a, b))
            if len(failures) >= 10:
                break
        if pair_lower < best_lower:
            best_lower = pair_lower
            best_pair = (a, b)
    if failures:
        msg = "failed to certify separation for disjoint simplex pairs, first failures: " + repr(failures[:10])
        raise RuntimeError(msg)
    return best_lower, best_pair[0], best_pair[1], f"checked {checked} disjoint simplex pairs"



def vec_sub_i(a: Sequence[I], b: Sequence[I]) -> List[I]:
    return [a[k] - b[k] for k in range(3)]


def dot_i(a: Sequence[I], b: Sequence[I]) -> I:
    return isum(a[k] * b[k] for k in range(3))


def cross_i(a: Sequence[I], b: Sequence[I]) -> List[I]:
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]


def norm_i(a: Sequence[I]) -> I:
    return isum(x.square() for x in a).sqrt()


def normalize_i(a: Sequence[I], what: str) -> List[I]:
    n = norm_i(a)
    if not (n.lo > 0.0):
        raise ArithmeticError(f"cannot normalize interval vector for {what}: norm interval [{n.lo}, {n.hi}]")
    return [x / n for x in a]


def vec_sub_f(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return a - b


def normalize_f(a: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(a))
    if not (n > 0.0 and math.isfinite(n)):
        raise ArithmeticError("cannot normalize zero/nonfinite float vector")
    return a / n


def turning_angle_float(uprev: np.ndarray, ucur: np.ndarray, unext: np.ndarray) -> float:
    nprev = normalize_f(np.cross(uprev, ucur))
    nnext = normalize_f(np.cross(ucur, unext))
    tin = normalize_f(np.cross(nprev, ucur))
    tout = normalize_f(np.cross(nnext, ucur))
    y = float(np.dot(ucur, np.cross(tin, tout)))
    x = float(np.dot(tin, tout))
    return math.atan2(y, x)


def turning_sincos_interval(uprev: Sequence[I], ucur: Sequence[I], unext: Sequence[I]) -> Tuple[I, I]:
    nprev = normalize_i(cross_i(uprev, ucur), "previous great-circle normal")
    nnext = normalize_i(cross_i(ucur, unext), "next great-circle normal")
    tin = normalize_i(cross_i(nprev, ucur), "incoming tangent")
    tout = normalize_i(cross_i(nnext, ucur), "outgoing tangent")
    s = dot_i(ucur, cross_i(tin, tout))
    c = dot_i(tin, tout)
    return s, c


def complex_mul_i(a_re: I, a_im: I, b_re: I, b_im: I) -> Tuple[I, I]:
    return a_re * b_re - a_im * b_im, a_re * b_im + a_im * b_re

def vertex_cycles_from_oriented_faces(
    vertex_count: int,
    faces: Sequence[Tuple[int, int, int]],
    reverse_link: bool = True,
) -> List[List[int]]:
    nxt: List[dict[int, int]] = [dict() for _ in range(vertex_count)]
    for a, b, c in faces:
        directed = [(a, b, c), (b, c, a), (c, a, b)]
        for v, x, y in directed:
            if reverse_link:
                x, y = y, x
            old = nxt[v].get(x)
            if old is not None and old != y:
                raise ValueError(f"inconsistent oriented link at vertex {v + 1}: {x + 1}->{old + 1} and {x + 1}->{y + 1}")
            nxt[v][x] = y
    cycles: List[List[int]] = []
    for v in range(vertex_count):
        if not nxt[v]:
            raise ValueError(f"isolated vertex {v + 1}")
        start = min(nxt[v])
        cycle = [start]
        seen = {start}
        cur = start
        while True:
            if cur not in nxt[v]:
                raise ValueError(f"broken oriented link at vertex {v + 1}, neighbor {cur + 1}")
            cur = nxt[v][cur]
            if cur == start:
                break
            if cur in seen:
                raise ValueError(f"oriented link at vertex {v + 1} is not a single cycle")
            seen.add(cur)
            cycle.append(cur)
        if len(cycle) != len(nxt[v]):
            raise ValueError(f"oriented link at vertex {v + 1} is disconnected")
        cycles.append(cycle)
    return cycles


def total_turning_float(vertices: Sequence[Tuple[float, float, float]], v: int, cycle: Sequence[int]) -> Tuple[float, List[float]]:
    arr = np.array(vertices, dtype=float)
    dirs = [normalize_f(arr[a] - arr[v]) for a in cycle]
    angles = []
    for i in range(len(dirs)):
        theta = turning_angle_float(dirs[(i - 1) % len(dirs)], dirs[i], dirs[(i + 1) % len(dirs)])
        angles.append(theta)
    return float(sum(angles)), angles


def total_turning_sin_half_interval(
    vertices: Sequence[Tuple[float, float, float]],
    v: int,
    cycle: Sequence[int],
    motion_radius: float,
) -> I:
    boxes: List[List[I]] = []
    for x, y, z in vertices:
        boxes.append([
            I(down(x - motion_radius), up(x + motion_radius)),
            I(down(y - motion_radius), up(y + motion_radius)),
            I(down(z - motion_radius), up(z + motion_radius)),
        ])
    dirs = []
    for a in cycle:
        dirs.append(normalize_i(vec_sub_i(boxes[a], boxes[v]), f"direction {v + 1}->{a + 1}"))

    # The local exterior dihedral/geodesic turn theta_i is represented by
    # c_i + I s_i = exp(I theta_i).  Since each theta_i is taken in the
    # principal range (-pi, pi), cos(theta_i/2) is positive and can be
    # certified as sqrt((1+c_i)/2).  Multiplying the half-angle factors
    # gives exp(I T/2), where T is the total geodesic turning around the
    # spherical link.  After the embedding certificate, the link is embedded,
    # and its enclosed spherical area is 2*pi - T, so -2*pi < T < 2*pi.
    # Hence T > 0 iff sin(T/2) > 0.
    re = I.point(1.0)
    im = I.point(0.0)
    one = I.point(1.0)
    two = I.point(2.0)
    for i in range(len(dirs)):
        s, c = turning_sincos_interval(
            dirs[(i - 1) % len(dirs)],
            dirs[i],
            dirs[(i + 1) % len(dirs)],
        )
        half_cos_squared = (one + c) / two
        if not (half_cos_squared.lo > 0.0):
            raise ArithmeticError(
                f"cannot certify positive half-angle cosine at link corner {i + 1}: "
                f"(1+c)/2 interval=[{half_cos_squared.lo}, {half_cos_squared.hi}]"
            )
        ch = half_cos_squared.sqrt()
        sh = s / (two * ch)
        re, im = complex_mul_i(re, im, ch, sh)
    return im

def certify_undented(
    vertices: Sequence[Tuple[float, float, float]],
    faces: Sequence[Tuple[int, int, int]],
    motion_radius: float,
    reverse_link: bool = True,
    verbose: bool = False,
) -> Tuple[bool, float, int, str, List[Tuple[int, int, float, float, float]]]:
    cycles = vertex_cycles_from_oriented_faces(len(vertices), faces, reverse_link=reverse_link)
    records: List[Tuple[int, int, float, float, float]] = []
    min_lower = POS_INF
    min_vertex = -1
    failures = []
    for v, cyc in enumerate(cycles):
        approx = math.nan
        try:
            sin_half = total_turning_sin_half_interval(vertices, v, cyc, motion_radius)
            lower = sin_half.lo
            upper = sin_half.hi
        except ArithmeticError as exc:
            failures.append(f"vertex {v + 1}: sin(T/2) interval evaluation failed: {exc}")
            lower = NEG_INF
            upper = POS_INF

        if verbose:
            try:
                approx, _angle_centers = total_turning_float(vertices, v, cyc)
            except Exception:
                approx = math.nan

        records.append((v + 1, len(cyc), approx, lower, upper))
        if lower < min_lower:
            min_lower = lower
            min_vertex = v + 1
        if not (lower > 0.0):
            failures.append(
                f"vertex {v + 1}: sin(T/2) lower bound not positive; degree={len(cyc)}, "
                f"sin_half_interval=[{lower:.17g}, {upper:.17g}]"
            )
    if failures:
        return False, min_lower, min_vertex, "; ".join(failures[:8]), records
    return True, min_lower, min_vertex, (
        f"all {len(vertices)} vertex sin(T/2) lower bounds are positive under ccw-outside face convention"
    ), records

def upper_motion_bound(sigma_lower: float, rho_upper: float, edge_count: int) -> float:
    sqrt_e = up(math.sqrt(float(edge_count)))
    disc = down(sigma_lower * sigma_lower - up(16.0 * rho_upper * sqrt_e))
    if disc <= 0.0:
        return POS_INF
    root = down(math.sqrt(disc))
    numerator = up(sigma_lower - root)
    denominator = down(8.0 * sqrt_e)
    return up(numerator / denominator)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Witness-checking IEEE interval verifier for OBJ unit-edge embedded realization existence")
    parser.add_argument("obj", help="triangulated OBJ file")
    parser.add_argument("--sigma-method", choices=["cholesky", "ldl", "interval-ldl"], default="cholesky")
    parser.add_argument("--sigma-shrink", type=float, default=0.98)
    parser.add_argument("--max-sigma-tries", type=int, default=80)
    parser.add_argument("--sharp-collision", action="store_true", help="use the sharper 2*motion < collision_lower embedding test instead of Ellison's default sqrt(V)*motion < collision_lower test")
    parser.add_argument("--print-turning", action="store_true", help="print per-vertex approximate and certified total turning bounds")
    args = parser.parse_args(argv)

    vertices, faces = parse_obj(args.obj)
    edges = edges_from_faces(faces)
    j = build_jacobian(vertices, edges)

    sv = np.linalg.svd(j, compute_uv=False)
    sigma_guess = float(sv[-1])
    rho_upper = certify_rho_upper(vertices, edges)
    sigma_lower, sigma_msg = certify_sigma_lower(j, sigma_guess, args.sigma_method, args.sigma_shrink, args.max_sigma_tries)
    collision_lower, ca, cb, collision_msg = certify_collision_distance_lower(vertices, edges, faces)
    motion_upper = upper_motion_bound(sigma_lower, rho_upper, len(edges))
    undented_ok, turning_min_lower, turning_min_vertex, turning_msg, turning_records = certify_undented(
        vertices,
        faces,
        motion_upper,
        reverse_link=True,
        verbose=args.print_turning,
    )

    sqrt_e = up(math.sqrt(float(len(edges))))
    existence_rhs = down((sigma_lower * sigma_lower) / up(16.0 * sqrt_e))
    existence_ok = rho_upper < existence_rhs
    two_motion_value = up(2.0 * motion_upper)
    two_motion_ok = two_motion_value < collision_lower
    ellison_value = up(up(math.sqrt(float(len(vertices)))) * motion_upper)
    ellison_motion_ok = ellison_value < collision_lower
    embedding_ok = two_motion_ok if args.sharp_collision else ellison_motion_ok
    embedding_test_used = "sharp_two_motion" if args.sharp_collision else "ellison_sqrtV"
    accepted = existence_ok and embedding_ok and undented_ok

    print(f"vertices: {len(vertices)}")
    print(f"faces: {len(faces)}")
    print(f"edges: {len(edges)}")
    print(f"rho_upper: {rho_upper:.17g}")
    print(f"sigma_guess_uncertified: {sigma_guess:.17g}")
    print(f"sigma_method: {args.sigma_method}")
    print(f"sigma_lower_certified: {sigma_lower:.17g}")
    print(f"sigma_certificate: {sigma_msg}")
    print(f"existence_rhs_sigma2_over_16sqrtE: {existence_rhs:.17g}")
    print(f"existence_test: {'PASS' if existence_ok else 'FAIL'}")
    print(f"motion_upper: {motion_upper:.17g}")
    print(f"collision_lower: {collision_lower:.17g}")
    print(f"closest_certified_pair: {ca} {cb}")
    print(f"collision_certificate: {collision_msg}")
    print(f"two_motion_test: {'PASS' if two_motion_ok else 'FAIL'}  value={two_motion_value:.17g} < {collision_lower:.17g}")
    print(f"ellison_sqrtV_motion_test: {'PASS' if ellison_motion_ok else 'FAIL'}  value={ellison_value:.17g} < {collision_lower:.17g}")
    print(f"embedding_test_used: {embedding_test_used}")
    print("turning_convention: ccw_outside_faces")
    print("turning_certificate: sin_half_total_turning_positive_using_embedding_bound_-2pi_lt_T_lt_2pi")
    if args.print_turning:
        for vertex_id, degree, approx, lower, upper in turning_records:
            print(f"turning_vertex {vertex_id}: degree={degree} approx_total_turning_diagnostic={approx:.17g} sin_half_total_turning_interval=[{lower:.17g}, {upper:.17g}] status={'PASS' if lower > 0.0 else 'FAIL'}")
    print(f"sin_half_total_turning_min_lower: {turning_min_lower:.17g}")
    print(f"turning_min_vertex: {turning_min_vertex}")
    print(f"undented_test: {'PASS' if undented_ok else 'FAIL'}  {turning_msg}")
    print(f"final: {'ACCEPT' if accepted else 'REJECT'}")
    return 0 if accepted else 2


if __name__ == "__main__":
    raise SystemExit(main())
