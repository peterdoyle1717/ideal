# Closed-form base-bend imputation in puffup

`puffup` solves for the per-edge dihedral bends of a unit-equilateral
Euclidean realization of a triangulated sphere by homotopy in the
corner angle α from 0 (the ideal-horoball limit, where bends are read
off `horou`'s `u`-values via `dihedral`) to π/3 (Euclidean unit
triangles).

The homotopy fixes a base face's three vertices at a canonical gauge
and Newton-solves the bends on the **non-base** edges (`3v − 9`
variables) against the holonomy residual at the **non-base** vertices
(`3 · (v − 3) = 3v − 9` equations). The three bends on the base-face
edges are not touched during the homotopy. After the loose-tol
homotopy and the tight final Newton at α = π/3, those three base
bends are recovered in closed form by the construction below.

## Setup

Let `Y(a)` denote yaw about the z-axis and `P(b)` the pitch
matx(−b):

```
Y(a) = | cos a   −sin a   0 |        P(b) = |  1     0      0    |
       | sin a    cos a   0 |               |  0   cos b  sin b  |
       |   0       0      1 |               |  0  −sin b  cos b  |
```

`P(b)` is exactly `matx(-b)` from the holonomy word, so `movemat(α, β) = Y(α) · P(β)`.

## Holonomy at a base vertex

Walk the flower around base vertex `v ∈ base = (v0, v1, v2)`. The
flower starts at the base face and ends at the base face's wrap-around
neighbor; the **first** factor uses the bend on the edge from `v` to
its preceding base neighbor; the **last** factor uses the bend on the
edge from `v` to its succeeding base neighbor. The middle `d − 2`
factors use bends on non-base edges. Calling those bends
`b_first, b_n1, …, b_n(d-2), b_last`:

```
M_v = Y(α) P(b_first) · [ Y(α) P(b_n1) · … · Y(α) P(b_n(d-2)) ] · Y(α) P(b_last)
    = Y(α) P(b_first) · K · Y(α) P(b_last)        ← define K
    = I                                            ← closure
```

where `K` is the product of the d − 2 non-base movemats, exactly the
ones whose bends were already solved during the homotopy.

## The closing block

Pre-multiply both sides by `Y(−α)`:

```
P(b_first) · K · Y(α) · P(b_last) = Y(−α)
P(b_first) · X · P(b_last) = Y(−α)            where X = K · Y(α)
```

From this we read off `X` algebraically:

```
X = P(−b_first) · Y(−α) · P(−b_last)
```

Taking the transpose (= inverse, since `X` is a rotation):

```
X^T = P(b_last) · Y(α) · P(b_first)
```

So `X^T` is in the canonical form `P(A) Y(g) P(B)` with `A = b_last`,
`g = α`, `B = b_first`.

## Recovering A from `P(A) Y(g) P(B)`

Direct expansion of `P(A) Y(g) P(B)` gives:

```
M[0,0] = cos g
M[1,0] =  sin g · cos A      M[2,0] = −sin g · sin A
M[0,1] = −sin g · cos B      M[0,2] = −sin g · sin B
```

For `0 < g < π` (i.e. `sin g > 0`):

```
A = atan2(−M[2,0],  M[1,0])
B = atan2(−M[0,2], −M[0,1])
```

Applied with `M = X^T` (so `M[2,0] = X[0,2]`, `M[1,0] = X[0,1]`, etc.):

```
b_last  = atan2(−X[0,2],  X[0,1])
b_first = atan2(−X[2,0], −X[1,0])
```

`b_last` is the bend on the edge from base vertex `v` to its
**succeeding** base neighbor. Each base bend is recovered exactly
once at the appropriate base vertex (`base[i] → b_last` becomes
`bend(base[i], base[(i+1) mod 3])`). The corresponding `b_first`
recovered at the same vertex is the bend on the **preceding** base
edge — the same value that vertex `(i − 1) mod 3` will produce as its
own `b_last`. Comparing the two is a free consistency check.

## Algorithm

```
for each base vertex v in (base[0], base[1], base[2]):
    K ← product of the (d − 2) non-base movemats around v
    X ← K · Y(α)
    A ← atan2(−X[0,2], X[0,1])
    bend on edge (v, next_base_vertex) ← A
```

No iteration. No tolerance. The recovered base bends close every
base-vertex holonomy to within the trig-cancellation floor of the
3×3 matrix product (~1e−14 in double precision for v ≤ 50).

## References

- C implementation: `ideal/src/puffup_c.c` (function
  `complete_base_bends`).
- Python implementation: `ideal/python/puffup.py` (function
  `_complete_base_bends`).
