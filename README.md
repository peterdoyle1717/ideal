# neo/ideal — Ideal Horoball Packing

Tools for computing ideal horoball packings of triangulated spheres and
proving their existence via a sub/supersolution bracket.

**Input:** face lists `a,b,c;d,e,f;...` (one triangulation per line),
as produced by `neo/clers/python/clers_decode.py`.

---

## What is an ideal horoball packing?

A triangulated sphere has an *ideal horoball packing* when each vertex carries
a weight `u[v] > 0` such that the angle sum around every interior vertex equals
2π. The edge length between adjacent vertices `i` and `j` is `u[i]*u[j]` (the
product metric). One vertex is sent to infinity; its neighbors are pinned to
`u = 1` (boundary). The remaining (interior) weights are the solution.

The packing always exists (proven here for all prime 6-nets through v = 100).

---

## Layout

```
ideal/
├── Makefile
├── data/
│   └── eps_needed.txt       # smallest eps proving all prime 6-nets, v=4..80
├── python/
│   ├── horou.py             # Newton solver for horoball weights u[v]
│   ├── horoz.py             # BFS placement of vertices in the upper half-plane
│   └── proof.py             # sub/supersolution existence proof (CLI + library)
├── scripts/
│   ├── run_proof_doob.sh    # parallel proof runner for doob
│   └── find_eps.sh          # find smallest eps that proves all nets at given v
└── src/
    ├── horou_c.c            # fast Newton solver (C)
    ├── horoz_c.c            # fast solver + BFS placement (C)
    └── proof_c.c            # fast prover (C)
```

---

## Build

```bash
make          # builds src/horou_c, src/horoz_c, src/proof_c
make clean
```

Requires: a C compiler and libm (standard on macOS and Linux).

---

## Tools

### horou — horoball weights

Solve for weights `u[v]` given a triangulation.

**C (fast):**
```bash
python3 ../clers/python/clers_decode.py < ../clers/prime/20.txt | src/horou_c > horou/20.bin
```
Output: one record of `NV` float64 values per net.
- `record[0]` = NaN (vertex 1 = point at infinity)
- `record[v-1]` = 1.0 for boundary vertices, solved weight for interior

**Python (readable):**
```bash
echo "CCAE" | python3 ../clers/python/clers_decode.py | python3 python/horou.py
```

---

### horoz — weights + flat positions

Solve for `u[v]` and place each vertex at `(x[v], y[v])` in the upper
half-plane via BFS over directed edges.

**C (fast):**
```bash
python3 ../clers/python/clers_decode.py < ../clers/prime/20.txt | src/horoz_c > horoz/20.bin
```
Output: one record of `3*NV` float64 values per net.
- `record[3*(v-1)+0]` = u[v]
- `record[3*(v-1)+1]` = x[v]
- `record[3*(v-1)+2]` = y[v]

Vertex 1 (infinity) has `u = NaN`, `x = NaN`, `y = NaN`.
Vertices v1, v2 (first two boundary neighbors) are pinned at (0,0) and (1,0).

**Python:**
```python
from horou import horou
from horoz import horoz
u   = horou(poly)
pos = horoz(poly, u)   # dict: vertex -> (x, y)
```

---

### proof — sub/supersolution existence proof

Proves that an ideal horoball packing exists by bracketing the solution
between a sub-solution (angle sums > 2π) and a super-solution (angle sums < 2π).

**The five checks** (all must pass with positive slack):

1. **mono** — super-solution weight > sub-solution weight at every interior vertex
2. **excess** — sub-solution angle sums exceed 2π; super-solution sums fall below 2π
3. **triangle** — triangle inequality holds for all weights in the bracket
4. **convex** — Delaunay/butterfly condition for every interior edge
5. **boundary** — boundary horoball subtended angles are less than π

A bracket `[umin, umax]` with `umin = horou(defect=-eps)` and
`umax = horou(defect=+eps)` is tried starting at `eps = 1/500`, halving
until all five checks pass or `eps` falls below `1/500000`.

**Python (readable proof):**
```bash
# Prove a single net
echo "CCAE" | python3 ../clers/python/clers_decode.py | python3 python/proof.py --verbose

# Prove all v=20 nets
python3 ../clers/python/clers_decode.py < ../clers/prime/20.txt | python3 python/proof.py
```

Sample output (verbose):
```
    mono:     0.00372
    excess:   0.00199
    triangle: 0.0814
    convex:   0.284
    boundary: 0.0403
PROVED   slack=0.0020  eps=1/500
```

**C (fast prover):**
```bash
python3 ../clers/python/clers_decode.py < ../clers/prime/20.txt | src/proof_c > proof_20.bin
```
Output: one float64 per net — `slack/eps > 0` if proved, `NaN` if solver failed.

---

## Results

All prime 6-nets through **v = 80** (73,920,746 nets at v = 80) are proved.
The eps needed stabilizes at **1/8000** from v = 60 onward.

```
 v    eps_needed      nets
 4    1/500              1
 6    1/500              1
...
25    1/1000          4711
...
44    1/4000        483170
...
60    1/8000       6,464,974
80    1/8000      73,920,746
```

See `data/eps_needed.txt` for the full table.

---

## Parallel proof on a server (doob)

```bash
# Prove a single vertex count
./scripts/run_proof_doob.sh 81

# Prove a range
./scripts/run_proof_doob.sh 81 100

# Find the smallest eps needed for a given v
./scripts/find_eps.sh 81
```

`run_proof_doob.sh` splits the input into 80 shards, runs them in parallel with
`nice -n 19`, and concatenates results into `proof_N.bin`.  Requires GNU
`parallel`.

---

## Dependencies

- **C tools:** standard C compiler, libm
- **Python tools:** Python 3, numpy (`horou.py`, `horoz.py`)
- **Parallel scripts:** GNU parallel, `neo/clers/python/clers_decode.py`
- **Prime net data:** `neo/clers/prime/N.txt` (bundled for v = 4..25;
  larger files generated by `neo/clers/scripts/grow.sh`)
