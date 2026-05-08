# C port: all-bends + quat + full-jump α-march α-carry

Implementation spec for the C version of the experiment currently
implemented in Python. Reference Python:

- `python/quat_residual.py` (residual + analytic Jacobian)
- `python/experiment_lm/lm_alpha_march_carry.py` (controller)
- `python/puffup.py::solver_lm` (LM inner)

Existing C scaffolding to reuse: `src/puffup_c_lm.c` (topology parser,
horou solver, vertexwish, dense+sparse LM with SuperLU, dent gate).

## Files

- **new** `src/lm_march_c.c` — standalone C binary, copies topology +
  horou + vertexwish from `puffup_c_lm.c`, replaces residual /
  Jacobian / controller.
- `Makefile` — add `lm-march-sparse` target (`-DHAVE_SUPERLU`).

## Python → C function map

| Python | C in `lm_march_c.c` |
|---|---|
| `q_mul`, `q_step`, `q_step_dbeta` | `q_mul`, `q_step`, `q_step_dbeta` (static inline) |
| `vertex_holonomy_quat` | `vertex_holonomy_quat` |
| `holonomy_residual_quat(vertices=tri.vertices)` | `holonomy_residual_quat_all` (3V residual rows) |
| `analytical_jacobian_quat_sparse(vertices=tri.vertices)` | `quat_AtA_fill` (per-vertex outer product into CSC `A=JᵀJ`) |
| `solver_lm` (acceptance: strict residual decrease ∧ dent gate) | `solve_lm_inner` |
| `lm_alpha_march_carry.py` controller (`--alpha-jump-policy full`) | `march_full_jump` |
| `vertex_turn` | `vertex_turn` |

## Conventions (frozen, must match Python exactly)

Quaternion `(w, x, y, z)`, Hamilton product, real-first.

Step quaternion `q_step(α, β) = q_z(α) · q_x(−β)`:

    ca = cos(α/2),  sa = sin(α/2)
    cb = cos(β/2),  sb = sin(β/2)
    q_step = ( ca·cb, −ca·sb, −sa·sb, sa·cb )

`∂q_step/∂β`:

    (−½ ca·sb, −½ ca·cb, −½ sa·cb, −½ sa·sb)

Flower order: `face_with_v_first(face, v)` returns `(v, b, c)`; step
edge is `edge_key(v, c)`. Walked via `vertex_flower[v]` from
`Tri.from_faces`.

Vertex holonomy: `Q_v = ∏_{t=0..k−1} q_step(α, β_{e_t})`.

Residual at vertex v: vector part `(Q_v[1], Q_v[2], Q_v[3])`. At a
true solution `Q_v ≈ −identity`. **Do not silently flip qw sign.**

Jacobian column for var-edge e at vertex v's flower position t:

    ∂Q_v/∂β_e = P[t] · (∂q_step/∂β)|_{α, β_{e_t}} · S[t+1]
    P[t] = ∏_{s<t} q_s,   S[t] = ∏_{s≥t} q_s,   P[0]=S[k]=identity

Three Jacobian entries (vec of `∂Q_v/∂β_e`) at rows `3·vi+(0,1,2)`,
column `var_idx[e]`.

## All-bends system

- `var_edges = ALL NE edges`. `NVAR = NE = 3V − 6`. No base freeze.
- `residual_vertices = ALL V vertices`. `NRES = 3V`.
- Over-determined by 6.
- Dent gate: `vertex_turn(v) ≥ 0` for **all V vertices**.
- Initial bends from horou (α=0 ideal).

Square mode (`--system square`) kept available as optional fallback,
unchanged from `puffup_c_lm.c`'s layout.

## Full-jump α-policy

    alpha_final = target · (1 − 1e−12)
    cap_fraction = 1.0
    try_fraction = cap_fraction
    α_curr = 0; n_step = n_retreat = 0
    history = []
    accepted_seq, rejected_seq, alpha_seq = []
    target_reached = false
    for attempt in range(max_attempts):       # exactly max_attempts iterations
        if α_curr ≥ alpha_final:
            target_reached = true; break
        gap = target − α_curr
        used = gap · try_fraction
        α_try = α_curr + used
        initial_resid = ‖F(x; α_try)‖              (stale-bend residual at new α)
        run quick LM at α_try, max_iter=max_quick_iter, tol=quick_tol, dent gate on
        resid = ‖F(out.x; α_try)‖
        accepted = (resid ≤ quick_tol) OR (resid ≤ quick_min_drop · max(initial_resid, 1e−30))
        on accept: x ← out.x; α_curr ← α_try; accepted_seq ← try_fraction;
                   alpha_seq ← degrees(α_curr); try_fraction ← cap_fraction
        on reject: rejected_seq ← try_fraction; n_retreat++;
                   try_fraction ← try_fraction · 0.5
                   if try_fraction · gap < min_step → break (STEP_UNDERFLOW)
    # Post-loop catch (lm_alpha_march_carry.py:421-422): if the for-loop
    # exits via max_attempts on the same iteration that reached the
    # threshold, the top-of-loop guard never fires but α_curr is at target.
    if α_curr ≥ alpha_final:
        target_reached = true

    if target_reached:
        final tight LM at α = target with tol = final_tol, max_iter = 200

Defaults match Python: `max_quick_iter = 10`, `quick_tol = 1e−3`,
`final_tol = 1e−12`, `lambda_init = 1.0`, `quick_min_drop = 1e−2`,
`init_step_deg = 1.0` (only used by half policy), `min_step_deg = 1e−4`,
`max_attempts = 2000`.

LM inner (`solve_lm_inner`): identical algorithm to existing
`solve_lm` in `puffup_c_lm.c` lines 1024-1171 (Marquardt
`(A + λD)δ = −g`, accept iff `‖r_trial‖ < ‖r‖` ∧ `!dent`,
`λ_down=0.3`, `λ_up=10`, `λ_max=1e12`, `tiny_rel=1e−12`,
`max_lm_retries=20`). The only differences for the all-bends quat
path are:

- residual evaluator → `holonomy_residual_quat_all`
- Jacobian-system fill → `quat_AtA_fill` (sparse) builds `A=JᵀJ` and
  `g=Jᵀr` directly from per-vertex outer products of `vec(P[t]·dq·S[t+1])`,
  symmetric to the matrix path's `sparse_lm_AtA_fill`.

`A`'s sparsity pattern is structurally identical to the matrix
path: per-vertex contribution is the outer product of the column
list `{var_idx[e_t] : t=0..k−1}`. The CSC pattern build code in
`puffup_c_lm.c::sparse_lm_setup` (lines 715–814) ports with two
edits for the all-bends path:

- replace the loop driver `INT_VS[0..N_INT)` with `1..NV` (i.e.
  iterate over **all V vertices** for both pattern build and
  per-iter fill — `NRES = 3V`, not 3·N_INT);
- `var_idx[e]` is identity on the full edge list `0..NE−1`
  (`VAR_OF_E[e] = e`), since no base-edge freeze.

Both edits are gated by a `CFG_SYSTEM == ALLBENDS` switch so the
square path stays available behind `--system square`.

## Sparse `JᵀJ` assembly (per iter, `quat_AtA_fill`)

For each vertex v with flower length k:

    Compute qs[t]    = q_step(α, β_{e_t}),  t = 0..k−1
    Compute dqs[t]   = q_step_dbeta(α, β_{e_t})
    Prefix:  P[0] = id; P[t+1] = q_mul(P[t], qs[t])
    Suffix:  S[k] = id; S[t]   = q_mul(qs[t], S[t+1])
    For each t in 0..k−1:
        dQ_t = q_mul(q_mul(P[t], dqs[t]), S[t+1])
        v_t  = (dQ_t[1], dQ_t[2], dQ_t[3])      ← vec part, 3-vector
        col  = var_idx[e_t]   (always ≥0 in all-bends mode)
        g[col] += v_t · r[3·vi : 3·vi+3]
    For each (s, t) in 0..k−1 × 0..k−1:
        cs = var_idx[e_s], ct = var_idx[e_t]
        A[cs, ct] += v_s[0]·v_t[0] + v_s[1]·v_t[1] + v_s[2]·v_t[2]

The (i, s, t) → CSC slot table is built once in `quat_AtA_setup`
exactly as in the matrix path.

## Chunk / list mode

Single-process binary supports both:

- `lm_march_c CLERS_FACELIST` (single case on stdin)
- `lm_march_c --input-list FILE --output-tsv OUT [...]` (loop over
  many CLERS-already-decoded netcodes; one TSV row per)

Format of `--input-list` file: one decoded netcode per line
(`a,b,c;a,b,d;...`). The driver script does the CLERS→facelist
decode once with `bin/clers decode` outside the C process.

Topology globals are reset between cases via `clear_topology()` (zero
all `EDGE_IDX_LOOKUP`, `DIRECTED_FACE`, `FLOWER_LEN`, etc., and free
the cached sparse setup). The `quat_AtA_setup` is rebuilt per case
(NRES/NVAR vary).

## Output TSV columns (one row per CLERS)

    v  clers  netcode_md5  status  final_alpha_deg  final_resid  attempts  retreats
    accepted_fractions  rejected_fractions  alpha_seq_deg  lm_iters_total
    lambda_final  wall_secs  failure_reason

Statuses (top three match Python; bottom three are C-sweep additions):

- `SOLVER_TOL` — target reached, `final_resid ≤ final_tol`
- `LAMBDA_SATURATED_POSITIVE_RESIDUAL` — target reached but final LM
  saturated λ
- `STALLED_POSITIVE_RESIDUAL` — α-march did not reach target, OR
  final LM stalled
- `HOROU_FAIL` — horou solver did not converge (C sweep only)
- `TOPOLOGY_FAIL` — netcode parse / closed-surface check failed
  (C sweep only)
- `EXCEPTION` — unexpected error (catch-all)

`failure_reason` distinguishes how a non-`SOLVER_TOL` case ended:

- `STEP_UNDERFLOW` — `try_fraction · gap < min_step` triggered break
- `MAX_ATTEMPTS` — exhausted `max_attempts` iterations without reach
- `LM_LAMBDA` — final tight LM saturated λ
- `LM_RETRIES` — final LM exhausted retries without accept
- `LM_STALLED` — final LM tiny-step-stall before tol
- `LM_MAX_ITER` — final LM hit `max_iter` without tol or saturation
- `HOROU` / `TOPOLOGY` / `EXCEPTION` — early-reject reasons
- empty `-` for `SOLVER_TOL` cases

## Strongest validation invariant

**On CCAE** (V=4, octahedron's smallest 6-net) with same flags
(`--system all-bends --residual quat --alpha-jump-policy full
--target-deg 60`):

- `final_resid ≤ 1e−12`
- accepted-fraction sequence in C **identical** to Python's, to ≥
  1e−10 absolute on each fraction
- final tight LM runs at exactly α = `target` (post-loop entry into
  the final-LM block; α_seq_deg’s last accepted entry is what the
  march reached, not necessarily 60.0)
- post-check: every vertex has `|vec(Q)| < 1e−6` and `qw < 0`

**On a v=10 known-good prime** (first line of `~/Dropbox/neo/data/primes/10.txt`):

- C reaches `SOLVER_TOL`
- accepted-fraction sequence in C identical to Python's

**Analytic-vs-FD spot check** at α = 30° on CCAE and on the chosen
v=10 net:

- C analytic Jacobian column j vs forward-difference at h=1e−7
  agrees to ≤ 1e−6 relative on every nonzero entry. (Built-in
  `--fd-spot-check`.) Done on a v=10 case, not just CCAE, because
  CCAE’s tetrahedral symmetry can hide indexing bugs that a generic
  v=10 case exposes.

**Dent-gate diagnostic** (`--dent-gate-trace`):

- For each attempted α step that triggers a dent rejection, the
  binary prints `dent_reject vertex=V flower_turn=…`. Required so we
  can verify the gate considers ALL V vertices in all-bends mode (not
  just the N_INT non-base set). Run on CCAE and on the v=10 case;
  cross-check that any dent rejection in the C log appears at the
  same step in the Python log.

## Failure modes the invariant is meant to catch

- sign flip anywhere in `q_step`, `q_x`, or Hamilton product
- flower walk direction reversed (CW vs CCW)
- prefix/suffix indexing off-by-one (P or S shifted)
- column indexing wrong in all-bends mode (square layout leaking
  through)
- α-policy retreat asymmetry with Python (`try_fraction *= 0.5` vs
  `step *= 0.5`)
- dent gate evaluating wrong vertex set in all-bends mode (must be
  ALL V, not interior-only)
- chunk mode leaking topology state between cases

## Smoke-test commands (in order)

    cd ~/Dropbox/neo/ideal
    make lm-march-sparse

    # 1. CCAE single case, C, with dent-gate trace
    bin/clers decode <<<CCAE | src/lm_march_c \
        --system all-bends --residual quat --alpha-jump-policy full \
        --target-deg 60 --bends-out /tmp/ccae_c.bends --dent-gate-trace \
        > /tmp/ccae_c.log 2>&1

    # 2. CCAE single case, Python (reference)
    python3 python/experiment_lm/lm_alpha_march_carry.py CCAE \
        --residual quat --system all-bends --alpha-jump-policy full \
        --target-deg 60 --bends-out /tmp/ccae_py.bends \
        > /tmp/ccae_py.log 2>&1

    # 3. Diff accepted-fraction sequences
    diff <(grep '^stats: accepted_fraction_seq' /tmp/ccae_c.log) \
         <(grep '^stats: accepted_fraction_seq' /tmp/ccae_py.log)

    # 4. Final residual: both ≤ 1e−12
    grep -E 'final_resid|status' /tmp/ccae_{c,py}.log

    # 5. v=10 known-good (first line of primes/10.txt), C with trace
    line=$(head -n1 ~/Dropbox/neo/data/primes/10.txt)
    bin/clers decode <<<"$line" | src/lm_march_c \
        --system all-bends --residual quat --alpha-jump-policy full \
        --target-deg 60 --dent-gate-trace > /tmp/v10_c.log 2>&1
    python3 python/experiment_lm/lm_alpha_march_carry.py "$line" \
        --residual quat --system all-bends --alpha-jump-policy full \
        --target-deg 60 > /tmp/v10_py.log 2>&1
    diff <(grep '^stats: accepted_fraction_seq' /tmp/v10_c.log) \
         <(grep '^stats: accepted_fraction_seq' /tmp/v10_py.log)

    # 6. FD spot checks (CCAE and v=10), with all-bends + quat flags
    src/lm_march_c --fd-spot-check --alpha-deg 30 \
        --system all-bends --residual quat < <(bin/clers decode <<<CCAE)
    src/lm_march_c --fd-spot-check --alpha-deg 30 \
        --system all-bends --residual quat < <(bin/clers decode <<<"$line")

    # 7. Chunk mode
    {
      bin/clers decode <<<CCAE
      bin/clers decode <<<CCCACAE   # v=6
      bin/clers decode <<<CCCACACCAACAAE   # v=8
    } > /tmp/chunk3.txt
    src/lm_march_c --input-list /tmp/chunk3.txt --output-tsv /tmp/chunk3.tsv \
        --system all-bends --residual quat --alpha-jump-policy full --target-deg 60
    wc -l /tmp/chunk3.tsv   # expect 4 (1 header + 3 data rows)

Only after all seven smoke tests pass: launch v=4..30 sweep on doob.

## Resolved by Codex audit

- Base-face residual in all-bends mode: every vertex contributes,
  including the three base-face vertices. Over-determination by 6
  handles the gauge. Confirmed, no suppression.
- Controller post-loop catch (`if alpha_curr >= alpha_final`) is
  required after the for-loop ends.
- Loop bound is exactly `max_attempts` iterations.
- Dent-gate fixture added (`--dent-gate-trace`) to reach the
  failure modes the CCAE-only invariant doesn't catch.
