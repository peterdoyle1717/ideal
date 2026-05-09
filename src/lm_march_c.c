/*
 * lm_march_c.c — all-bends + quaternion-residual + full-jump α-march α-carry
 *               LM solver in C, with chunk/list mode for sweeps.
 *
 * Spec: docs/all_bends_quat_c_port.md (Codex-audited).
 * Python reference:
 *   python/quat_residual.py
 *   python/experiment_lm/lm_alpha_march_carry.py
 *   python/puffup.py::solver_lm
 *
 * Build:
 *   make src/lm_march_c           (dense, no SuperLU)
 *   make src/lm_march_c_sparse    (-DHAVE_SUPERLU + libsuperlu)
 *
 * CLI (single case):
 *   lm_march_c < netcode_facelist [flags]
 *
 * CLI (chunk mode, one C process / many CLERS):
 *   lm_march_c --input-list FILE --output-tsv OUT [flags]
 *   FILE: one decoded netcode per line ("a,b,c;a,b,d;..."); preceded
 *         by an optional "# clers=..." comment line whose value is
 *         echoed into the TSV's "clers" column. If absent, "clers"
 *         defaults to "-".
 *
 * Conventions (frozen by the spec):
 *   q_step(α, β) = q_z(α) · q_x(−β)
 *                = (cos(α/2)cos(β/2),
 *                   −cos(α/2)sin(β/2),
 *                   −sin(α/2)sin(β/2),
 *                    sin(α/2)cos(β/2))
 *   ∂q_step/∂β  = (−½ ca·sb, −½ ca·cb, −½ sa·cb, −½ sa·sb)
 *   Q_v = ∏_{t=0..k−1} q_step(α, β_{e_t})  (flower order matches Python)
 *   residual at v = vec(Q_v) = (Q_v.x, Q_v.y, Q_v.z)
 *   ∂Q_v/∂β_e = P[t] · (∂q_step/∂β) · S[t+1]
 */

#define _POSIX_C_SOURCE 200809L
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef HAVE_SUPERLU
#include <slu_ddefs.h>
#endif

#ifndef PI
#define PI 3.14159265358979323846
#endif

/* ---------- size limits ---------------------------------------------------- */
#define MAXV       200
#define MAXE       (3 * MAXV)
#define MAXF       (2 * MAXV)
#define MAXFLOWER  20

/* ---------- topology globals ---------------------------------------------- */
static int NV, NE, NF;
static int FACES[MAXF][3];
static int EDGE_A[MAXE], EDGE_B[MAXE];
static int EDGE_IDX_LOOKUP[MAXV + 1][MAXV + 1];
static int VERT_DEG[MAXV + 1];
static int FLOWER_LEN[MAXV + 1];
static int FLOWER_THIRD[MAXV + 1][MAXFLOWER];
static int FLOWER_E[MAXV + 1][MAXFLOWER];
static int DIRECTED_FACE[MAXV + 1][MAXV + 1];

static int BASE[3];
static bool IS_BASE_E[MAXE];
static int NVAR;
static int VAR_OF_E[MAXE];
static int N_INT;
static int INT_VS[MAXV];

/* All-bends residual + gate vertex sets. In all-bends mode RES_VS and
 * GATE_VS span all V vertices; in square mode they cover N_INT non-base
 * vertices. */
static int N_RES;
static int RES_VS[MAXV + 1];
static int N_GATE;
static int GATE_VS[MAXV + 1];

enum { SYSTEM_ALLBENDS = 0, SYSTEM_SQUARE = 1 };
static int CFG_SYSTEM = SYSTEM_ALLBENDS;

enum { RESIDUAL_QUAT = 0, RESIDUAL_MATRIX = 1 };
static int CFG_RESIDUAL = RESIDUAL_QUAT;

enum { SOLVER_DENSE = 0, SOLVER_SPARSE = 1 };
#ifdef HAVE_SUPERLU
static int CFG_SOLVER = SOLVER_SPARSE;
#else
static int CFG_SOLVER = SOLVER_DENSE;
#endif

static int CFG_DENT_TRACE = 0;

static double tau_const;

static void die(const char *msg) {
    fprintf(stderr, "error: %s\n", msg);
    exit(2);
}

static void *xcalloc(size_t n, size_t sz) {
    void *p = calloc(n, sz);
    if (!p) die("calloc failed");
    return p;
}

static double now_secs(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* ---------- netcode parser ------------------------------------------------- */
/* Reads facelist from a string. Tokens separated by ';' or whitespace; each
 * face is "a,b,c". Returns 0 on success, -1 on parse error (so chunk-mode
 * can skip the row instead of dying on the whole sweep). */
static int parse_netcode_str(const char *s) {
    NF = 0;
    int a = 0, b = 0, c = 0, slot = 0, val = 0, have = 0;
    int max_v = 0;
    for (const char *p = s; ; p++) {
        int ch = (unsigned char)*p;
        int eos = (ch == '\0');
        if (!eos && isdigit(ch)) {
            val = val * 10 + (ch - '0');
            have = 1;
        } else if (!eos && ch == ',') {
            if (!have) return -1;
            if      (slot == 0) a = val;
            else if (slot == 1) b = val;
            else                return -1;
            slot++; val = 0; have = 0;
        } else if (eos || ch == ';' || ch == '\n' || ch == '\r' ||
                   ch == ' ' || ch == '\t') {
            if (have) {
                if (slot != 2) return -1;
                c = val;
                if (NF >= MAXF) return -1;
                if (a < 1 || b < 1 || c < 1) return -1;
                FACES[NF][0] = a; FACES[NF][1] = b; FACES[NF][2] = c;
                if (a > max_v) max_v = a;
                if (b > max_v) max_v = b;
                if (c > max_v) max_v = c;
                NF++;
            }
            slot = 0; val = 0; have = 0;
            if (eos) break;
        } else if (!eos) {
            return -1;
        }
    }
    if (NF == 0) return -1;
    NV = max_v;
    if (NV > MAXV) return -1;
    return 0;
}

/* ---------- topology build (returns 0 ok, -1 on bad netcode) -------------- */
static int add_edge_or_lookup(int u, int v) {
    int a = u < v ? u : v;
    int b = u < v ? v : u;
    int idx = EDGE_IDX_LOOKUP[a][b];
    if (idx >= 0) return idx;
    if (NE >= MAXE) return -1;
    EDGE_A[NE] = a; EDGE_B[NE] = b;
    EDGE_IDX_LOOKUP[a][b] = NE;
    return NE++;
}

static int build_topology(void) {
    NE = 0;
    for (int v = 0; v <= NV; v++)
        for (int u = 0; u <= NV; u++) {
            EDGE_IDX_LOOKUP[v][u] = -1;
            DIRECTED_FACE[v][u] = -1;
        }

    for (int fi = 0; fi < NF; fi++) {
        int a = FACES[fi][0], b = FACES[fi][1], c = FACES[fi][2];
        if (a == b || b == c || a == c) return -1;
        if (DIRECTED_FACE[a][b] >= 0) return -1;
        if (DIRECTED_FACE[b][c] >= 0) return -1;
        if (DIRECTED_FACE[c][a] >= 0) return -1;
        DIRECTED_FACE[a][b] = fi;
        DIRECTED_FACE[b][c] = fi;
        DIRECTED_FACE[c][a] = fi;
        if (add_edge_or_lookup(a, b) < 0) return -1;
        if (add_edge_or_lookup(b, c) < 0) return -1;
        if (add_edge_or_lookup(c, a) < 0) return -1;
    }
    for (int e = 0; e < NE; e++) {
        int a = EDGE_A[e], b = EDGE_B[e];
        if (DIRECTED_FACE[a][b] < 0 || DIRECTED_FACE[b][a] < 0) return -1;
    }
    for (int v = 0; v <= NV; v++) VERT_DEG[v] = 0;
    for (int e = 0; e < NE; e++) {
        VERT_DEG[EDGE_A[e]]++;
        VERT_DEG[EDGE_B[e]]++;
    }
    for (int v = 1; v <= NV; v++) {
        int start = -1;
        for (int fi = 0; fi < NF; fi++) {
            if (FACES[fi][0] == v || FACES[fi][1] == v || FACES[fi][2] == v) {
                start = fi; break;
            }
        }
        if (start < 0) return -1;
        FLOWER_LEN[v] = 0;
        int cur = start;
        int seen_count = 0;
        while (1) {
            if (seen_count >= MAXFLOWER) return -1;
            FLOWER_LEN[v]++;
            int third;
            if      (FACES[cur][0] == v) third = FACES[cur][2];
            else if (FACES[cur][1] == v) third = FACES[cur][0];
            else                          third = FACES[cur][1];
            FLOWER_THIRD[v][seen_count] = third;
            FLOWER_E[v][seen_count] =
                EDGE_IDX_LOOKUP[v < third ? v : third][v < third ? third : v];
            seen_count++;
            int nxt = DIRECTED_FACE[v][third];
            if (nxt < 0) return -1;
            if (nxt == start) break;
            cur = nxt;
        }
        if (FLOWER_LEN[v] != VERT_DEG[v]) return -1;
    }
    return 0;
}

/* ---------- system layout: all-bends or square ---------------------------- */
static void choose_base_face(void) {
    BASE[0] = FACES[0][0];
    BASE[1] = FACES[0][1];
    BASE[2] = FACES[0][2];
    for (int e = 0; e < NE; e++) IS_BASE_E[e] = false;
}

static void setup_all_bends(void) {
    /* var_edges = ALL NE edges (no base freeze); residuals at all V verts. */
    NVAR = NE;
    for (int e = 0; e < NE; e++) VAR_OF_E[e] = e;
    N_RES = NV;
    for (int i = 0; i < NV; i++) RES_VS[i] = i + 1;
    N_GATE = NV;
    for (int i = 0; i < NV; i++) GATE_VS[i] = i + 1;
    /* INT_VS / N_INT unused in this mode; leave consistent for safety. */
    N_INT = 0;
}

static void setup_square(void) {
    /* var_edges = NE − 3 non-base edges; residuals at N_INT non-base verts. */
    int e01_a = BASE[0] < BASE[1] ? BASE[0] : BASE[1];
    int e01_b = BASE[0] < BASE[1] ? BASE[1] : BASE[0];
    int e12_a = BASE[1] < BASE[2] ? BASE[1] : BASE[2];
    int e12_b = BASE[1] < BASE[2] ? BASE[2] : BASE[1];
    int e20_a = BASE[2] < BASE[0] ? BASE[2] : BASE[0];
    int e20_b = BASE[2] < BASE[0] ? BASE[0] : BASE[2];
    int e01 = EDGE_IDX_LOOKUP[e01_a][e01_b];
    int e12 = EDGE_IDX_LOOKUP[e12_a][e12_b];
    int e20 = EDGE_IDX_LOOKUP[e20_a][e20_b];
    IS_BASE_E[e01] = IS_BASE_E[e12] = IS_BASE_E[e20] = true;
    NVAR = 0;
    for (int e = 0; e < NE; e++) VAR_OF_E[e] = IS_BASE_E[e] ? -1 : NVAR++;
    N_INT = 0;
    for (int v = 1; v <= NV; v++) {
        if (v != BASE[0] && v != BASE[1] && v != BASE[2]) {
            INT_VS[N_INT++] = v;
        }
    }
    N_RES = N_INT;
    for (int i = 0; i < N_INT; i++) RES_VS[i] = INT_VS[i];
    N_GATE = N_INT;
    for (int i = 0; i < N_INT; i++) GATE_VS[i] = INT_VS[i];
}

/* ---------- horou solver (port of horou_c.c, lifted from puffup_c_lm.c) -- */
static double horou_petal(double ui, double uj, double uk) {
    double a=ui*uj, b=ui*uk, c=uj*uk;
    double cos_t = (a*a + b*b - c*c) / (2.0*a*b);
    if (cos_t >  1.0) cos_t =  1.0;
    if (cos_t < -1.0) cos_t = -1.0;
    return acos(cos_t);
}
static void horou_petal_grad(double ui, double uj, double uk,
                              double *dui, double *duj, double *duk) {
    double a=ui*uj, b=ui*uk, c=uj*uk;
    double cos_t = (a*a + b*b - c*c) / (2.0*a*b);
    if (cos_t >  1.0) cos_t =  1.0;
    if (cos_t < -1.0) cos_t = -1.0;
    double sin_t = sqrt(1.0 - cos_t*cos_t);
    if (sin_t < 1e-15) { *dui = *duj = *duk = 0.0; return; }
    double s = -1.0 / sin_t;
    *dui = s * (uj*uk / (ui*ui*ui));
    *duj = s * (1.0/(2*uk) - uk/(2*uj*uj) - uk/(2*ui*ui));
    *duk = s * (1.0/(2*uj) - uj/(2*uk*uk) - uj/(2*ui*ui));
}

static int horou_lu_solve_dense(double *J, double *b, int n) {
    for (int col = 0; col < n; col++) {
        int piv = col;
        double best = fabs(J[col*n + col]);
        for (int r = col+1; r < n; r++) {
            double v = fabs(J[r*n + col]);
            if (v > best) { best = v; piv = r; }
        }
        if (best < 1e-14) return -1;
        if (piv != col) {
            for (int k = col; k < n; k++) {
                double t = J[col*n + k]; J[col*n + k] = J[piv*n + k]; J[piv*n + k] = t;
            }
            double t = b[col]; b[col] = b[piv]; b[piv] = t;
        }
        double inv = 1.0 / J[col*n + col];
        for (int r = col+1; r < n; r++) {
            double fac = J[r*n + col] * inv;
            for (int k = col; k < n; k++) J[r*n + k] -= fac * J[col*n + k];
            b[r] -= fac * b[col];
        }
    }
    for (int i = n-1; i >= 0; i--) {
        double s = b[i];
        for (int j = i+1; j < n; j++) s -= J[i*n + j] * b[j];
        b[i] = s / J[i*n + i];
    }
    return 0;
}

static int horou_solve(double *u_out) {
    const double TAU = 2.0 * PI;
    int *bndry = (int *)xcalloc(NV + 2, sizeof(int));
    int *int_idx = (int *)xcalloc(NV + 2, sizeof(int));
    for (int v = 0; v <= NV + 1; v++) int_idx[v] = -1;

    int v1_deg = FLOWER_LEN[1];
    for (int t = 0; t < v1_deg; t++) {
        int nb = FLOWER_THIRD[1][t];
        bndry[nb] = 1;
    }
    int n_int = 0;
    int *interior = (int *)xcalloc(NV + 2, sizeof(int));
    for (int v = 2; v <= NV; v++) {
        if (!bndry[v]) { int_idx[v] = n_int; interior[n_int++] = v; }
    }
    int *ringlen = (int *)xcalloc(NV + 2, sizeof(int));
    for (int i = 0; i < n_int; i++) ringlen[interior[i]] = FLOWER_LEN[interior[i]];

    double *xvec = (double *)xcalloc(NV + 2, sizeof(double));
    for (int i = 0; i < n_int; i++) xvec[i] = 1.0;

    int nff = 0;
    int *ff_a = (int *)xcalloc(NF, sizeof(int));
    int *ff_b = (int *)xcalloc(NF, sizeof(int));
    int *ff_c = (int *)xcalloc(NF, sizeof(int));
    for (int i = 0; i < NF; i++) {
        int a = FACES[i][0], b = FACES[i][1], c = FACES[i][2];
        if (a != 1 && b != 1 && c != 1) {
            ff_a[nff] = a; ff_b[nff] = b; ff_c[nff] = c; nff++;
        }
    }

    if (n_int == 0) goto done;

    double *HJ = (double *)xcalloc((size_t)n_int * (size_t)n_int, sizeof(double));
    double *HF = (double *)xcalloc(n_int, sizeof(double));
    double *Hdx = (double *)xcalloc(n_int, sizeof(double));

    #define U_AT(vv) (bndry[vv] ? 1.0 : xvec[int_idx[vv]])

    for (int it = 0; it < 200; it++) {
        double res = 0.0;
        for (int i = 0; i < n_int; i++) {
            int v = interior[i]; int k = ringlen[v]; double ui = xvec[i];
            double s = 0.0;
            for (int j = 0; j < k; j++) {
                int rj = FLOWER_THIRD[v][j];
                int rk = FLOWER_THIRD[v][(j+1) % k];
                s += horou_petal(ui, U_AT(rj), U_AT(rk));
            }
            HF[i] = s - TAU;
            double af = fabs(HF[i]); if (af > res) res = af;
        }
        if (res < 1e-12) break;

        memset(HJ, 0, (size_t)n_int * (size_t)n_int * sizeof(double));
        for (int i = 0; i < n_int; i++) {
            int v = interior[i]; int k = ringlen[v]; double ui = xvec[i];
            for (int j = 0; j < k; j++) {
                int vj = FLOWER_THIRD[v][j];
                int vk = FLOWER_THIRD[v][(j+1) % k];
                double uj = U_AT(vj), uk_ = U_AT(vk);
                double dui, duj, duk;
                horou_petal_grad(ui, uj, uk_, &dui, &duj, &duk);
                HJ[i*n_int + i] += dui;
                if (int_idx[vj] >= 0) HJ[i*n_int + int_idx[vj]] += duj;
                if (int_idx[vk] >= 0) HJ[i*n_int + int_idx[vk]] += duk;
            }
        }
        for (int i = 0; i < n_int; i++) Hdx[i] = -HF[i];
        if (horou_lu_solve_dense(HJ, Hdx, n_int) < 0) {
            free(HJ); free(HF); free(Hdx);
            free(ff_a); free(ff_b); free(ff_c);
            free(xvec); free(ringlen); free(interior); free(int_idx); free(bndry);
            return -1;
        }
        double step = 1.0; int found = 0;
        for (int bt = 0; bt < 60; bt++, step *= 0.5) {
            int ok = 1;
            for (int i = 0; i < n_int; i++)
                if (xvec[i] + step*Hdx[i] <= 0.0) { ok = 0; break; }
            if (!ok) continue;
            for (int f = 0; f < nff && ok; f++) {
                #define UU(vv) (int_idx[vv] >= 0 \
                    ? xvec[int_idx[vv]] + step*Hdx[int_idx[vv]] : 1.0)
                double ua = UU(ff_a[f]), ub = UU(ff_b[f]), uc = UU(ff_c[f]);
                double p = ua*ub, q = ua*uc, r = ub*uc;
                if (p+q <= r || p+r <= q || q+r <= p) ok = 0;
                #undef UU
            }
            if (!ok) continue;
            double res2 = 0.0;
            for (int i = 0; i < n_int; i++) {
                int v = interior[i]; int k = ringlen[v];
                double ui = xvec[i] + step*Hdx[i];
                double s = 0.0;
                for (int j = 0; j < k; j++) {
                    int rj = FLOWER_THIRD[v][j];
                    int rk = FLOWER_THIRD[v][(j+1) % k];
                    double uj = (int_idx[rj] >= 0 ? xvec[int_idx[rj]] + step*Hdx[int_idx[rj]] : 1.0);
                    double uk_ = (int_idx[rk] >= 0 ? xvec[int_idx[rk]] + step*Hdx[int_idx[rk]] : 1.0);
                    s += horou_petal(ui, uj, uk_);
                }
                double af = fabs(s - TAU); if (af > res2) res2 = af;
            }
            if (res2 < res) { found = 1; break; }
        }
        if (!found) break;
        for (int i = 0; i < n_int; i++) xvec[i] += step * Hdx[i];
    }
    free(HJ); free(HF); free(Hdx);

    #undef U_AT

done:
    u_out[0] = NAN;
    for (int v = 2; v <= NV; v++) {
        u_out[v - 1] = bndry[v] ? 1.0 : xvec[int_idx[v]];
    }
    free(ff_a); free(ff_b); free(ff_c);
    free(xvec); free(ringlen); free(interior); free(int_idx); free(bndry);
    return 0;
}

/* α=0 dihedral per edge from u: bend = π − interior dihedral. Lifted from
 * puffup_c_lm.c::compute_bends_at_zero. */
static double horou_angleat(int v, int p, int q, const double u[]) {
    double x = 1.0/u[v-1], y = 1.0/u[p-1], z = 1.0/u[q-1];
    double cos_t = (y*y + z*z - x*x) / (2.0*y*z);
    if (cos_t >  1.0) cos_t =  1.0;
    if (cos_t < -1.0) cos_t = -1.0;
    return acos(cos_t);
}
static double horou_boundaryangleat(int v, const double u[]) {
    double tot = 0.0;
    int k = FLOWER_LEN[v];
    for (int t = 0; t < k; t++) {
        int third = FLOWER_THIRD[v][t];
        int fi = DIRECTED_FACE[v][third];
        int a = FACES[fi][0], b = FACES[fi][1], c = FACES[fi][2];
        if (a == 1 || b == 1 || c == 1) continue;
        int p, q;
        if      (v == a) { p = b; q = c; }
        else if (v == b) { p = c; q = a; }
        else             { p = a; q = b; }
        tot += horou_angleat(v, p, q, u);
    }
    return tot;
}
static int compute_bends_at_zero(const double u[], double bend_out[]) {
    for (int e = 0; e < NE; e++) {
        int a = EDGE_A[e], b = EDGE_B[e];
        int f1 = DIRECTED_FACE[a][b];
        int f2 = DIRECTED_FACE[b][a];
        int third1, third2;
        {
            int x = FACES[f1][0], y = FACES[f1][1], z = FACES[f1][2];
            third1 = (x != a && x != b) ? x : ((y != a && y != b) ? y : z);
        }
        {
            int x = FACES[f2][0], y = FACES[f2][1], z = FACES[f2][2];
            third2 = (x != a && x != b) ? x : ((y != a && y != b) ? y : z);
        }
        int touches_inf_1 = (FACES[f1][0] == 1 || FACES[f1][1] == 1 || FACES[f1][2] == 1);
        int touches_inf_2 = (FACES[f2][0] == 1 || FACES[f2][1] == 1 || FACES[f2][2] == 1);

        double bend;
        if (a == 1 || b == 1) {
            int v = (a == 1) ? b : a;
            bend = PI - horou_boundaryangleat(v, u);
        } else if (!touches_inf_1 && !touches_inf_2) {
            double pa = horou_angleat(third1, a, b, u);
            double pd = horou_angleat(third2, a, b, u);
            bend = PI - pa - pd;
        } else {
            int finite_third = touches_inf_1 ? third2 : third1;
            bend = PI - horou_angleat(finite_third, a, b, u);
        }
        bend_out[e] = bend;
    }
    return 0;
}

/* ---------- quaternion utilities ----------------------------------------- */
typedef struct { double w, x, y, z; } Quat;

static inline Quat q_mul(Quat a, Quat b) {
    Quat r;
    r.w = a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z;
    r.x = a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y;
    r.y = a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x;
    r.z = a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w;
    return r;
}

static inline Quat q_step(double alpha, double beta) {
    /* q_z(α) · q_x(−β) */
    double ca = cos(alpha * 0.5), sa = sin(alpha * 0.5);
    double cb = cos(beta  * 0.5), sb = sin(beta  * 0.5);
    Quat q = { ca*cb, -ca*sb, -sa*sb, sa*cb };
    return q;
}

static inline Quat q_step_dbeta(double alpha, double beta) {
    /* ∂q_step/∂β: matches python/quat_residual.py:86–94 */
    double ca = cos(alpha * 0.5), sa = sin(alpha * 0.5);
    double cb = cos(beta  * 0.5), sb = sin(beta  * 0.5);
    Quat q = { -0.5*ca*sb, -0.5*ca*cb, -0.5*sa*cb, -0.5*sa*sb };
    return q;
}

static inline Quat q_id(void) {
    Quat q = { 1.0, 0.0, 0.0, 0.0 };
    return q;
}

/* ---------- residuals ----------------------------------------------------- */
/* Stacked vector parts of per-vertex holonomy quaternions (3 entries per
 * residual vertex). r length = 3 * N_RES. */
static void holonomy_residual_quat(double alpha, const double bend[], double r[]) {
    int rk = 0;
    for (int t = 0; t < N_RES; t++) {
        int v = RES_VS[t];
        Quat Q = q_id();
        int kk = FLOWER_LEN[v];
        for (int s = 0; s < kk; s++) {
            int e = FLOWER_E[v][s];
            Q = q_mul(Q, q_step(alpha, bend[e]));
        }
        r[rk++] = Q.x;
        r[rk++] = Q.y;
        r[rk++] = Q.z;
    }
}

/* ---------- analytic Jacobian (dense; for FD spot-check) ----------------- */
/* J shape: (3*N_RES) × NVAR, row-major. Computes per-vertex prefix/suffix
 * quaternion products and accumulates contributions. */
static void quat_analytic_jacobian_dense(double alpha, const double bend[],
                                          double *J) {
    int rows = 3 * N_RES;
    memset(J, 0, (size_t)rows * (size_t)NVAR * sizeof(double));
    Quat qs[MAXFLOWER], dqs[MAXFLOWER];
    Quat P[MAXFLOWER + 1], S[MAXFLOWER + 1];
    for (int i = 0; i < N_RES; i++) {
        int v = RES_VS[i];
        int k = FLOWER_LEN[v];
        for (int t = 0; t < k; t++) {
            int e = FLOWER_E[v][t];
            qs[t]  = q_step(alpha, bend[e]);
            dqs[t] = q_step_dbeta(alpha, bend[e]);
        }
        P[0] = q_id();
        for (int t = 0; t < k; t++) P[t+1] = q_mul(P[t], qs[t]);
        S[k] = q_id();
        for (int t = k - 1; t >= 0; t--) S[t] = q_mul(qs[t], S[t+1]);
        for (int t = 0; t < k; t++) {
            int e = FLOWER_E[v][t];
            int col = VAR_OF_E[e];
            if (col < 0) continue;
            Quat dQ = q_mul(q_mul(P[t], dqs[t]), S[t+1]);
            int row = 3 * i;
            J[(row + 0) * NVAR + col] += dQ.x;
            J[(row + 1) * NVAR + col] += dQ.y;
            J[(row + 2) * NVAR + col] += dQ.z;
        }
    }
}

/* Forward-difference reference (central differences). For spot-checks. */
static void quat_fd_jacobian_dense(double alpha, const double bend[],
                                     double *J) {
    int rows = 3 * N_RES;
    memset(J, 0, (size_t)rows * (size_t)NVAR * sizeof(double));
    double *bend_tmp = (double *)xcalloc(NE, sizeof(double));
    double *r_plus = (double *)xcalloc(rows, sizeof(double));
    double *r_minus = (double *)xcalloc(rows, sizeof(double));
    memcpy(bend_tmp, bend, NE * sizeof(double));
    for (int e = 0; e < NE; e++) {
        int j = VAR_OF_E[e];
        if (j < 0) continue;
        double original = bend_tmp[e];
        double h = 1e-7 * (fabs(original) > 1.0 ? fabs(original) : 1.0);
        bend_tmp[e] = original + h;
        holonomy_residual_quat(alpha, bend_tmp, r_plus);
        bend_tmp[e] = original - h;
        holonomy_residual_quat(alpha, bend_tmp, r_minus);
        bend_tmp[e] = original;
        for (int i = 0; i < rows; i++) {
            J[i * NVAR + j] = (r_plus[i] - r_minus[i]) / (2.0 * h);
        }
    }
    free(bend_tmp); free(r_plus); free(r_minus);
}

/* ---------- sparse JᵀJ assembly ------------------------------------------ */
#ifdef HAVE_SUPERLU
static int_t  *SQ_colptr = NULL;
static int_t  *SQ_rowidx = NULL;
static double *SQ_val    = NULL;
static double *SQ_val_pl = NULL;
static int_t   SQ_NNZ    = 0;
static int_t  *SQ_diag_pos = NULL;
/* Per-(i, s, t) offset into SQ_val. i ∈ 0..N_RES, s,t ∈ 0..MAXFLOWER. */
static int_t  *SQ_off_st = NULL;
static int    *SQ_perm_c = NULL;
static int    *SQ_perm_r = NULL;
static int     SQ_setup_done = 0;

static int **g_col_rows;
static int  *g_col_cap;
static int  *g_col_len;

static void sq_insert_row(int col, int row) {
    int n = g_col_len[col];
    int lo = 0, hi = n;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (g_col_rows[col][mid] < row) lo = mid + 1; else hi = mid;
    }
    if (lo < n && g_col_rows[col][lo] == row) return;
    if (n + 1 > g_col_cap[col]) {
        int newcap = g_col_cap[col] ? g_col_cap[col] * 2 : 8;
        g_col_rows[col] = (int *)realloc(g_col_rows[col],
                                         (size_t)newcap * sizeof(int));
        if (!g_col_rows[col]) die("realloc col_rows");
        g_col_cap[col] = newcap;
    }
    for (int p = n; p > lo; p--) g_col_rows[col][p] = g_col_rows[col][p-1];
    g_col_rows[col][lo] = row;
    g_col_len[col] = n + 1;
}

static void sparse_quat_free(void) {
    free(SQ_colptr); SQ_colptr = NULL;
    free(SQ_rowidx); SQ_rowidx = NULL;
    free(SQ_val);    SQ_val = NULL;
    free(SQ_val_pl); SQ_val_pl = NULL;
    free(SQ_diag_pos); SQ_diag_pos = NULL;
    free(SQ_off_st); SQ_off_st = NULL;
    if (SQ_perm_c) { SUPERLU_FREE(SQ_perm_c); SQ_perm_c = NULL; }
    if (SQ_perm_r) { SUPERLU_FREE(SQ_perm_r); SQ_perm_r = NULL; }
    SQ_NNZ = 0;
    SQ_setup_done = 0;
}

static int sparse_quat_setup(void) {
    sparse_quat_free();

    g_col_rows = (int **)xcalloc((size_t)NVAR, sizeof(int *));
    g_col_cap  = (int  *)xcalloc((size_t)NVAR, sizeof(int));
    g_col_len  = (int  *)xcalloc((size_t)NVAR, sizeof(int));

    for (int i = 0; i < N_RES; i++) {
        int v = RES_VS[i];
        int k = FLOWER_LEN[v];
        for (int s = 0; s < k; s++) {
            int es = FLOWER_E[v][s];
            int cs = VAR_OF_E[es];
            if (cs < 0) continue;
            for (int t = 0; t < k; t++) {
                int et = FLOWER_E[v][t];
                int ct = VAR_OF_E[et];
                if (ct < 0) continue;
                sq_insert_row(ct, cs);
            }
        }
    }

    SQ_colptr = (int_t *)xcalloc((size_t)NVAR + 1, sizeof(int_t));
    int_t total = 0;
    for (int c = 0; c < NVAR; c++) {
        SQ_colptr[c] = total;
        total += g_col_len[c];
    }
    SQ_colptr[NVAR] = total;
    SQ_NNZ = total;
    SQ_rowidx = (int_t *)xcalloc((size_t)SQ_NNZ, sizeof(int_t));
    SQ_val    = (double *)xcalloc((size_t)SQ_NNZ, sizeof(double));
    SQ_val_pl = (double *)xcalloc((size_t)SQ_NNZ, sizeof(double));
    for (int c = 0; c < NVAR; c++) {
        int_t off = SQ_colptr[c];
        for (int j = 0; j < g_col_len[c]; j++) {
            SQ_rowidx[off + j] = g_col_rows[c][j];
        }
    }
    SQ_diag_pos = (int_t *)xcalloc((size_t)NVAR, sizeof(int_t));
    for (int c = 0; c < NVAR; c++) {
        SQ_diag_pos[c] = -1;
        int lo = SQ_colptr[c], hi = SQ_colptr[c+1];
        while (lo < hi) {
            int mid = (lo + hi) / 2;
            if (SQ_rowidx[mid] < c) lo = mid + 1; else hi = mid;
        }
        if (lo < SQ_colptr[c+1] && SQ_rowidx[lo] == c) SQ_diag_pos[c] = lo;
    }

    SQ_off_st = (int_t *)xcalloc((size_t)N_RES * MAXFLOWER * MAXFLOWER, sizeof(int_t));
    for (size_t z = 0; z < (size_t)N_RES * MAXFLOWER * MAXFLOWER; z++) SQ_off_st[z] = -1;
    for (int i = 0; i < N_RES; i++) {
        int v = RES_VS[i];
        int k = FLOWER_LEN[v];
        for (int s = 0; s < k; s++) {
            int cs = VAR_OF_E[FLOWER_E[v][s]];
            if (cs < 0) continue;
            for (int t = 0; t < k; t++) {
                int ct = VAR_OF_E[FLOWER_E[v][t]];
                if (ct < 0) continue;
                int lo = SQ_colptr[ct], hi = SQ_colptr[ct+1];
                while (lo < hi) {
                    int mid = (lo + hi) / 2;
                    if (SQ_rowidx[mid] < cs) lo = mid + 1; else hi = mid;
                }
                SQ_off_st[((size_t)i * MAXFLOWER + s) * MAXFLOWER + t] = lo;
            }
        }
    }

    for (int c = 0; c < NVAR; c++) free(g_col_rows[c]);
    free(g_col_rows); free(g_col_cap); free(g_col_len);
    g_col_rows = NULL; g_col_cap = NULL; g_col_len = NULL;

    /* Column ordering on the pattern. */
    for (int_t i = 0; i < SQ_NNZ; i++) SQ_val[i] = 1.0;
    SuperMatrix Apat;
    dCreate_CompCol_Matrix(&Apat, NVAR, NVAR, SQ_NNZ,
                            SQ_val, SQ_rowidx, SQ_colptr,
                            SLU_NC, SLU_D, SLU_GE);
    SQ_perm_c = intMalloc(NVAR);
    SQ_perm_r = intMalloc(NVAR);
    if (!SQ_perm_c || !SQ_perm_r) {
        Destroy_SuperMatrix_Store(&Apat);
        sparse_quat_free();
        return -1;
    }
    get_perm_c(3, &Apat, SQ_perm_c);
    Destroy_SuperMatrix_Store(&Apat);

    SQ_setup_done = 1;
    return 0;
}

/* Per-iter assembly of A = JᵀJ + g = Jᵀr from per-vertex outer products
 * of the 3-vector v_t = vec(P[t] · dqs[t] · S[t+1]). */
static void sparse_quat_AtA_fill(double alpha, const double *bend,
                                   const double *r, double *g) {
    memset(SQ_val, 0, (size_t)SQ_NNZ * sizeof(double));
    memset(g, 0, (size_t)NVAR * sizeof(double));

    Quat qs[MAXFLOWER], dqs[MAXFLOWER];
    Quat P[MAXFLOWER + 1], S[MAXFLOWER + 1];
    double v01[MAXFLOWER], v02[MAXFLOWER], v12[MAXFLOWER];
    int    cs_arr[MAXFLOWER];

    for (int i = 0; i < N_RES; i++) {
        int v = RES_VS[i];
        int k = FLOWER_LEN[v];
        for (int t = 0; t < k; t++) {
            int e = FLOWER_E[v][t];
            qs[t]  = q_step(alpha, bend[e]);
            dqs[t] = q_step_dbeta(alpha, bend[e]);
        }
        P[0] = q_id();
        for (int t = 0; t < k; t++) P[t+1] = q_mul(P[t], qs[t]);
        S[k] = q_id();
        for (int t = k - 1; t >= 0; t--) S[t] = q_mul(qs[t], S[t+1]);

        for (int s = 0; s < k; s++) {
            int e = FLOWER_E[v][s];
            int cs = VAR_OF_E[e];
            cs_arr[s] = cs;
            if (cs < 0) {
                v01[s] = v02[s] = v12[s] = 0.0;
            } else {
                Quat dQ = q_mul(q_mul(P[s], dqs[s]), S[s+1]);
                v01[s] = dQ.x;  /* row 3i+0 */
                v02[s] = dQ.y;  /* row 3i+1 */
                v12[s] = dQ.z;  /* row 3i+2 */
                g[cs] += v01[s]*r[3*i+0] + v02[s]*r[3*i+1] + v12[s]*r[3*i+2];
            }
        }
        for (int s = 0; s < k; s++) {
            if (cs_arr[s] < 0) continue;
            for (int t = 0; t < k; t++) {
                if (cs_arr[t] < 0) continue;
                int_t off = SQ_off_st[((size_t)i * MAXFLOWER + s) * MAXFLOWER + t];
                if (off < 0) continue;
                SQ_val[off] += v01[s]*v01[t] + v02[s]*v02[t] + v12[s]*v12[t];
            }
        }
    }
}

static int sparse_quat_step_solve(double lam, const double *D, double *bwork) {
    memcpy(SQ_val_pl, SQ_val, (size_t)SQ_NNZ * sizeof(double));
    for (int c = 0; c < NVAR; c++) {
        int_t dpos = SQ_diag_pos[c];
        if (dpos < 0) return -2;
        SQ_val_pl[dpos] += lam * D[c];
    }
    SuperMatrix A, B;
    SuperMatrix L = {0};
    SuperMatrix U = {0};
    dCreate_CompCol_Matrix(&A, NVAR, NVAR, SQ_NNZ,
                            SQ_val_pl, SQ_rowidx, SQ_colptr,
                            SLU_NC, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&B, NVAR, 1, bwork, NVAR, SLU_DN, SLU_D, SLU_GE);

    superlu_options_t options;
    set_default_options(&options);
    options.ColPerm = MY_PERMC;
    SuperLUStat_t stat;
    StatInit(&stat);

    int_t info = 0;
    dgssv(&options, &A, SQ_perm_c, SQ_perm_r, &L, &U, &B, &stat, &info);

    Destroy_SuperMatrix_Store(&A);
    Destroy_SuperMatrix_Store(&B);
    if (L.Store) Destroy_SuperNode_Matrix(&L);
    if (U.Store) Destroy_CompCol_Matrix(&U);
    StatFree(&stat);
    return (info == 0) ? 0 : -1;
}
#endif /* HAVE_SUPERLU */

/* ---------- dense LM linear solve ---------------------------------------- */
static int dense_solve_inplace(double *A, double *b, int n) {
    for (int k = 0; k < n; k++) {
        int piv = k;
        double best = fabs(A[k * n + k]);
        for (int i = k + 1; i < n; i++) {
            double v = fabs(A[i * n + k]);
            if (v > best) { best = v; piv = i; }
        }
        if (best < 1e-30) return 1;
        if (piv != k) {
            for (int j = k; j < n; j++) {
                double t = A[k * n + j]; A[k * n + j] = A[piv * n + j]; A[piv * n + j] = t;
            }
            double t = b[k]; b[k] = b[piv]; b[piv] = t;
        }
        double p = A[k * n + k];
        for (int i = k + 1; i < n; i++) {
            double f = A[i * n + k] / p;
            A[i * n + k] = 0.0;
            for (int j = k + 1; j < n; j++) A[i * n + j] -= f * A[k * n + j];
            b[i] -= f * b[k];
        }
    }
    double *x = (double *)xcalloc(n, sizeof(double));
    for (int i = n - 1; i >= 0; i--) {
        double s = b[i];
        for (int j = i + 1; j < n; j++) s -= A[i * n + j] * x[j];
        x[i] = s / A[i * n + i];
    }
    memcpy(b, x, n * sizeof(double));
    free(x);
    return 0;
}

/* ---------- dent gate ---------------------------------------------------- */
static double vertex_turn(int v, const double bend[]) {
    double s = 0;
    int k = FLOWER_LEN[v];
    /* match python: dedupe by edge (an edge appears twice in v's flower
     * only if v has self-edges, which our triangulations don't, but be
     * safe). FLOWER_E lists v's edge per flower step; in a simple
     * triangulation each edge appears exactly once. */
    for (int t = 0; t < k; t++) s += bend[FLOWER_E[v][t]];
    return s;
}

static int has_dent_with_first(const double bend[], int *first_v) {
    for (int t = 0; t < N_GATE; t++) {
        int v = GATE_VS[t];
        if (vertex_turn(v, bend) < 0.0) {
            if (first_v) *first_v = v;
            return 1;
        }
    }
    if (first_v) *first_v = 0;
    return 0;
}

static double vec_norm(const double *v, int n) {
    double s = 0;
    for (int i = 0; i < n; i++) s += v[i] * v[i];
    return sqrt(s);
}

/* ---------- LM inner ----------------------------------------------------- */
typedef struct {
    bool   success;
    int    iters;
    double final_resid;
    double final_lambda;
    int    total_retries_residual;
    int    total_retries_dent;
    char   msg[32];   /* "tol", "lambda_saturated", "lambda_saturated_linsolve",
                       * "lm_retries_exhausted", "stalled", "max_iter",
                       * "max_iter_tol" */
} LMOut;

static LMOut solve_lm_inner(double alpha, double bend_full[],
                              double tol, int max_iter, int max_lm_retries,
                              double lambda_init, double lambda_down, double lambda_up,
                              double lambda_max, double tiny_rel,
                              int *out_iters_used) {
    LMOut out = {0};
    snprintf(out.msg, sizeof(out.msg), "uninit");
    int M = 3 * N_RES;
    int N = NVAR;

    double *r        = (double *)xcalloc(M, sizeof(double));
    double *r_trial  = (double *)xcalloc(M, sizeof(double));
    double *bend_trial = (double *)xcalloc(NE, sizeof(double));
    double *J        = NULL;
    double *A        = NULL;
    double *Awork    = NULL;
    double *g        = (double *)xcalloc(N, sizeof(double));
    double *D        = (double *)xcalloc(N, sizeof(double));
    double *delta    = (double *)xcalloc(N, sizeof(double));
    double *bwork    = (double *)xcalloc(N, sizeof(double));
    if (CFG_SOLVER == SOLVER_DENSE) {
        J     = (double *)xcalloc((size_t)M * (size_t)N, sizeof(double));
        A     = (double *)xcalloc((size_t)N * (size_t)N, sizeof(double));
        Awork = (double *)xcalloc((size_t)N * (size_t)N, sizeof(double));
    }

    holonomy_residual_quat(alpha, bend_full, r);
    double norm = vec_norm(r, M);
    double lam = lambda_init;

    int it;
    for (it = 0; it < max_iter; it++) {
        if (norm <= tol) {
            out.success = true; out.iters = it; out.final_resid = norm;
            out.final_lambda = lam; snprintf(out.msg, sizeof(out.msg), "tol");
            goto done;
        }
        if (CFG_SOLVER == SOLVER_DENSE) {
            quat_analytic_jacobian_dense(alpha, bend_full, J);
            for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
                double s = 0;
                for (int k = 0; k < M; k++) s += J[k * N + i] * J[k * N + j];
                A[i * N + j] = s;
            }
            for (int i = 0; i < N; i++) {
                double s = 0;
                for (int k = 0; k < M; k++) s += J[k * N + i] * r[k];
                g[i] = s;
            }
            double max_diag = 0;
            for (int i = 0; i < N; i++) if (A[i*N+i] > max_diag) max_diag = A[i*N+i];
            double floor_v = tiny_rel * max_diag;
            if (floor_v < 1e-30) floor_v = 1e-30;
            for (int i = 0; i < N; i++) {
                D[i] = A[i*N+i] > floor_v ? A[i*N+i] : floor_v;
            }
        }
#ifdef HAVE_SUPERLU
        else {
            sparse_quat_AtA_fill(alpha, bend_full, r, g);
            double max_diag = 0;
            for (int c = 0; c < N; c++) {
                int_t dpos = SQ_diag_pos[c];
                double diag = (dpos >= 0) ? SQ_val[dpos] : 0.0;
                if (diag > max_diag) max_diag = diag;
            }
            double floor_v = tiny_rel * max_diag;
            if (floor_v < 1e-30) floor_v = 1e-30;
            for (int c = 0; c < N; c++) {
                int_t dpos = SQ_diag_pos[c];
                double diag = (dpos >= 0) ? SQ_val[dpos] : 0.0;
                D[c] = diag > floor_v ? diag : floor_v;
            }
        }
#endif

        bool accepted = false;
        for (int rt = 0; rt < max_lm_retries; rt++) {
            int sing;
            if (CFG_SOLVER == SOLVER_DENSE) {
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < N; j++) Awork[i * N + j] = A[i * N + j];
                    Awork[i * N + i] += lam * D[i];
                    bwork[i] = -g[i];
                }
                sing = dense_solve_inplace(Awork, bwork, N);
            }
#ifdef HAVE_SUPERLU
            else {
                for (int i = 0; i < N; i++) bwork[i] = -g[i];
                sing = sparse_quat_step_solve(lam, D, bwork);
            }
#else
            else { sing = -1; }
#endif
            if (sing) {
                lam = lam * lambda_up;
                if (lam > lambda_max) {
                    out.iters = it; out.final_resid = norm; out.final_lambda = lam;
                    snprintf(out.msg, sizeof(out.msg), "lambda_saturated_linsolve");
                    goto done;
                }
                continue;
            }
            for (int i = 0; i < N; i++) delta[i] = bwork[i];

            for (int e = 0; e < NE; e++) bend_trial[e] = bend_full[e];
            for (int e = 0; e < NE; e++) {
                int j = VAR_OF_E[e];
                if (j >= 0) bend_trial[e] += delta[j];
            }
            holonomy_residual_quat(alpha, bend_trial, r_trial);
            double n_trial = vec_norm(r_trial, M);
            int dent_v = 0;
            int dent_bad = has_dent_with_first(bend_trial, &dent_v);

            if (n_trial < norm && !dent_bad) {
                /* tiny-step stall detection (matches python solver_lm) */
                double dnorm = 0; double xnorm = 0;
                for (int i = 0; i < N; i++) dnorm += delta[i]*delta[i];
                dnorm = sqrt(dnorm);
                /* x is implicit in bend_full; use ‖bend_var‖ */
                for (int e = 0; e < NE; e++) {
                    int j = VAR_OF_E[e];
                    if (j >= 0) xnorm += bend_full[e]*bend_full[e];
                }
                xnorm = sqrt(xnorm);
                if (dnorm < 1e-15 * (1.0 + xnorm) && n_trial >= norm) {
                    for (int e = 0; e < NE; e++) bend_full[e] = bend_trial[e];
                    memcpy(r, r_trial, M * sizeof(double));
                    norm = n_trial;
                    out.iters = it; out.final_resid = norm; out.final_lambda = lam;
                    snprintf(out.msg, sizeof(out.msg), "stalled");
                    goto done;
                }
                for (int e = 0; e < NE; e++) bend_full[e] = bend_trial[e];
                memcpy(r, r_trial, M * sizeof(double));
                norm = n_trial;
                lam = lam * lambda_down;
                if (lam < 1e-30) lam = 1e-30;
                accepted = true;
                break;
            }
            if (n_trial >= norm) {
                out.total_retries_residual++;
            } else {
                out.total_retries_dent++;
                if (CFG_DENT_TRACE) {
                    fprintf(stderr,
                        "dent_reject vertex=%d alpha_deg=%.6f flower_turn=%.6e iter=%d\n",
                        dent_v, alpha * 180.0 / PI,
                        vertex_turn(dent_v, bend_trial), it);
                }
            }
            lam = lam * lambda_up;
            if (lam > lambda_max) {
                out.iters = it; out.final_resid = norm; out.final_lambda = lam;
                snprintf(out.msg, sizeof(out.msg), "lambda_saturated");
                goto done;
            }
        }
        if (!accepted) {
            out.iters = it; out.final_resid = norm; out.final_lambda = lam;
            snprintf(out.msg, sizeof(out.msg), "lm_retries_exhausted");
            goto done;
        }
    }
    out.iters = max_iter; out.final_resid = norm; out.final_lambda = lam;
    out.success = (norm <= tol);
    snprintf(out.msg, sizeof(out.msg), out.success ? "max_iter_tol" : "max_iter");

done:
    if (out_iters_used) *out_iters_used = out.iters;
    free(r); free(r_trial); free(bend_trial);
    free(g); free(D); free(delta); free(bwork);
    if (J) free(J);
    if (A) free(A);
    if (Awork) free(Awork);
    return out;
}

/* ---------- Full-jump α-march α-carry controller ------------------------- */
typedef struct {
    /* status / failure_reason: short strings for TSV output */
    char   status[40];
    char   failure_reason[24];
    /* per-case stats */
    double final_alpha_deg;
    double final_resid;
    int    attempts;
    int    retreats;
    int    lm_iters_total;
    double final_lambda;
    double wall_secs;
    /* sequences (truncated at SEQ_BUF; '-' if empty) */
    char   accepted_fracs[8192];
    char   rejected_fracs[8192];
    char   alpha_seq_deg[8192];
} MarchResult;

static void seq_append(char *buf, size_t cap, double val) {
    /* Append "v" or ",v" depending on whether buf is "-" or not. */
    char tmp[64];
    snprintf(tmp, sizeof(tmp), "%.6g", val);
    size_t L = strlen(buf);
    if (L == 1 && buf[0] == '-') {
        if (cap > 1) {
            strncpy(buf, tmp, cap - 1);
            buf[cap - 1] = '\0';
        }
        return;
    }
    size_t need = strlen(tmp) + 1; /* leading comma */
    if (L + need + 1 > cap) return;  /* silently drop on overflow */
    buf[L] = ',';
    strncpy(buf + L + 1, tmp, cap - (L + 1) - 1);
    buf[cap - 1] = '\0';
}

/* alpha_jump_policy: 0 = full, 1 = half. (Default: full per spec.) */
typedef struct {
    double target_deg;
    double init_step_deg;
    double min_step_deg;
    int    max_quick_iter;
    double quick_tol;
    double final_tol;
    double lambda_init;
    int    max_attempts;
    double quick_min_drop;
    int    alpha_jump_policy;  /* 0 full, 1 half */
} MarchOpts;

static MarchResult march_run(double bend[], const MarchOpts *opts) {
    MarchResult res;
    memset(&res, 0, sizeof(res));
    snprintf(res.accepted_fracs, sizeof(res.accepted_fracs), "-");
    snprintf(res.rejected_fracs, sizeof(res.rejected_fracs), "-");
    snprintf(res.alpha_seq_deg,  sizeof(res.alpha_seq_deg),  "-");
    snprintf(res.failure_reason, sizeof(res.failure_reason), "-");
    snprintf(res.status, sizeof(res.status), "STALLED_POSITIVE_RESIDUAL");

    double target = opts->target_deg * PI / 180.0;
    double min_step = opts->min_step_deg * PI / 180.0;
    double step = opts->init_step_deg * PI / 180.0;
    double alpha_final = target * (1.0 - 1e-12);
    double alpha_curr = 0.0;
    double final_resid = INFINITY;
    int target_reached = 0;
    int n_step = 0, n_retreat = 0, n_attempts = 0;
    int lm_iters_total = 0;
    double final_lambda = opts->lambda_init;
    double cap_fraction = (opts->alpha_jump_policy == 0) ? 1.0 : 0.5;
    double try_fraction = cap_fraction;

    double t_start = now_secs();

    int *underflow = NULL; (void)underflow;
    int hit_underflow = 0;
    int hit_max_attempts = 0;

    for (int attempt = 0; attempt < opts->max_attempts; attempt++) {
        if (alpha_curr >= alpha_final) {
            target_reached = 1;
            break;
        }
        double gap = target - alpha_curr;
        double used;
        if (opts->alpha_jump_policy == 0) {
            used = gap * try_fraction;
        } else {
            double cap = gap * 0.5;
            used = (step < cap) ? step : cap;
            try_fraction = used / gap;
        }
        double alpha_try = alpha_curr + used;

        /* Stale-bend residual at α_try (input to the LM). */
        int M = 3 * N_RES;
        double *r_init = (double *)xcalloc(M, sizeof(double));
        holonomy_residual_quat(alpha_try, bend, r_init);
        double initial_resid = vec_norm(r_init, M);
        free(r_init);

        int iters_used = 0;
        LMOut out = solve_lm_inner(alpha_try, bend,
                                     opts->quick_tol, opts->max_quick_iter,
                                     20, /* max_lm_retries */
                                     opts->lambda_init,
                                     0.3, 10.0, 1e12, 1e-12,
                                     &iters_used);
        lm_iters_total += iters_used;
        n_attempts++;
        n_step++;
        final_lambda = out.final_lambda;

        double resid = out.final_resid;
        int accepted;
        if (opts->alpha_jump_policy == 0) {
            double denom = (initial_resid > 1e-30) ? initial_resid : 1e-30;
            int accepted_tol = (resid <= opts->quick_tol);
            int accepted_drop = (resid <= opts->quick_min_drop * denom);
            accepted = accepted_tol || accepted_drop;
        } else {
            accepted = (out.success && resid < opts->quick_tol);
        }

        if (accepted) {
            alpha_curr = alpha_try;
            seq_append(res.accepted_fracs, sizeof(res.accepted_fracs), try_fraction);
            seq_append(res.alpha_seq_deg, sizeof(res.alpha_seq_deg),
                       alpha_curr * 180.0 / PI);
            if (opts->alpha_jump_policy == 0) {
                try_fraction = cap_fraction;
            }
        } else {
            seq_append(res.rejected_fracs, sizeof(res.rejected_fracs), try_fraction);
            n_retreat++;
            if (opts->alpha_jump_policy == 0) {
                try_fraction *= 0.5;
                if (try_fraction * gap < min_step) {
                    hit_underflow = 1;
                    break;
                }
            } else {
                step *= 0.5;
                if (step < min_step) {
                    hit_underflow = 1;
                    break;
                }
            }
        }
    }

    /* Post-loop catch (lm_alpha_march_carry.py:421-422). */
    if (alpha_curr >= alpha_final) target_reached = 1;
    if (!target_reached && !hit_underflow) hit_max_attempts = 1;

    if (target_reached) {
        int iters_used = 0;
        LMOut out_final = solve_lm_inner(target, bend,
                                          opts->final_tol, 200,
                                          20, opts->lambda_init,
                                          0.3, 10.0, 1e12, 1e-12,
                                          &iters_used);
        lm_iters_total += iters_used;
        alpha_curr = target;
        final_resid = out_final.final_resid;
        final_lambda = out_final.final_lambda;
        if (final_resid <= opts->final_tol) {
            snprintf(res.status, sizeof(res.status), "SOLVER_TOL");
            snprintf(res.failure_reason, sizeof(res.failure_reason), "-");
        } else if (strcmp(out_final.msg, "lambda_saturated") == 0 ||
                   strcmp(out_final.msg, "lambda_saturated_linsolve") == 0) {
            snprintf(res.status, sizeof(res.status), "LAMBDA_SATURATED_POSITIVE_RESIDUAL");
            snprintf(res.failure_reason, sizeof(res.failure_reason), "LM_LAMBDA");
        } else if (strcmp(out_final.msg, "stalled") == 0) {
            snprintf(res.status, sizeof(res.status), "STALLED_POSITIVE_RESIDUAL");
            snprintf(res.failure_reason, sizeof(res.failure_reason), "LM_STALLED");
        } else if (strcmp(out_final.msg, "lm_retries_exhausted") == 0) {
            snprintf(res.status, sizeof(res.status), "STALLED_POSITIVE_RESIDUAL");
            snprintf(res.failure_reason, sizeof(res.failure_reason), "LM_RETRIES");
        } else if (strcmp(out_final.msg, "max_iter") == 0) {
            snprintf(res.status, sizeof(res.status), "STALLED_POSITIVE_RESIDUAL");
            snprintf(res.failure_reason, sizeof(res.failure_reason), "LM_MAX_ITER");
        } else {
            /* "tol"/"max_iter_tol" with resid > final_tol shouldn't happen,
             * but classify defensively. */
            snprintf(res.status, sizeof(res.status), "STALLED_POSITIVE_RESIDUAL");
            snprintf(res.failure_reason, sizeof(res.failure_reason), "LM_OTHER");
        }
    } else {
        snprintf(res.status, sizeof(res.status), "STALLED_POSITIVE_RESIDUAL");
        snprintf(res.failure_reason, sizeof(res.failure_reason),
                  hit_underflow ? "STEP_UNDERFLOW" : "MAX_ATTEMPTS");
    }

    res.final_alpha_deg = alpha_curr * 180.0 / PI;
    res.final_resid = final_resid;
    res.attempts = n_attempts;
    res.retreats = n_retreat;
    res.lm_iters_total = lm_iters_total;
    res.final_lambda = final_lambda;
    res.wall_secs = now_secs() - t_start;
    return res;
}

/* ---------- topology reset between cases (chunk mode) -------------------- */
static void clear_topology(void) {
    NV = NE = NF = 0;
    /* Lookups will be cleared at the top of build_topology(). */
    for (int v = 0; v <= MAXV; v++) {
        FLOWER_LEN[v] = 0;
        VERT_DEG[v] = 0;
    }
    NVAR = 0; N_INT = 0; N_RES = 0; N_GATE = 0;
#ifdef HAVE_SUPERLU
    sparse_quat_free();
#endif
}

/* ---------- end-to-end per case ------------------------------------------ */
/* Returns 0 on solver completion (regardless of status), -1 on early-reject
 * (HOROU_FAIL / TOPOLOGY_FAIL). On early reject, fills res.status. */
static int run_one_case(const char *netcode_str, const MarchOpts *opts,
                          MarchResult *res, double *bend_out, int *out_NV) {
    memset(res, 0, sizeof(*res));
    snprintf(res->accepted_fracs, sizeof(res->accepted_fracs), "-");
    snprintf(res->rejected_fracs, sizeof(res->rejected_fracs), "-");
    snprintf(res->alpha_seq_deg,  sizeof(res->alpha_seq_deg),  "-");
    snprintf(res->failure_reason, sizeof(res->failure_reason), "-");

    clear_topology();
    if (parse_netcode_str(netcode_str) < 0) {
        snprintf(res->status, sizeof(res->status), "TOPOLOGY_FAIL");
        snprintf(res->failure_reason, sizeof(res->failure_reason), "TOPOLOGY");
        return -1;
    }
    if (build_topology() < 0) {
        snprintf(res->status, sizeof(res->status), "TOPOLOGY_FAIL");
        snprintf(res->failure_reason, sizeof(res->failure_reason), "TOPOLOGY");
        return -1;
    }
    if (out_NV) *out_NV = NV;
    choose_base_face();
    if (CFG_SYSTEM == SYSTEM_ALLBENDS) setup_all_bends();
    else                                 setup_square();

    double *u = (double *)xcalloc(NV + 1, sizeof(double));
    if (horou_solve(u) < 0) {
        free(u);
        snprintf(res->status, sizeof(res->status), "HOROU_FAIL");
        snprintf(res->failure_reason, sizeof(res->failure_reason), "HOROU");
        return -1;
    }
    if (compute_bends_at_zero(u, bend_out) < 0) {
        free(u);
        snprintf(res->status, sizeof(res->status), "HOROU_FAIL");
        snprintf(res->failure_reason, sizeof(res->failure_reason), "HOROU");
        return -1;
    }
    free(u);

#ifdef HAVE_SUPERLU
    if (CFG_SOLVER == SOLVER_SPARSE) {
        if (sparse_quat_setup() != 0) {
            snprintf(res->status, sizeof(res->status), "EXCEPTION");
            snprintf(res->failure_reason, sizeof(res->failure_reason), "EXCEPTION");
            return -1;
        }
    }
#endif

    *res = march_run(bend_out, opts);
    return 0;
}

/* ---------- TSV writer --------------------------------------------------- */
static void write_tsv_header(FILE *fh) {
    fprintf(fh, "v\tclers\tnetcode_md5\tstatus\tfinal_alpha_deg\tfinal_resid\t"
                 "attempts\tretreats\taccepted_fractions\trejected_fractions\t"
                 "alpha_seq_deg\tlm_iters_total\tlambda_final\twall_secs\t"
                 "failure_reason\n");
}

/* Tiny MD5 — implemented inline so we don't add an external dep. RFC 1321. */
static void md5_to_hex(const unsigned char *in, size_t len, char hex[33]) {
    /* Public-domain compact MD5; sufficient as a stable column tag. */
    unsigned int K[64] = {
        0xd76aa478,0xe8c7b756,0x242070db,0xc1bdceee,0xf57c0faf,0x4787c62a,
        0xa8304613,0xfd469501,0x698098d8,0x8b44f7af,0xffff5bb1,0x895cd7be,
        0x6b901122,0xfd987193,0xa679438e,0x49b40821,0xf61e2562,0xc040b340,
        0x265e5a51,0xe9b6c7aa,0xd62f105d,0x02441453,0xd8a1e681,0xe7d3fbc8,
        0x21e1cde6,0xc33707d6,0xf4d50d87,0x455a14ed,0xa9e3e905,0xfcefa3f8,
        0x676f02d9,0x8d2a4c8a,0xfffa3942,0x8771f681,0x6d9d6122,0xfde5380c,
        0xa4beea44,0x4bdecfa9,0xf6bb4b60,0xbebfbc70,0x289b7ec6,0xeaa127fa,
        0xd4ef3085,0x04881d05,0xd9d4d039,0xe6db99e5,0x1fa27cf8,0xc4ac5665,
        0xf4292244,0x432aff97,0xab9423a7,0xfc93a039,0x655b59c3,0x8f0ccc92,
        0xffeff47d,0x85845dd1,0x6fa87e4f,0xfe2ce6e0,0xa3014314,0x4e0811a1,
        0xf7537e82,0xbd3af235,0x2ad7d2bb,0xeb86d391
    };
    unsigned int s[64] = {
        7,12,17,22, 7,12,17,22, 7,12,17,22, 7,12,17,22,
        5, 9,14,20, 5, 9,14,20, 5, 9,14,20, 5, 9,14,20,
        4,11,16,23, 4,11,16,23, 4,11,16,23, 4,11,16,23,
        6,10,15,21, 6,10,15,21, 6,10,15,21, 6,10,15,21
    };
    unsigned int a0=0x67452301, b0=0xefcdab89, c0=0x98badcfe, d0=0x10325476;

    /* Pre-processing: pad. */
    size_t orig = len;
    size_t newlen = ((len + 8) / 64 + 1) * 64;
    unsigned char *msg = (unsigned char *)calloc(newlen, 1);
    if (!msg) { hex[0] = '\0'; return; }
    memcpy(msg, in, len);
    msg[len] = 0x80;
    unsigned long long bits = (unsigned long long)orig * 8ULL;
    for (int i = 0; i < 8; i++) msg[newlen - 8 + i] = (unsigned char)(bits >> (8*i));

    for (size_t off = 0; off < newlen; off += 64) {
        unsigned int M[16];
        for (int i = 0; i < 16; i++) {
            M[i] = (unsigned int)msg[off + i*4]
                 | ((unsigned int)msg[off + i*4 + 1] << 8)
                 | ((unsigned int)msg[off + i*4 + 2] << 16)
                 | ((unsigned int)msg[off + i*4 + 3] << 24);
        }
        unsigned int A=a0, B=b0, C=c0, D=d0;
        for (int i = 0; i < 64; i++) {
            unsigned int F, g;
            if (i < 16) { F = (B & C) | (~B & D); g = i; }
            else if (i < 32) { F = (D & B) | (~D & C); g = (5*i + 1) % 16; }
            else if (i < 48) { F = B ^ C ^ D; g = (3*i + 5) % 16; }
            else { F = C ^ (B | ~D); g = (7*i) % 16; }
            unsigned int temp = D;
            D = C; C = B;
            unsigned int sum = A + F + K[i] + M[g];
            B = B + ((sum << s[i]) | (sum >> (32 - s[i])));
            A = temp;
        }
        a0 += A; b0 += B; c0 += C; d0 += D;
    }
    free(msg);
    unsigned char dig[16];
    for (int i = 0; i < 4; i++) {
        dig[i]    = (unsigned char)(a0 >> (8*i));
        dig[i+4]  = (unsigned char)(b0 >> (8*i));
        dig[i+8]  = (unsigned char)(c0 >> (8*i));
        dig[i+12] = (unsigned char)(d0 >> (8*i));
    }
    static const char *hexd = "0123456789abcdef";
    for (int i = 0; i < 16; i++) {
        hex[2*i]   = hexd[dig[i] >> 4];
        hex[2*i+1] = hexd[dig[i] & 0xf];
    }
    hex[32] = '\0';
}

static void write_tsv_row(FILE *fh, int v, const char *clers,
                            const char *netcode_str, const MarchResult *res) {
    char md5[33];
    md5_to_hex((const unsigned char *)netcode_str, strlen(netcode_str), md5);
    fprintf(fh,
        "%d\t%s\t%s\t%s\t%.10f\t%.6e\t%d\t%d\t%s\t%s\t%s\t%d\t%.6e\t%.6f\t%s\n",
        v, (clers && *clers) ? clers : "-", md5,
        res->status, res->final_alpha_deg, res->final_resid,
        res->attempts, res->retreats,
        res->accepted_fracs, res->rejected_fracs, res->alpha_seq_deg,
        res->lm_iters_total, res->final_lambda, res->wall_secs,
        res->failure_reason);
    fflush(fh);
}

/* ---------- bends-out (puffup-bends 1) ----------------------------------- */
static void write_puffup_bends(FILE *fh, const char *netcode_str,
                                  const double bend[], double alpha_rad) {
    fprintf(fh, "puffup-bends 1\n");
    fprintf(fh, "NV %d NE %d alpha_deg %.15g\n", NV, NE, alpha_rad * 180.0 / PI);
    fprintf(fh, "faces %s\n", netcode_str);
    fprintf(fh, "bends\n");
    /* Canonical homotopy_stage edge order: for u=1..NV walk NBR[u] in
     * face-scan / nbr_add insertion order, emit (u,w) iff u<w. Rebuilt
     * locally to match euclid_realize's parser. */
    static int canon_nbr[MAXV+1][MAXFLOWER];
    static int canon_nnbr[MAXV+1];
    static char canon_seen[MAXV+1][MAXV+1];
    memset(canon_nnbr, 0, sizeof(canon_nnbr));
    memset(canon_seen, 0, sizeof(canon_seen));
    #define CANON_ADD(u,w) do { \
        if (!canon_seen[(u)][(w)]) { \
            canon_nbr[(u)][canon_nnbr[(u)]++] = (w); \
            canon_seen[(u)][(w)] = 1; \
        } \
    } while (0)
    for (int i = 0; i < NF; i++) {
        int a = FACES[i][0], b = FACES[i][1], c = FACES[i][2];
        CANON_ADD(a,b); CANON_ADD(b,a);
        CANON_ADD(b,c); CANON_ADD(c,b);
        CANON_ADD(a,c); CANON_ADD(c,a);
    }
    for (int u = 1; u <= NV; u++) {
        for (int j = 0; j < canon_nnbr[u]; j++) {
            int w = canon_nbr[u][j];
            if (u < w) {
                int e = EDGE_IDX_LOOKUP[u][w];
                fprintf(fh, "%d %d %.15g\n", u, w, bend[e]);
            }
        }
    }
    #undef CANON_ADD
}

/* ---------- FD spot-check ----------------------------------------------- */
static int fd_spot_check(double alpha) {
    /* For each var-edge e, compare analytic Jacobian column to FD column.
     * Print max-relative-error per column and overall pass/fail at 1e-6. */
    double *u = (double *)xcalloc(NV + 1, sizeof(double));
    double *bend = (double *)xcalloc(NE, sizeof(double));
    if (horou_solve(u) < 0) { free(u); free(bend); return -1; }
    if (compute_bends_at_zero(u, bend) < 0) { free(u); free(bend); return -1; }
    free(u);
    int rows = 3 * N_RES;
    double *Ja = (double *)xcalloc((size_t)rows * (size_t)NVAR, sizeof(double));
    double *Jf = (double *)xcalloc((size_t)rows * (size_t)NVAR, sizeof(double));
    quat_analytic_jacobian_dense(alpha, bend, Ja);
    quat_fd_jacobian_dense     (alpha, bend, Jf);

    double max_abs_err = 0.0;
    double max_rel_err = 0.0;
    int worst_row = -1, worst_col = -1;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < NVAR; j++) {
            double a = Ja[i * NVAR + j];
            double f = Jf[i * NVAR + j];
            double err = fabs(a - f);
            double scale = fabs(a) + fabs(f);
            double rel = (scale > 1e-12) ? err / scale : err;
            if (err > max_abs_err) max_abs_err = err;
            if (rel > max_rel_err) { max_rel_err = rel; worst_row = i; worst_col = j; }
        }
    }
    fprintf(stdout,
        "fd_spot_check: NV=%d NE=%d N_RES=%d NVAR=%d alpha_deg=%.6f "
        "max_abs_err=%.3e max_rel_err=%.3e worst=(row=%d,col=%d) verdict=%s\n",
        NV, NE, N_RES, NVAR, alpha * 180.0 / PI,
        max_abs_err, max_rel_err, worst_row, worst_col,
        max_rel_err <= 1e-6 ? "PASS" : "FAIL");
    free(Ja); free(Jf); free(bend);
    return (max_rel_err <= 1e-6) ? 0 : 1;
}

/* ---------- main --------------------------------------------------------- */
static void usage(FILE *fh) {
    fprintf(fh,
        "usage: lm_march_c [flags] < netcode_facelist\n"
        "       lm_march_c --input-list FILE --output-tsv OUT [flags]\n"
        "flags:\n"
        "  --system all-bends|square        (default: all-bends)\n"
        "  --residual quat|matrix           (default: quat; matrix not implemented)\n"
        "  --alpha-jump-policy full|half    (default: full)\n"
        "  --target-deg N                   (default: 60)\n"
        "  --init-step-deg N                (default: 1.0)\n"
        "  --min-step-deg N                 (default: 1e-4)\n"
        "  --max-quick-iter N               (default: 10)\n"
        "  --quick-tol R                    (default: 1e-3)\n"
        "  --final-tol R                    (default: 1e-12)\n"
        "  --quick-min-drop R               (default: 1e-2)\n"
        "  --lambda-init R                  (default: 1.0)\n"
        "  --max-attempts N                 (default: 2000)\n"
        "  --solver dense|sparse            (default: sparse if HAVE_SUPERLU)\n"
        "  --bends-out FILE                 (single-case mode only)\n"
        "  --bends-out-dir DIR              (chunk mode: dump <md5>.bends per SOLVER_TOL case)\n"
        "  --dent-gate-trace                (print dent rejections to stderr)\n"
        "  --fd-spot-check --alpha-deg N    (analytic vs FD on horou-ideal bends)\n"
        "  --input-list FILE --output-tsv F (chunk mode; CLERS lines via --clers-list optional)\n"
        "  --clers-list FILE                (parallel to --input-list; one CLERS string per line)\n"
        "  --case-tag NAME                  (single-case mode TSV row label)\n"
        "  --output-tsv FILE                (single-case mode: write one TSV row to FILE)\n"
    );
}

int main(int argc, char **argv) {
    tau_const = 2.0 * PI;
    MarchOpts opts;
    opts.target_deg = 60.0;
    opts.init_step_deg = 1.0;
    opts.min_step_deg = 1e-4;
    opts.max_quick_iter = 10;
    opts.quick_tol = 1e-3;
    opts.final_tol = 1e-12;
    opts.lambda_init = 1.0;
    opts.max_attempts = 2000;
    opts.quick_min_drop = 1e-2;
    opts.alpha_jump_policy = 0;  /* full */
    const char *bends_out = NULL;
    const char *bends_out_dir = NULL;
    const char *input_list = NULL;
    const char *output_tsv = NULL;
    const char *clers_list = NULL;
    const char *case_tag = NULL;
    int do_fd_spot = 0;
    double fd_alpha_deg = 30.0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--system") == 0 && i + 1 < argc) {
            const char *s = argv[++i];
            if      (strcmp(s, "all-bends") == 0) CFG_SYSTEM = SYSTEM_ALLBENDS;
            else if (strcmp(s, "square") == 0)    CFG_SYSTEM = SYSTEM_SQUARE;
            else die("--system must be 'all-bends' or 'square'");
        } else if (strcmp(argv[i], "--residual") == 0 && i + 1 < argc) {
            const char *s = argv[++i];
            if      (strcmp(s, "quat") == 0)   CFG_RESIDUAL = RESIDUAL_QUAT;
            else if (strcmp(s, "matrix") == 0) die("--residual matrix not implemented in lm_march_c (use puffup_c_lm)");
            else die("--residual must be 'quat' or 'matrix'");
        } else if (strcmp(argv[i], "--alpha-jump-policy") == 0 && i + 1 < argc) {
            const char *s = argv[++i];
            if      (strcmp(s, "full") == 0) opts.alpha_jump_policy = 0;
            else if (strcmp(s, "half") == 0) opts.alpha_jump_policy = 1;
            else die("--alpha-jump-policy must be 'full' or 'half'");
        } else if (strcmp(argv[i], "--target-deg") == 0 && i + 1 < argc) {
            opts.target_deg = atof(argv[++i]);
        } else if (strcmp(argv[i], "--init-step-deg") == 0 && i + 1 < argc) {
            opts.init_step_deg = atof(argv[++i]);
        } else if (strcmp(argv[i], "--min-step-deg") == 0 && i + 1 < argc) {
            opts.min_step_deg = atof(argv[++i]);
        } else if (strcmp(argv[i], "--max-quick-iter") == 0 && i + 1 < argc) {
            opts.max_quick_iter = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--quick-tol") == 0 && i + 1 < argc) {
            opts.quick_tol = atof(argv[++i]);
        } else if (strcmp(argv[i], "--final-tol") == 0 && i + 1 < argc) {
            opts.final_tol = atof(argv[++i]);
        } else if (strcmp(argv[i], "--quick-min-drop") == 0 && i + 1 < argc) {
            opts.quick_min_drop = atof(argv[++i]);
        } else if (strcmp(argv[i], "--lambda-init") == 0 && i + 1 < argc) {
            opts.lambda_init = atof(argv[++i]);
        } else if (strcmp(argv[i], "--max-attempts") == 0 && i + 1 < argc) {
            opts.max_attempts = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--solver") == 0 && i + 1 < argc) {
            const char *s = argv[++i];
            if      (strcmp(s, "dense")  == 0) CFG_SOLVER = SOLVER_DENSE;
            else if (strcmp(s, "sparse") == 0) {
#ifdef HAVE_SUPERLU
                CFG_SOLVER = SOLVER_SPARSE;
#else
                die("--solver sparse requires building with -DHAVE_SUPERLU");
#endif
            } else die("--solver must be 'dense' or 'sparse'");
        } else if (strcmp(argv[i], "--bends-out") == 0 && i + 1 < argc) {
            bends_out = argv[++i];
        } else if (strcmp(argv[i], "--bends-out-dir") == 0 && i + 1 < argc) {
            bends_out_dir = argv[++i];
        } else if (strcmp(argv[i], "--dent-gate-trace") == 0) {
            CFG_DENT_TRACE = 1;
        } else if (strcmp(argv[i], "--fd-spot-check") == 0) {
            do_fd_spot = 1;
        } else if (strcmp(argv[i], "--alpha-deg") == 0 && i + 1 < argc) {
            fd_alpha_deg = atof(argv[++i]);
        } else if (strcmp(argv[i], "--input-list") == 0 && i + 1 < argc) {
            input_list = argv[++i];
        } else if (strcmp(argv[i], "--output-tsv") == 0 && i + 1 < argc) {
            output_tsv = argv[++i];
        } else if (strcmp(argv[i], "--clers-list") == 0 && i + 1 < argc) {
            clers_list = argv[++i];
        } else if (strcmp(argv[i], "--case-tag") == 0 && i + 1 < argc) {
            case_tag = argv[++i];
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            usage(stdout); return 0;
        } else {
            usage(stderr);
            die("unknown command-line option");
        }
    }

    /* FD spot-check mode: read facelist on stdin, build, run check. */
    if (do_fd_spot) {
        /* slurp stdin */
        size_t cap = 4096, len = 0;
        char *buf = (char *)xcalloc(cap, 1);
        int ch;
        while ((ch = fgetc(stdin)) != EOF) {
            if (len + 1 >= cap) { cap *= 2; buf = (char *)realloc(buf, cap); if (!buf) die("realloc"); }
            buf[len++] = (char)ch;
        }
        buf[len] = '\0';
        clear_topology();
        if (parse_netcode_str(buf) < 0) die("fd_spot_check: bad netcode");
        if (build_topology() < 0) die("fd_spot_check: bad topology");
        choose_base_face();
        if (CFG_SYSTEM == SYSTEM_ALLBENDS) setup_all_bends();
        else                                 setup_square();
        int rc = fd_spot_check(fd_alpha_deg * PI / 180.0);
        free(buf);
        return rc;
    }

    /* Chunk mode */
    if (input_list) {
        if (!output_tsv) die("--input-list requires --output-tsv");
        FILE *fin = fopen(input_list, "r");
        if (!fin) die("cannot open --input-list");
        FILE *fclers = NULL;
        if (clers_list) {
            fclers = fopen(clers_list, "r");
            if (!fclers) die("cannot open --clers-list");
        }
        FILE *fout = fopen(output_tsv, "w");
        if (!fout) die("cannot open --output-tsv");
        write_tsv_header(fout);

        char line[8192];
        char clers_line[1024];
        long n_done = 0, n_pass = 0;
        while (fgets(line, sizeof(line), fin)) {
            /* strip newline */
            size_t L = strlen(line);
            while (L > 0 && (line[L-1] == '\n' || line[L-1] == '\r')) line[--L] = '\0';
            if (L == 0) continue;
            if (line[0] == '#') continue;

            const char *clers = "-";
            if (fclers) {
                if (!fgets(clers_line, sizeof(clers_line), fclers)) clers_line[0] = '\0';
                size_t M = strlen(clers_line);
                while (M > 0 && (clers_line[M-1] == '\n' || clers_line[M-1] == '\r'))
                    clers_line[--M] = '\0';
                clers = clers_line;
            }

            MarchResult res;
            double *bend = (double *)xcalloc(MAXE, sizeof(double));
            int nv_out = 0;
            int rc = run_one_case(line, &opts, &res, bend, &nv_out);
            if (bends_out_dir && rc == 0 && strcmp(res.status, "SOLVER_TOL") == 0) {
                char md5[33];
                md5_to_hex((const unsigned char *)line, strlen(line), md5);
                char path[2048];
                snprintf(path, sizeof(path), "%s/%s.bends", bends_out_dir, md5);
                FILE *bf = fopen(path, "w");
                if (!bf) {
                    fprintf(stderr, "lm_march_c: fopen %s: %s\n", path, strerror(errno));
                } else {
                    write_puffup_bends(bf, line, bend, opts.target_deg * PI / 180.0);
                    if (fclose(bf) != 0) {
                        fprintf(stderr, "lm_march_c: fclose %s: %s\n", path, strerror(errno));
                    }
                }
            }
            free(bend);
            write_tsv_row(fout, nv_out, clers, line, &res);
            n_done++;
            if (strcmp(res.status, "SOLVER_TOL") == 0) n_pass++;
        }
        fclose(fin);
        if (fclers) fclose(fclers);
        fclose(fout);
        fprintf(stderr, "lm_march_c: chunk done — %ld cases, %ld SOLVER_TOL\n",
                  n_done, n_pass);
        return 0;
    }

    /* Single-case mode */
    {
        size_t cap = 4096, len = 0;
        char *buf = (char *)xcalloc(cap, 1);
        int ch;
        while ((ch = fgetc(stdin)) != EOF) {
            if (len + 1 >= cap) { cap *= 2; buf = (char *)realloc(buf, cap); if (!buf) die("realloc"); }
            buf[len++] = (char)ch;
        }
        buf[len] = '\0';

        MarchResult res;
        double *bend = (double *)xcalloc(MAXE, sizeof(double));
        int nv_out = 0;
        int rc = run_one_case(buf, &opts, &res, bend, &nv_out);

        /* Stdout: stats: lines (mirrors python format) for diff-friendliness. */
        printf("clers=%s V=%d E=%d system=%s residual=%s\n",
                case_tag ? case_tag : "-", NV, NE,
                CFG_SYSTEM == SYSTEM_ALLBENDS ? "all-bends" : "square",
                CFG_RESIDUAL == RESIDUAL_QUAT ? "quat" : "matrix");
        printf("summary: status=%s final_alpha_deg=%.10f final_resid=%.6e "
               "attempts=%d retreats=%d wall=%.4fs\n",
                res.status, res.final_alpha_deg, res.final_resid,
                res.attempts, res.retreats, res.wall_secs);
        printf("stats: accepted_fraction_seq=%s\n", res.accepted_fracs);
        printf("stats: rejected_fraction_seq=%s\n", res.rejected_fracs);
        printf("stats: alpha_seq_deg=%s\n", res.alpha_seq_deg);
        printf("stats: failure_reason=%s\n", res.failure_reason);

        if (output_tsv) {
            FILE *fh = fopen(output_tsv, "w");
            if (!fh) die("cannot open --output-tsv");
            write_tsv_header(fh);
            write_tsv_row(fh, nv_out, case_tag ? case_tag : "-", buf, &res);
            fclose(fh);
        }

        if (bends_out && rc == 0 && strcmp(res.status, "SOLVER_TOL") == 0) {
            FILE *fh = fopen(bends_out, "w");
            if (!fh) die("cannot open --bends-out");
            write_puffup_bends(fh, buf, bend, opts.target_deg * PI / 180.0);
            fclose(fh);
            fprintf(stderr, "lm_march_c: wrote bends -> %s\n", bends_out);
        }

        free(buf); free(bend);
        return (rc == 0 && strcmp(res.status, "SOLVER_TOL") == 0) ? 0 : 1;
    }
}
