/*
  puffup_c_lm.c — experimental LM solver with vertexwish start.

  Mirrors python/puffup.py + python/experiment_lm/vertexwish.py + the
  solver_lm + dent gate from the experiment branch. NOT for production:
  parallel to src/puffup_c.c, intentionally standalone (no shared code).

  Build:
      cc -O3 -std=c11 -Wall -Wextra -pedantic -o src/puffup_c_lm src/puffup_c_lm.c -lm

  CLI:
      puffup_c_lm < netcode_file [--alpha RAD] [--lambda-init L] [--tol T] [--max-iter N]

  Reads one netcode line on stdin: "1,2,3;1,3,4;...". Computes vertexwish
  start, runs LM at the given α, prints status to stdout.

  Math conventions match python/puffup.py:
    movemat(α, β) = matz(α) · matx(-β)
    vertex residual = three off-diagonals of cumulative holomat at non-base verts
    dent gate = vertex_turn(v) ≥ 0 at non-base verts; vertex_turn = sum of bends
    LM step = (J^T J + λ · diag(D)) δ = -J^T r,
              D_i = max(diag(J^T J)_i, tiny_rel · max(diag), 1e-30)
    accept if strict residual decrease AND no dent introduced

  No homotopy. Single shot at α.
*/

#define _POSIX_C_SOURCE 200809L
#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ---------- size limits (raise if you need v > 100) ----------------------- */
#define MAXV       200
#define MAXE       (3 * MAXV)   /* 3V - 6 < 3V */
#define MAXF       (2 * MAXV)   /* 2V - 4 < 2V */
#define MAXFLOWER  20

/* ---------- topology globals ---------------------------------------------- */
static int NV, NE, NF;
static int FACES[MAXF][3];
static int EDGE_A[MAXE], EDGE_B[MAXE];
static int EDGE_IDX_LOOKUP[MAXV + 1][MAXV + 1];  /* edge_idx[u][v]; -1 if not an edge */
static int VERT_DEG[MAXV + 1];
static int FLOWER_LEN[MAXV + 1];
static int FLOWER_THIRD[MAXV + 1][MAXFLOWER];     /* third vertex of t-th flower face */
static int FLOWER_E[MAXV + 1][MAXFLOWER];         /* edge_idx for that step */
static int DIRECTED_FACE[MAXV + 1][MAXV + 1];     /* directed edge (a,b) → face_idx, or -1 */

static int BASE[3];
static bool IS_BASE_E[MAXE];
static int NVAR;
static int VAR_OF_E[MAXE];   /* non-base edges → var index; base edges → -1 */
static int N_INT;            /* non-base vertex count = NV - 3 */
static int INT_VS[MAXV];

static double tau_const;     /* 2 π */

static void die(const char *msg) { fprintf(stderr, "error: %s\n", msg); exit(2); }

static void *xcalloc(size_t n, size_t sz) {
    void *p = calloc(n, sz);
    if (!p) die("calloc failed");
    return p;
}

#ifndef PI
#define PI 3.14159265358979323846
#endif

/* ---------- netcode parser ------------------------------------------------- */
static int parse_netcode_stream(FILE *fh) {
    /* Read one face per token "a,b,c"; tokens separated by ';' or whitespace. */
    NF = 0;
    int a = 0, b = 0, c = 0, slot = 0, val = 0, have = 0;
    int ch;
    int max_v = 0;
    while ((ch = fgetc(fh)) != EOF) {
        if (isdigit(ch)) {
            val = val * 10 + (ch - '0');
            have = 1;
        } else if (ch == ',') {
            if (!have) die("netcode: empty field");
            if (slot == 0) a = val;
            else if (slot == 1) b = val;
            else die("netcode: too many commas in face");
            slot++; val = 0; have = 0;
        } else if (ch == ';' || ch == '\n' || ch == '\r' || ch == ' ' || ch == '\t') {
            if (have) {
                if (slot != 2) die("netcode: face with wrong number of fields");
                c = val;
                if (NF >= MAXF) die("netcode: too many faces");
                if (a < 1 || b < 1 || c < 1) die("netcode: vertex IDs must be >= 1");
                FACES[NF][0] = a; FACES[NF][1] = b; FACES[NF][2] = c;
                if (a > max_v) max_v = a;
                if (b > max_v) max_v = b;
                if (c > max_v) max_v = c;
                NF++;
            }
            slot = 0; val = 0; have = 0;
        } else {
            die("netcode: unexpected character");
        }
    }
    if (have) {
        if (slot != 2) die("netcode: trailing face with wrong field count");
        c = val;
        if (NF >= MAXF) die("netcode: too many faces");
        if (a < 1 || b < 1 || c < 1) die("netcode: vertex IDs must be >= 1");
        FACES[NF][0] = a; FACES[NF][1] = b; FACES[NF][2] = c;
        if (a > max_v) max_v = a;
        if (b > max_v) max_v = b;
        if (c > max_v) max_v = c;
        NF++;
    }
    if (NF == 0) die("netcode: empty");
    NV = max_v;
    if (NV > MAXV) die("netcode: NV exceeds MAXV (raise MAXV)");
    return 0;
}

/* ---------- topology build ------------------------------------------------- */
static int add_edge_or_lookup(int u, int v) {
    int a = u < v ? u : v;
    int b = u < v ? v : u;
    int idx = EDGE_IDX_LOOKUP[a][b];
    if (idx >= 0) return idx;
    if (NE >= MAXE) die("too many edges");
    EDGE_A[NE] = a; EDGE_B[NE] = b;
    EDGE_IDX_LOOKUP[a][b] = NE;
    return NE++;
}

static void build_topology(void) {
    /* edges, edge index lookup */
    NE = 0;
    for (int v = 0; v <= NV; v++)
        for (int u = 0; u <= NV; u++)
            EDGE_IDX_LOOKUP[v][u] = -1;
    for (int v = 0; v <= NV; v++)
        for (int u = 0; u <= NV; u++)
            DIRECTED_FACE[v][u] = -1;

    for (int fi = 0; fi < NF; fi++) {
        int a = FACES[fi][0], b = FACES[fi][1], c = FACES[fi][2];
        if (a == b || b == c || a == c) die("degenerate face");
        if (DIRECTED_FACE[a][b] >= 0) die("orientation inconsistent (directed edge appears twice)");
        if (DIRECTED_FACE[b][c] >= 0) die("orientation inconsistent");
        if (DIRECTED_FACE[c][a] >= 0) die("orientation inconsistent");
        DIRECTED_FACE[a][b] = fi;
        DIRECTED_FACE[b][c] = fi;
        DIRECTED_FACE[c][a] = fi;
        add_edge_or_lookup(a, b);
        add_edge_or_lookup(b, c);
        add_edge_or_lookup(c, a);
    }

    /* check every undirected edge has both directions covered */
    for (int e = 0; e < NE; e++) {
        int a = EDGE_A[e], b = EDGE_B[e];
        if (DIRECTED_FACE[a][b] < 0 || DIRECTED_FACE[b][a] < 0)
            die("not a closed surface (some directed edge unmatched)");
    }

    /* per-vertex degree (count of edges) */
    for (int v = 0; v <= NV; v++) VERT_DEG[v] = 0;
    for (int e = 0; e < NE; e++) {
        VERT_DEG[EDGE_A[e]]++;
        VERT_DEG[EDGE_B[e]]++;
    }

    /* per-vertex flower: cyclic list of faces around v.
       Convention from puffup.py: from face (v, b, c) step to face containing (v, c). */
    for (int v = 1; v <= NV; v++) {
        /* find any face containing v */
        int start = -1;
        for (int fi = 0; fi < NF; fi++) {
            if (FACES[fi][0] == v || FACES[fi][1] == v || FACES[fi][2] == v) {
                start = fi; break;
            }
        }
        if (start < 0) die("vertex with no incident face");
        /* walk */
        FLOWER_LEN[v] = 0;
        int cur = start;
        int seen_count = 0;
        while (1) {
            if (seen_count >= MAXFLOWER) die("flower exceeds MAXFLOWER (raise it)");
            FLOWER_LEN[v]++;
            int third;
            if (FACES[cur][0] == v) third = FACES[cur][2];
            else if (FACES[cur][1] == v) third = FACES[cur][0];
            else third = FACES[cur][1];
            FLOWER_THIRD[v][seen_count] = third;
            FLOWER_E[v][seen_count] = EDGE_IDX_LOOKUP[v < third ? v : third][v < third ? third : v];
            seen_count++;
            int nxt = DIRECTED_FACE[v][third];
            if (nxt < 0) die("flower walk failed: missing directed edge");
            if (nxt == start) break;
            cur = nxt;
        }
        if (FLOWER_LEN[v] != VERT_DEG[v]) {
            fprintf(stderr, "vertex %d: flower length %d != degree %d\n",
                    v, FLOWER_LEN[v], VERT_DEG[v]);
            die("non-manifold or bad-link vertex");
        }
    }
}

/* ---------- base face / vars ---------------------------------------------- */
static void choose_base_face(void) {
    BASE[0] = FACES[0][0];
    BASE[1] = FACES[0][1];
    BASE[2] = FACES[0][2];

    for (int e = 0; e < NE; e++) IS_BASE_E[e] = false;
    int e01 = EDGE_IDX_LOOKUP[BASE[0] < BASE[1] ? BASE[0] : BASE[1]][BASE[0] < BASE[1] ? BASE[1] : BASE[0]];
    int e12 = EDGE_IDX_LOOKUP[BASE[1] < BASE[2] ? BASE[1] : BASE[2]][BASE[1] < BASE[2] ? BASE[2] : BASE[1]];
    int e20 = EDGE_IDX_LOOKUP[BASE[2] < BASE[0] ? BASE[2] : BASE[0]][BASE[2] < BASE[0] ? BASE[0] : BASE[2]];
    IS_BASE_E[e01] = IS_BASE_E[e12] = IS_BASE_E[e20] = true;

    NVAR = 0;
    for (int e = 0; e < NE; e++) {
        VAR_OF_E[e] = IS_BASE_E[e] ? -1 : NVAR++;
    }

    N_INT = 0;
    for (int v = 1; v <= NV; v++) {
        if (v != BASE[0] && v != BASE[1] && v != BASE[2]) {
            INT_VS[N_INT++] = v;
        }
    }
}

/* ---------- horou (port of src/horou_c.c, adapted to local topology) ------ */
/* Reference: src/horou_c.c and the horou block in src/puffup_c.c (lines
   242-417). All algebra preserved; only the topology accesses are swapped
   to use FLOWER_THIRD / DIRECTED_FACE / VERT_DEG. */
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

/* Generic in-place LU + back-sub. row-major n×n in J, RHS in b. */
static int horou_lu_solve(double *J, double *b, int n) {
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

/* horou solver: returns u_out[v-1] for v=2..NV. u_out[0] = NaN
   (vertex 1 is the point at infinity). Boundary vertices (neighbors of 1)
   are pinned to u=1; interior vertices solved by Newton with petal-sum
   = 2π constraint per interior vertex.
   Returns 0 on success, -1 on failure. */
static int horou_solve(double *u_out) {
    const double TAU = 2.0 * PI;
    int *bndry = (int *)xcalloc(NV + 2, sizeof(int));
    int *int_idx = (int *)xcalloc(NV + 2, sizeof(int));
    for (int v = 0; v <= NV + 1; v++) int_idx[v] = -1;

    /* boundary = neighbors of vertex 1 (cyclic_nbrs of 1) */
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

    /* per-interior-vertex cyclic neighbor ring; reuse FLOWER_THIRD for that. */
    int *ringlen = (int *)xcalloc(NV + 2, sizeof(int));
    for (int i = 0; i < n_int; i++) {
        int v = interior[i];
        ringlen[v] = FLOWER_LEN[v];
    }

    double *xvec = (double *)xcalloc(NV + 2, sizeof(double));
    for (int i = 0; i < n_int; i++) xvec[i] = 1.0;

    /* face list excluding any face containing vertex 1 (used to enforce
       triangle inequality during line search). */
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
        if (horou_lu_solve(HJ, Hdx, n_int) < 0) {
            free(HJ); free(HF); free(Hdx);
            free(ff_a); free(ff_b); free(ff_c);
            free(xvec); free(ringlen); free(interior); free(int_idx); free(bndry);
            return -1;
        }

        double step = 1.0;
        int found = 0;
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

/* α=0 dihedral per edge from u: bend = π − interior dihedral.
   Adapted from src/puffup_c.c's compute_bends_at_zero. */
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
        /* face containing v in flower at t. We need (a,b,c) with a==v
           (FLOWER_THIRD gives c, but we also need b). Find the face
           explicitly via DIRECTED_FACE walk: c = FLOWER_THIRD[v][t]. The
           corresponding face contains directed edge (v, c). */
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

/* ---------- math: matz(α) · matx(-β) -------------------------------------- */
typedef double M3[3][3];

static inline void mat_eye(M3 A) {
    for (int i=0;i<3;i++) for (int j=0;j<3;j++) A[i][j] = (i==j) ? 1.0 : 0.0;
}
static inline void matmul(const M3 A, const M3 B, M3 C) {
    M3 T;
    for (int i=0;i<3;i++)
        for (int j=0;j<3;j++) {
            double s = 0;
            for (int k=0;k<3;k++) s += A[i][k] * B[k][j];
            T[i][j] = s;
        }
    memcpy(C, T, sizeof(M3));
}
static inline void movemat(double alpha, double beta, M3 M) {
    /* matz(α) · matx(-β) */
    double ca=cos(alpha), sa=sin(alpha);
    double cb=cos(beta),  sb=sin(beta);  /* matx(-β) uses cos β, sin β with sign flips */
    M[0][0] =  ca;  M[0][1] = -sa*cb;  M[0][2] = -sa*sb;
    M[1][0] =  sa;  M[1][1] =  ca*cb;  M[1][2] =  ca*sb;
    M[2][0] = 0.0;  M[2][1] = -sb;     M[2][2] =  cb;
}
/* dM/dβ at given α,β. Lifted from canonical homotopy_stage puffup_c.c
 * dmovemat_dbeta. d/dβ matx(-β) = -Jx · matx(-β), so
 *   d/dβ [matz(α) · matx(-β)] = matz(α) · [[0,0,0],[0,-sb,cb],[0,-cb,-sb]]. */
static inline void dmovemat_dbeta(double alpha, double beta, M3 D) {
    double ca=cos(alpha), sa=sin(alpha);
    double cb=cos(beta),  sb=sin(beta);
    D[0][0] = 0.0;
    D[0][1] =  sa*sb;
    D[0][2] = -sa*cb;
    D[1][0] = 0.0;
    D[1][1] = -ca*sb;
    D[1][2] =  ca*cb;
    D[2][0] = 0.0;
    D[2][1] = -cb;
    D[2][2] = -sb;
}

/* ---------- residual / Jacobian ------------------------------------------- */
static void holonomy_residual(double alpha, const double bend[], double r[]) {
    /* r has length 3 * N_INT */
    int rk = 0;
    for (int t = 0; t < N_INT; t++) {
        int v = INT_VS[t];
        M3 M; mat_eye(M);
        int kk = FLOWER_LEN[v];
        for (int s = 0; s < kk; s++) {
            int e = FLOWER_E[v][s];
            M3 step; movemat(alpha, bend[e], step);
            matmul(M, step, M);
        }
        r[rk++] = M[0][1];
        r[rk++] = M[0][2];
        r[rk++] = M[1][2];
    }
}

/* finite-difference Jacobian. J is row-major (3*N_INT) × NVAR.
   bend_tmp, r_plus, r_minus are caller-supplied scratch buffers. None
   may alias caller's residual buffer. Kept for fallback/comparison;
   analytical_jacobian is the production path. */
static void fd_jacobian(double alpha, const double bend[], double *J,
                         double *bend_tmp,
                         double *r_plus, double *r_minus) {
    int M = 3 * N_INT;
    for (int e = 0; e < NE; e++) bend_tmp[e] = bend[e];
    for (int e = 0; e < NE; e++) {
        int j = VAR_OF_E[e];
        if (j < 0) continue;
        double original = bend_tmp[e];
        double h = 1e-7 * (fabs(original) > 1.0 ? fabs(original) : 1.0);
        bend_tmp[e] = original + h;
        holonomy_residual(alpha, bend_tmp, r_plus);
        bend_tmp[e] = original - h;
        holonomy_residual(alpha, bend_tmp, r_minus);
        bend_tmp[e] = original;
        for (int i = 0; i < M; i++) {
            J[i * NVAR + j] = (r_plus[i] - r_minus[i]) / (2.0 * h);
        }
    }
}

/* Analytical Jacobian — lifted from canonical homotopy_stage puffup_c.c
 * (functions vertex_flower_prefixes / jac_contrib_at / analytical_jacobian).
 * Produces the same (3*N_INT) × NVAR row-major dense matrix that
 * fd_jacobian builds, but in one pass with closed-form derivatives. */
static void vertex_flower_prefixes(int v, double alpha, const double bend[],
                                    int *k_out, M3 *Ms, M3 *P, M3 *S) {
    int k = FLOWER_LEN[v];
    *k_out = k;
    for (int t = 0; t < k; t++) movemat(alpha, bend[FLOWER_E[v][t]], Ms[t]);
    mat_eye(P[0]);
    for (int t = 0; t < k; t++) matmul(P[t], Ms[t], P[t+1]);
    mat_eye(S[k]);
    for (int t = k-1; t >= 0; t--) matmul(Ms[t], S[t+1], S[t]);
}
static void jac_contrib_at(double alpha, double beta_e,
                            const M3 P_t, const M3 S_tnext,
                            double *v01, double *v02, double *v12) {
    M3 dM, tmp, contrib;
    dmovemat_dbeta(alpha, beta_e, dM);
    matmul(P_t, dM, tmp);
    matmul(tmp, S_tnext, contrib);
    *v01 = contrib[0][1];
    *v02 = contrib[0][2];
    *v12 = contrib[1][2];
}
static void analytical_jacobian(double alpha, const double bend[], double *J) {
    int rows = 3 * N_INT;
    memset(J, 0, (size_t)rows * (size_t)NVAR * sizeof(double));
    M3 Ms[MAXFLOWER], P[MAXFLOWER+1], S[MAXFLOWER+1];
    for (int i = 0; i < N_INT; i++) {
        int v = INT_VS[i];
        int k;
        vertex_flower_prefixes(v, alpha, bend, &k, Ms, P, S);
        for (int t = 0; t < k; t++) {
            int e = FLOWER_E[v][t];
            int col = VAR_OF_E[e];
            if (col < 0) continue;
            double v01, v02, v12;
            jac_contrib_at(alpha, bend[e], P[t], S[t+1], &v01, &v02, &v12);
            int row = 3 * i;
            J[(row+0)*NVAR + col] += v01;
            J[(row+1)*NVAR + col] += v02;
            J[(row+2)*NVAR + col] += v12;
        }
    }
}

/* ---------- dent gate ----------------------------------------------------- */
static double vertex_turn(int v, const double bend[]) {
    double s = 0;
    int k = FLOWER_LEN[v];
    for (int t = 0; t < k; t++) s += bend[FLOWER_E[v][t]];
    return s;
}
static bool has_dent(const double bend[]) {
    /* skip base-face vertices to match python solve_homotopy's gate */
    for (int t = 0; t < N_INT; t++) {
        if (vertex_turn(INT_VS[t], bend) < 0.0) return true;
    }
    return false;
}

/* ---------- linear algebra: dense Gaussian elimination -------------------- */
static int dense_solve(double *A, double *b, int n) {
    /* Solve A x = b in place. Partial pivoting. Returns 0 on success, 1 on singular. */
    for (int k = 0; k < n; k++) {
        /* find pivot */
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
    /* back-substitute (b ← solution) */
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

/* ---------- vertexwish start ---------------------------------------------- */
/* Solve the QP from python/experiment_lm/vertexwish.py:
     min sum_e [(x_e - 1/d_i)^2 + (x_e - 1/d_j)^2]   e = {i,j}
     s.t. sum_{e ∋ v} x_e = 1
   via the V×V system (B B^T) λ = B g - 2·1, then x = (g - B^T λ) / 2.
   Returns x in revolutions. */
static int vertexwish_revs(double *x_rev) {
    int V = NV;
    /* g_e = 1/d_i + 1/d_j */
    double *g = (double *)xcalloc(NE, sizeof(double));
    for (int e = 0; e < NE; e++) {
        int a = EDGE_A[e], b = EDGE_B[e];
        g[e] = 1.0 / VERT_DEG[a] + 1.0 / VERT_DEG[b];
    }
    /* (BB^T)_ij = degree if i==j, 1 if {i,j} is an edge, 0 otherwise.
       Build dense V×V (1-indexed; we'll use 0-indexed [0..V-1] mapping to verts 1..V). */
    double *BBt = (double *)xcalloc(V * V, sizeof(double));
    for (int v = 1; v <= V; v++) BBt[(v-1) * V + (v-1)] = (double)VERT_DEG[v];
    for (int e = 0; e < NE; e++) {
        int a = EDGE_A[e] - 1, b = EDGE_B[e] - 1;
        BBt[a * V + b] += 1.0;
        BBt[b * V + a] += 1.0;
    }
    /* rhs = B·g - 2·1.  (B·g)_v = sum over edges incident to v of g_e */
    double *rhs = (double *)xcalloc(V, sizeof(double));
    for (int v = 0; v < V; v++) rhs[v] = -2.0;
    for (int e = 0; e < NE; e++) {
        int a = EDGE_A[e] - 1, b = EDGE_B[e] - 1;
        rhs[a] += g[e];
        rhs[b] += g[e];
    }
    int sing = dense_solve(BBt, rhs, V);
    if (sing) { free(g); free(BBt); free(rhs); return 1; }
    /* x_e = (g_e - λ_a - λ_b) / 2 */
    double max_violation = 0;
    for (int e = 0; e < NE; e++) {
        int a = EDGE_A[e] - 1, b = EDGE_B[e] - 1;
        x_rev[e] = (g[e] - rhs[a] - rhs[b]) / 2.0;
    }
    /* sanity: B x = 1 */
    for (int v = 1; v <= V; v++) {
        double s = 0;
        for (int e = 0; e < NE; e++) {
            if (EDGE_A[e] == v || EDGE_B[e] == v) s += x_rev[e];
        }
        double err = fabs(s - 1.0);
        if (err > max_violation) max_violation = err;
    }
    if (max_violation > 1e-9) {
        fprintf(stderr, "vertexwish: B x = 1 violated by %.3e\n", max_violation);
        free(g); free(BBt); free(rhs); return 2;
    }
    free(g); free(BBt); free(rhs);
    return 0;
}

/* ---------- LM solver ----------------------------------------------------- */
typedef struct {
    bool success;
    int iters;
    double final_resid;
    double final_lambda;
    int total_retries_residual;
    int total_retries_dent;
    const char *msg;
} LMOut;

static double vec_norm(const double *v, int n) {
    double s = 0;
    for (int i = 0; i < n; i++) s += v[i] * v[i];
    return sqrt(s);
}

static LMOut solve_lm(double alpha, double bend_full[],
                       double tol, int max_iter, int max_lm_retries,
                       double lambda_init, double lambda_down, double lambda_up,
                       double lambda_max, double tiny_rel) {
    LMOut out = { false, 0, 0.0, 0.0, 0, 0, "uninit" };
    int M = 3 * N_INT;
    int N = NVAR;

    double *r = (double *)xcalloc(M, sizeof(double));
    double *r_trial = (double *)xcalloc(M, sizeof(double));
    double *r_tmp_a = (double *)xcalloc(M, sizeof(double));
    double *r_tmp_b = (double *)xcalloc(M, sizeof(double));
    double *bend_tmp = (double *)xcalloc(NE, sizeof(double));
    double *bend_trial = (double *)xcalloc(NE, sizeof(double));
    double *J = (double *)xcalloc((size_t)M * (size_t)N, sizeof(double));
    double *A = (double *)xcalloc((size_t)N * (size_t)N, sizeof(double));
    double *g = (double *)xcalloc(N, sizeof(double));
    double *D = (double *)xcalloc(N, sizeof(double));
    double *delta = (double *)xcalloc(N, sizeof(double));
    double *Awork = (double *)xcalloc((size_t)N * (size_t)N, sizeof(double));
    double *bwork = (double *)xcalloc(N, sizeof(double));

    holonomy_residual(alpha, bend_full, r);
    double norm = vec_norm(r, M);
    double lam = lambda_init;

    for (int it = 0; it < max_iter; it++) {
        if (norm <= tol) {
            out.success = true; out.iters = it; out.final_resid = norm;
            out.final_lambda = lam; out.msg = "tol"; goto done;
        }
        /* J = ∂r/∂x analytical (lifted from canonical homotopy_stage). */
        analytical_jacobian(alpha, bend_full, J);
        (void)bend_tmp; (void)r_tmp_a; (void)r_tmp_b;  /* unused now; kept allocated for ABI continuity */
        /* A = J^T J, g = J^T r */
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
        /* D_i = max(A_ii, tiny_rel * max(diag), 1e-30) */
        double max_diag = 0;
        for (int i = 0; i < N; i++) if (A[i * N + i] > max_diag) max_diag = A[i * N + i];
        double floor = tiny_rel * max_diag;
        if (floor < 1e-30) floor = 1e-30;
        for (int i = 0; i < N; i++) {
            D[i] = A[i * N + i] > floor ? A[i * N + i] : floor;
        }

        bool accepted = false;
        for (int rt = 0; rt < max_lm_retries; rt++) {
            /* Awork = A + λ · diag(D); bwork = -g */
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) Awork[i * N + j] = A[i * N + j];
                Awork[i * N + i] += lam * D[i];
                bwork[i] = -g[i];
            }
            int sing = dense_solve(Awork, bwork, N);
            if (sing) {
                lam = lam * lambda_up;
                if (lam > lambda_max) {
                    out.iters = it; out.final_resid = norm; out.final_lambda = lam;
                    out.msg = "lambda_saturated_linsolve"; goto done;
                }
                continue;
            }
            for (int i = 0; i < N; i++) delta[i] = bwork[i];

            /* trial bend */
            for (int e = 0; e < NE; e++) bend_trial[e] = bend_full[e];
            for (int e = 0; e < NE; e++) {
                int j = VAR_OF_E[e];
                if (j >= 0) bend_trial[e] += delta[j];
            }
            holonomy_residual(alpha, bend_trial, r_trial);
            double n_trial = vec_norm(r_trial, M);
            bool dent_bad = has_dent(bend_trial);

            if (n_trial < norm && !dent_bad) {
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
            }
            lam = lam * lambda_up;
            if (lam > lambda_max) {
                out.iters = it; out.final_resid = norm;
                out.msg = "lambda_saturated"; goto done;
            }
        }
        if (!accepted) {
            out.iters = it; out.final_resid = norm;
            out.msg = "lm_retries_exhausted"; goto done;
        }
    }
    out.iters = max_iter; out.final_resid = norm;
    out.success = (norm <= tol);
    out.msg = out.success ? "max_iter_tol" : "max_iter";
done:
    free(r); free(r_trial); free(r_tmp_a); free(r_tmp_b);
    free(bend_tmp); free(bend_trial);
    free(J); free(A); free(g); free(D); free(delta);
    free(Awork); free(bwork);
    return out;
}

/* ---------- main ---------------------------------------------------------- */
static double now_secs(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

int main(int argc, char **argv) {
    tau_const = 2.0 * PI;

    double alpha = PI / 180.0;            /* 1° default */
    double tol = 1e-12;
    int max_iter = 300;
    double lambda_init = 1.0;
    double lambda_down = 0.3, lambda_up = 10.0;
    double lambda_max = 1e12;
    int max_lm_retries = 20;
    double tiny_rel = 1e-12;
    int print_bends = 0;
    const char *bends_out = NULL;
    int start_kind = 0;  /* 0 = vertexwish, 1 = ideal (horou-derived) */

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--alpha") == 0 && i + 1 < argc) {
            alpha = atof(argv[++i]);
            if (!isfinite(alpha)) die("--alpha must be finite");
        }
        else if (strcmp(argv[i], "--alpha-deg") == 0 && i + 1 < argc) {
            double v = atof(argv[++i]);
            if (!isfinite(v)) die("--alpha-deg must be finite");
            alpha = v * PI / 180.0;
        }
        else if (strcmp(argv[i], "--tol") == 0 && i + 1 < argc) {
            tol = atof(argv[++i]);
            if (!isfinite(tol) || tol <= 0) die("--tol must be finite > 0");
        }
        else if (strcmp(argv[i], "--max-iter") == 0 && i + 1 < argc) {
            max_iter = atoi(argv[++i]);
            if (max_iter <= 0) die("--max-iter must be > 0");
        }
        else if (strcmp(argv[i], "--lambda-init") == 0 && i + 1 < argc) {
            lambda_init = atof(argv[++i]);
            if (!isfinite(lambda_init) || lambda_init <= 0) die("--lambda-init must be finite > 0");
        }
        else if (strcmp(argv[i], "--lambda-down") == 0 && i + 1 < argc) {
            lambda_down = atof(argv[++i]);
            if (!isfinite(lambda_down) || lambda_down <= 0 || lambda_down >= 1)
                die("--lambda-down must be in (0,1)");
        }
        else if (strcmp(argv[i], "--lambda-up") == 0 && i + 1 < argc) {
            lambda_up = atof(argv[++i]);
            if (!isfinite(lambda_up) || lambda_up <= 1) die("--lambda-up must be > 1");
        }
        else if (strcmp(argv[i], "--lambda-max") == 0 && i + 1 < argc) {
            lambda_max = atof(argv[++i]);
            if (!isfinite(lambda_max) || lambda_max <= 0) die("--lambda-max must be finite > 0");
        }
        else if (strcmp(argv[i], "--max-lm-retries") == 0 && i + 1 < argc) {
            max_lm_retries = atoi(argv[++i]);
            if (max_lm_retries <= 0) die("--max-lm-retries must be > 0");
        }
        else if (strcmp(argv[i], "--print-bends") == 0) print_bends = 1;
        else if (strcmp(argv[i], "--bends-out") == 0 && i + 1 < argc) {
            bends_out = argv[++i];
        }
        else if (strcmp(argv[i], "--start") == 0 && i + 1 < argc) {
            const char *s = argv[++i];
            if      (strcmp(s, "vertexwish") == 0) start_kind = 0;
            else if (strcmp(s, "ideal") == 0)      start_kind = 1;
            else die("--start must be 'vertexwish' or 'ideal'");
        }
        else die("unknown command-line option");
    }

    parse_netcode_stream(stdin);
    build_topology();
    choose_base_face();

    double *bend = (double *)xcalloc(NE, sizeof(double));
    if (start_kind == 0) {
        double *x_rev = (double *)xcalloc(NE, sizeof(double));
        int sw = vertexwish_revs(x_rev);
        if (sw) { free(x_rev); die("vertexwish failed"); }
        for (int e = 0; e < NE; e++) bend[e] = tau_const * x_rev[e];
        free(x_rev);
    } else {
        double *u = (double *)xcalloc(NV + 1, sizeof(double));
        if (horou_solve(u) < 0) { free(u); die("horou failed"); }
        if (compute_bends_at_zero(u, bend) < 0) { free(u); die("compute_bends_at_zero failed"); }
        free(u);
    }

    double t0 = now_secs();
    LMOut out = solve_lm(alpha, bend, tol, max_iter, max_lm_retries,
                         lambda_init, lambda_down, lambda_up, lambda_max, tiny_rel);
    double dt = now_secs() - t0;

    printf("nv: %d\n", NV);
    printf("ne: %d\n", NE);
    printf("nvar: %d\n", NVAR);
    printf("base_face: %d %d %d\n", BASE[0], BASE[1], BASE[2]);
    printf("alpha_rad: %.17g\n", alpha);
    printf("alpha_deg: %.6f\n", alpha * 180.0 / PI);
    printf("lambda_init: %.6e\n", lambda_init);
    printf("max_iter: %d\n", max_iter);
    printf("tol: %.6e\n", tol);
    printf("solver: lm\n");
    printf("start: %s\n", start_kind == 0 ? "vertexwish" : "ideal");
    printf("success: %s\n", out.success ? "true" : "false");
    printf("iters: %d\n", out.iters);
    printf("final_resid: %.6e\n", out.final_resid);
    printf("final_lambda: %.6e\n", out.final_lambda);
    printf("retries_residual: %d\n", out.total_retries_residual);
    printf("retries_dent: %d\n", out.total_retries_dent);
    printf("message: %s\n", out.msg);
    printf("wall_secs: %.6f\n", dt);
    if (print_bends) {
        for (int e = 0; e < NE; e++) {
            printf("bend %d-%d %.17g\n", EDGE_A[e], EDGE_B[e], bend[e]);
        }
    }

    /* puffup-bends 1 writer — matches canonical homotopy_stage format
     * so euclid_realize can consume the file. Emits faces in input
     * order. Bends must be in canonical homotopy_stage edge order
     * (`for u=1..NV walk NBR[u] in nbr_add insertion order, emit (u,w)
     * iff u<w`). puffup_c_lm.c's EDGE_A/EDGE_B array is in face-scan
     * encounter order, which differs in general — so we rebuild the
     * canonical order locally. */
    if (bends_out) {
        FILE *fh = fopen(bends_out, "w");
        if (!fh) {
            fprintf(stderr, "ERROR: cannot open --bends-out %s\n", bends_out);
        } else {
            /* canonical-order rebuild */
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
            #undef CANON_ADD

            fprintf(fh, "puffup-bends 1\n");
            fprintf(fh, "NV %d NE %d alpha_deg %.15g\n",
                    NV, NE, alpha * 180.0 / PI);
            fprintf(fh, "faces ");
            for (int i = 0; i < NF; i++) {
                fprintf(fh, "%d,%d,%d%s",
                        FACES[i][0], FACES[i][1], FACES[i][2],
                        i + 1 < NF ? ";" : "");
            }
            fprintf(fh, "\n");
            fprintf(fh, "bends\n");
            int rows = 0;
            for (int u = 1; u <= NV; u++) {
                for (int j = 0; j < canon_nnbr[u]; j++) {
                    int w = canon_nbr[u][j];
                    if (u < w) {
                        int e = EDGE_IDX_LOOKUP[u][w];
                        if (e < 0) {
                            fprintf(stderr,
                                "ERROR: --bends-out: edge (%d,%d) has no index\n",
                                u, w);
                            fclose(fh);
                            free(bend);
                            return 1;
                        }
                        fprintf(fh, "%d %d %.15g\n", u, w, bend[e]);
                        rows++;
                    }
                }
            }
            if (rows != NE) {
                fprintf(stderr,
                    "ERROR: --bends-out emitted %d rows, expected NE=%d\n",
                    rows, NE);
                fclose(fh);
                free(bend);
                return 1;
            }
            fclose(fh);
        }
    }

    free(bend);
    return out.success ? 0 : 1;
}
