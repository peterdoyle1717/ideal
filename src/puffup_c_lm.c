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
   may alias caller's residual buffer. */
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
        /* J = ∂r/∂x via FD. r_tmp_a/b are scratch; do not alias r. */
        fd_jacobian(alpha, bend_full, J, bend_tmp, r_tmp_a, r_tmp_b);
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
        else die("unknown command-line option");
    }

    parse_netcode_stream(stdin);
    build_topology();
    choose_base_face();

    /* vertexwish start in revolutions */
    double *x_rev = (double *)calloc(NE, sizeof(double));
    int sw = vertexwish_revs(x_rev);
    if (sw) { free(x_rev); die("vertexwish failed"); }

    /* convert to radians */
    double *bend = (double *)calloc(NE, sizeof(double));
    for (int e = 0; e < NE; e++) bend[e] = tau_const * x_rev[e];
    free(x_rev);

    double t0 = now_secs();
    LMOut out = solve_lm(alpha, bend, tol, max_iter, max_lm_retries,
                         lambda_init, lambda_down, lambda_up, lambda_max, tiny_rel);
    double dt = now_secs() - t0;

    printf("nv: %d\n", NV);
    printf("ne: %d\n", NE);
    printf("nvar: %d\n", NVAR);
    printf("alpha_rad: %.17g\n", alpha);
    printf("alpha_deg: %.6f\n", alpha * 180.0 / PI);
    printf("lambda_init: %.6e\n", lambda_init);
    printf("max_iter: %d\n", max_iter);
    printf("tol: %.6e\n", tol);
    printf("solver: lm\n");
    printf("start: vertexwish\n");
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

    free(bend);
    return out.success ? 0 : 1;
}
