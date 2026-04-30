/*
 * puffup.c — neoplatonic dihedral solver via halfway-to-60° homotopy.
 *
 * Pipeline per net:
 *   parse facelist  →  horou (u-values)  →  dihedral (α=0 bends)
 *                  →  homotopy (α=0 → π/3) with halve-α + halfway cap.
 *
 * Input:  facelist on stdin, one net per line ("a,b,c;d,e,f;..." 1-indexed).
 *         Vertex 1 is the "point at infinity" in horou's convention.
 * Output: one stats line per input:
 *           STATUS V steps halvings newton_iters final_alpha_deg
 *         STATUS in {ok, dented, max_iter, lstsq, horou, parse}.
 *
 * Compile:  cc -O3 -o puffup puffup.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXV    1200      /* room for icosa{5,2}=392 + headroom for bigger subdivisions */
#define MAXF    (2*MAXV + 4)
#define MAXE    (3*MAXV - 6)
#define MAXRING 16
#define MAXLINE 65536     /* 16K was tight for V > 400 face lists */

/* number of homotopy variables = non-base bends = E - 3 = 3V-9 */
#define MAXN    (3*MAXV - 9)

/* ---------- graph state ---------------------------------------------------- */
typedef struct { int a, b, c; } Face;

static int  NV, NF;
static Face FACES[MAXF];
static int  EM[MAXV+1][MAXV+1];        /* third vertex of face on directed (a,b) */
static int  EM_F[MAXV+1][MAXV+1];      /* face index for directed (a,b) */
static int  DEG[MAXV+1];
static int  NBR[MAXV+1][MAXRING];
static int  NNBR[MAXV+1];
static short DU[MAXF*3], DW[MAXF*3];
static int  ND;

/* edges (canonical, EDGE_A[i] < EDGE_B[i]) */
static int  NE;
static int  EDGE_A[MAXE], EDGE_B[MAXE];
static int  EDGE_IDX[MAXV+1][MAXV+1];  /* edge index for unordered (u,v); -1 if none */

/* per-vertex flower of face indices (cyclic walk) */
static int  FLOWER[MAXV+1][MAXRING];
static int  FLOWER_LEN[MAXV+1];

/* ---------- input parser --------------------------------------------------- */
/* Returns NF on success, 0 on parse failure, -1 on size overflow
   (NF would exceed MAXF, or any vertex index exceeds MAXV).
   On -1 the caller should report size_too_big and skip; FACES[] is
   NOT written past its bound. */
static int parse_facelist(const char *line) {
    NF = 0; NV = 0;
    const char *p = line;
    while (*p) {
        int a, b, c;
        if (sscanf(p, "%d,%d,%d", &a, &b, &c) != 3) break;
        if (NF >= MAXF || a > MAXV || b > MAXV || c > MAXV) {
            if (a > NV) NV = a;
            if (b > NV) NV = b;
            if (c > NV) NV = c;
            return -1;
        }
        FACES[NF++] = (Face){a, b, c};
        if (a > NV) NV = a;
        if (b > NV) NV = b;
        if (c > NV) NV = c;
        while (*p && *p != ';') p++;
        if (*p == ';') p++;
    }
    return NF;
}

static int nbr_add(int u, int w) {
    for (int i = 0; i < NNBR[u]; i++) if (NBR[u][i] == w) return 0;
    if (NNBR[u] >= MAXRING) {
        fprintf(stderr, "ERROR: vertex %d degree exceeds MAXRING=%d "
                        "(input is not a valid 6-net)\n", u, MAXRING);
        return -1;
    }
    NBR[u][NNBR[u]++] = w;
    return 0;
}

static void build_clear(void) {
    for (int i = 0; i < ND; i++) {
        EM  [DU[i]][DW[i]] = 0;
        EM_F[DU[i]][DW[i]] = 0;
    }
    ND = 0;
    for (int i = 0; i < NE; i++) {
        EDGE_IDX[EDGE_A[i]][EDGE_B[i]] = -1;
        EDGE_IDX[EDGE_B[i]][EDGE_A[i]] = -1;
    }
    NE = 0;
}

static int build(void) {
    memset(DEG,  0, (NV+2)*sizeof(int));
    memset(NNBR, 0, (NV+2)*sizeof(int));
    ND = 0;
    NE = 0;
    /* directed edges + neighbors */
    for (int i = 0; i < NF; i++) {
        int a=FACES[i].a, b=FACES[i].b, c=FACES[i].c;
        EM[a][b]=c;  EM_F[a][b]=i;  DU[ND]=a; DW[ND]=b; ND++;
        EM[b][c]=a;  EM_F[b][c]=i;  DU[ND]=b; DW[ND]=c; ND++;
        EM[c][a]=b;  EM_F[c][a]=i;  DU[ND]=c; DW[ND]=a; ND++;
        DEG[a]++; DEG[b]++; DEG[c]++;
        if (nbr_add(a,b) || nbr_add(b,a) ||
            nbr_add(b,c) || nbr_add(c,b) ||
            nbr_add(a,c) || nbr_add(c,a)) return -1;
    }
    /* canonical undirected edges */
    for (int u = 1; u <= NV; u++) {
        for (int j = 0; j < NNBR[u]; j++) {
            int w = NBR[u][j];
            if (u < w) {
                EDGE_A[NE] = u; EDGE_B[NE] = w;
                EDGE_IDX[u][w] = NE;
                EDGE_IDX[w][u] = NE;
                NE++;
            }
        }
    }
    /* flowers: walk faces around each vertex, matching plat1000 / puffup.py:
       from face (v,b,c) step to face containing directed edge (v,c). */
    for (int v = 1; v <= NV; v++) {
        int start = -1;
        for (int i = 0; i < NF; i++) {
            if (FACES[i].a==v || FACES[i].b==v || FACES[i].c==v) { start = i; break; }
        }
        if (start < 0) { FLOWER_LEN[v] = 0; continue; }
        int cur = start, k = 0;
        for (;;) {
            if (k >= MAXRING) {
                fprintf(stderr, "ERROR: vertex %d flower exceeds MAXRING=%d "
                                "(input is not a valid 6-net)\n", v, MAXRING);
                return -1;
            }
            FLOWER[v][k++] = cur;
            int a=FACES[cur].a, b=FACES[cur].b, c=FACES[cur].c;
            int third;
            if      (v==a) third = c;
            else if (v==b) third = a;
            else           third = b;
            int nxt = EM_F[v][third];
            if (nxt == start) break;
            cur = nxt;
        }
        FLOWER_LEN[v] = k;
    }
    return 0;
}

static int cyclic_nbrs(int v, int ring[]) {
    if (!NNBR[v]) return 0;
    int start = NBR[v][0];
    ring[0] = start;
    int k = 1, cur = start;
    for (;;) {
        int nxt = EM[v][cur];
        if (nxt == start) break;
        if (k >= MAXRING) {
            fprintf(stderr, "ERROR: vertex %d cyclic-nbrs exceeds MAXRING=%d\n",
                    v, MAXRING);
            return -1;
        }
        ring[k++] = nxt;
        cur = nxt;
    }
    return k;
}

/* ---------- horou (port of ideal/src/horou_c.c) --------------------------- */
static double petal(double ui, double uj, double uk) {
    double a=ui*uj, b=ui*uk, c=uj*uk;
    double cos_t = (a*a + b*b - c*c) / (2.0*a*b);
    if (cos_t >  1.0) cos_t =  1.0;
    if (cos_t < -1.0) cos_t = -1.0;
    return acos(cos_t);
}
static void petal_grad(double ui, double uj, double uk,
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

/* generic in-place LU + back-sub. row-major n×n in J, RHS in b. */
static int lu_solve(double *J, double *b, int n, int stride) {
    #define A(i,j) J[(i)*stride + (j)]
    for (int col = 0; col < n; col++) {
        int piv = col;
        double best = fabs(A(col,col));
        for (int r = col+1; r < n; r++) {
            double v = fabs(A(r,col));
            if (v > best) { best = v; piv = r; }
        }
        if (best < 1e-14) return -1;
        if (piv != col) {
            for (int k = col; k < n; k++) { double t=A(col,k); A(col,k)=A(piv,k); A(piv,k)=t; }
            { double t=b[col]; b[col]=b[piv]; b[piv]=t; }
        }
        double inv = 1.0 / A(col,col);
        for (int r = col+1; r < n; r++) {
            double fac = A(r,col) * inv;
            for (int k = col; k < n; k++) A(r,k) -= fac * A(col,k);
            b[r] -= fac * b[col];
        }
    }
    for (int i = n-1; i >= 0; i--) {
        double s = b[i];
        for (int j = i+1; j < n; j++) s -= A(i,j) * b[j];
        b[i] = s / A(i,i);
    }
    return 0;
    #undef A
}

/* horou solver: u_out[v-1] = u[v]; u_out[0] = NaN (vertex 1 = ∞).
   Returns 0 on success, -1 on Newton failure. */
static double HJ[MAXV][MAXV];
static double HF[MAXV];
static double Hdx[MAXV];

static int horou(double u_out[]) {
    static int  bndry[MAXV+1];
    static int  int_idx[MAXV+1];
    static int  interior[MAXV];
    static int  ring[MAXV+1][MAXRING];
    static int  ringlen[MAXV+1];
    static int  ff_a[MAXF], ff_b[MAXF], ff_c[MAXF];
    static double xvec[MAXV];

    const double TAU = 2.0 * M_PI;
    const double target = TAU;

    memset(bndry, 0, (NV+2)*sizeof(int));
    memset(int_idx, -1, (NV+2)*sizeof(int));
    int br[MAXV+1], n_b = cyclic_nbrs(1, br);
    if (n_b < 0) return -1;
    for (int i = 0; i < n_b; i++) bndry[br[i]] = 1;

    int n_int = 0;
    for (int v = 2; v <= NV; v++) {
        if (!bndry[v]) { int_idx[v] = n_int; interior[n_int++] = v; }
    }
    #define U(vv) (bndry[vv] ? 1.0 : xvec[int_idx[vv]])

    for (int i = 0; i < n_int; i++) {
        int v = interior[i];
        ringlen[v] = cyclic_nbrs(v, ring[v]);
        if (ringlen[v] < 0) return -1;
    }
    int nff = 0;
    for (int i = 0; i < NF; i++) {
        int a=FACES[i].a, b=FACES[i].b, c=FACES[i].c;
        if (a!=1 && b!=1 && c!=1) { ff_a[nff]=a; ff_b[nff]=b; ff_c[nff]=c; nff++; }
    }

    if (n_int == 0) goto done;

    for (int i = 0; i < n_int; i++) xvec[i] = 1.0;

    for (int it = 0; it < 200; it++) {
        double res = 0.0;
        for (int i = 0; i < n_int; i++) {
            int v = interior[i]; int k = ringlen[v]; double ui = xvec[i];
            double s = 0.0;
            for (int j = 0; j < k; j++)
                s += petal(ui, U(ring[v][j]), U(ring[v][(j+1)%k]));
            HF[i] = s - target;
            double af = fabs(HF[i]); if (af > res) res = af;
        }
        if (res < 1e-12) break;

        for (int i = 0; i < n_int; i++)
            for (int j = 0; j < n_int; j++) HJ[i][j] = 0.0;
        for (int i = 0; i < n_int; i++) {
            int v = interior[i]; int k = ringlen[v]; double ui = xvec[i];
            for (int j = 0; j < k; j++) {
                int vj = ring[v][j], vk = ring[v][(j+1)%k];
                double uj = U(vj), uk = U(vk);
                double dui, duj, duk;
                petal_grad(ui, uj, uk, &dui, &duj, &duk);
                HJ[i][i] += dui;
                if (int_idx[vj] >= 0) HJ[i][int_idx[vj]] += duj;
                if (int_idx[vk] >= 0) HJ[i][int_idx[vk]] += duk;
            }
        }
        for (int i = 0; i < n_int; i++) Hdx[i] = -HF[i];
        if (lu_solve(&HJ[0][0], Hdx, n_int, MAXV) < 0) return -1;

        double step = 1.0;
        int found = 0;
        for (int bt = 0; bt < 60; bt++, step *= 0.5) {
            int ok = 1;
            for (int i = 0; i < n_int; i++) {
                if (xvec[i] + step*Hdx[i] <= 0.0) { ok=0; break; }
            }
            if (!ok) continue;
            for (int f = 0; f < nff && ok; f++) {
                #define UU(vv) (int_idx[vv] >= 0 \
                    ? xvec[int_idx[vv]] + step*Hdx[int_idx[vv]] : 1.0)
                double ua=UU(ff_a[f]), ub=UU(ff_b[f]), uc=UU(ff_c[f]);
                double p=ua*ub, q=ua*uc, r=ub*uc;
                if (p+q<=r || p+r<=q || q+r<=p) ok=0;
                #undef UU
            }
            if (!ok) continue;
            double res2 = 0.0;
            for (int i = 0; i < n_int; i++) {
                int v=interior[i]; int k=ringlen[v];
                double ui = xvec[i] + step*Hdx[i];
                double s = 0.0;
                for (int j = 0; j < k; j++) {
                    int rj=ring[v][j], rk=ring[v][(j+1)%k];
                    double uj = (int_idx[rj] >= 0 ? xvec[int_idx[rj]] + step*Hdx[int_idx[rj]] : 1.0);
                    double uk2= (int_idx[rk] >= 0 ? xvec[int_idx[rk]] + step*Hdx[int_idx[rk]] : 1.0);
                    s += petal(ui, uj, uk2);
                }
                double af = fabs(s - target); if (af > res2) res2 = af;
            }
            if (res2 < res) { found = 1; break; }
        }
        if (!found) break;
        for (int i = 0; i < n_int; i++) xvec[i] += step * Hdx[i];
    }

done:
    u_out[0] = NAN;  /* vertex 1 = ∞ */
    for (int v = 2; v <= NV; v++) {
        u_out[v-1] = bndry[v] ? 1.0 : xvec[int_idx[v]];
    }
    #undef U
    return 0;
}

/* ---------- α=0 dihedral per edge from u ---------------------------------- */
/* Returns angle at v in face (v,p,q) in horou metric. */
static double angleat(int v, int p, int q, const double u[]) {
    double x = 1.0/u[v-1], y = 1.0/u[p-1], z = 1.0/u[q-1];
    double cos_t = (y*y + z*z - x*x) / (2.0*y*z);
    if (cos_t >  1.0) cos_t =  1.0;
    if (cos_t < -1.0) cos_t = -1.0;
    return acos(cos_t);
}

/* Sum of petals at v over faces not touching vertex 1. */
static double boundaryangleat(int v, const double u[]) {
    double tot = 0.0;
    for (int t = 0; t < FLOWER_LEN[v]; t++) {
        int fi = FLOWER[v][t];
        int a=FACES[fi].a, b=FACES[fi].b, c=FACES[fi].c;
        if (a==1 || b==1 || c==1) continue;
        int p, q;
        if      (v==a) { p=b; q=c; }
        else if (v==b) { p=c; q=a; }
        else           { p=a; q=b; }
        tot += angleat(v, p, q, u);
    }
    return tot;
}

/* For each edge i (canonical (a,b)), set BEND[i] = α=0 bend (= π − interior dihedral). */
static int compute_bends_at_zero(const double u[], double bend_out[]) {
    for (int i = 0; i < NE; i++) {
        int a = EDGE_A[i], b = EDGE_B[i];
        /* the two faces adjacent to this edge: EM_F[a][b] and EM_F[b][a] */
        int f1 = EM_F[a][b], f2 = EM_F[b][a];
        int third1, third2;
        {
            int x=FACES[f1].a, y=FACES[f1].b, z=FACES[f1].c;
            third1 = (x!=a && x!=b) ? x : ((y!=a && y!=b) ? y : z);
        }
        {
            int x=FACES[f2].a, y=FACES[f2].b, z=FACES[f2].c;
            third2 = (x!=a && x!=b) ? x : ((y!=a && y!=b) ? y : z);
        }
        int touches_inf_1 = (FACES[f1].a==1 || FACES[f1].b==1 || FACES[f1].c==1);
        int touches_inf_2 = (FACES[f2].a==1 || FACES[f2].b==1 || FACES[f2].c==1);

        double bend;
        if (a == 1 || b == 1) {
            int v = (a == 1) ? b : a;
            bend = M_PI - boundaryangleat(v, u);
        } else if (!touches_inf_1 && !touches_inf_2) {
            double pa = angleat(third1, a, b, u);
            double pd = angleat(third2, a, b, u);
            bend = M_PI - pa - pd;
        } else {
            int finite_face = touches_inf_1 ? f2 : f1;
            int finite_third = touches_inf_1 ? third2 : third1;
            (void)finite_face;
            bend = M_PI - angleat(finite_third, a, b, u);
        }
        bend_out[i] = bend;
    }
    return 0;
}

/* ---------- holonomy machinery -------------------------------------------- */
typedef double M3[3][3];

static inline void matmul(const M3 A, const M3 B, M3 C) {
    M3 T;
    for (int i=0;i<3;i++)
        for (int j=0;j<3;j++) {
            double s=0;
            for (int k=0;k<3;k++) s += A[i][k]*B[k][j];
            T[i][j] = s;
        }
    memcpy(C, T, sizeof(M3));
}
static inline void mat_eye(M3 A) {
    for (int i=0;i<3;i++) for (int j=0;j<3;j++) A[i][j] = (i==j);
}
static inline void mat_copy(const M3 A, M3 B) { memcpy(B, A, sizeof(M3)); }

/* movemat(α, β) = matz(α) · matx(−β). */
static inline void movemat(double alpha, double beta, M3 M) {
    double ca=cos(alpha), sa=sin(alpha);
    double cb=cos(beta),  sb=sin(beta);   /* matx(-β) → cos(-β)=cb, sin(-β)=-sb */
    /* matz(α) = [[ca,-sa,0],[sa,ca,0],[0,0,1]]
       matx(-β) = [[1,0,0],[0,cb,sb],[0,-sb,cb]] */
    M[0][0] = ca;
    M[0][1] = -sa*cb;
    M[0][2] = -sa*sb;
    M[1][0] = sa;
    M[1][1] = ca*cb;
    M[1][2] = ca*sb;
    M[2][0] = 0.0;
    M[2][1] = -sb;
    M[2][2] = cb;
}

/* dM/dβ at given α,β: -matz(α) · Jx · matx(-β), where Jx generates x-rotations */
static inline void dmovemat_dbeta(double alpha, double beta, M3 D) {
    /* Jx · matx(-β) = [[0,0,0],[0,sb,-cb],[0,cb,sb]]   (because Jx swaps rows 1,2 with sign)
       More directly: dM/dβ = column-wise derivatives; we computed
       d/dβ matx(-β) = -Jx · matx(-β) = [[0,0,0],[0,-sb,cb],[0,-cb,-sb]].
       Then -matz(α) · Jx · matx(-β) = matz(α) · [[0,0,0],[0,-sb,cb],[0,-cb,-sb]] */
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

/* ---------- base face / vars ---------------------------------------------- */
static int  N_INT;             /* non-base vertex count = NV - 3 */
static int  INT_VS[MAXV];      /* the non-base vertices */
static int  INT_OF[MAXV+1];    /* vertex → INT index, -1 if base */
static int  BASE_VS[3];
static int  BASE_EDGES[3];
static int  IS_BASE_EDGE[MAXE];

static int  NVAR;
static int  VAR_OF_EDGE[MAXE]; /* edge_idx → var_idx, -1 if base edge */

/* For each non-base vertex, precompute the sequence of edge indices in flower order.
   FLOWER_E[v][t] is the edge index of the t-th edge in v's flower (the edge from v to
   the third vertex of FLOWER[v][t]). */
static int FLOWER_E[MAXV+1][MAXRING];

static void choose_base_face(void) {
    /* CLERS-decoded face[0] = (1,2,3) — use it as the canonical base. */
    BASE_VS[0] = FACES[0].a;
    BASE_VS[1] = FACES[0].b;
    BASE_VS[2] = FACES[0].c;
    /* base edges */
    for (int k = 0; k < 3; k++) {
        int u = BASE_VS[k], v = BASE_VS[(k+1)%3];
        BASE_EDGES[k] = EDGE_IDX[u][v];
    }
    memset(IS_BASE_EDGE, 0, NE*sizeof(int));
    IS_BASE_EDGE[BASE_EDGES[0]] = 1;
    IS_BASE_EDGE[BASE_EDGES[1]] = 1;
    IS_BASE_EDGE[BASE_EDGES[2]] = 1;

    /* interior (non-base) vertices */
    N_INT = 0;
    for (int v = 1; v <= NV; v++) INT_OF[v] = -1;
    for (int v = 1; v <= NV; v++) {
        if (v != BASE_VS[0] && v != BASE_VS[1] && v != BASE_VS[2]) {
            INT_OF[v] = N_INT;
            INT_VS[N_INT++] = v;
        }
    }

    /* var edges */
    NVAR = 0;
    for (int i = 0; i < NE; i++) {
        if (IS_BASE_EDGE[i]) { VAR_OF_EDGE[i] = -1; }
        else { VAR_OF_EDGE[i] = NVAR++; }
    }

    /* per-vertex flower-edge lookup */
    for (int v = 1; v <= NV; v++) {
        int k = FLOWER_LEN[v];
        for (int t = 0; t < k; t++) {
            int fi = FLOWER[v][t];
            int a=FACES[fi].a, b=FACES[fi].b, c=FACES[fi].c;
            int third;
            if      (v==a) third = c;
            else if (v==b) third = a;
            else           third = b;
            FLOWER_E[v][t] = EDGE_IDX[v][third];
        }
    }
}

/* ---------- residual + analytical Jacobian -------------------------------- */
/* Compute holomat at vertex v: M = M_0 · M_1 · ... · M_{k-1}.
   Returns 3 off-diagonals (M[0,1], M[0,2], M[1,2]) into out[3]. */
static void vertex_residual(int v, double alpha, const double bend[],
                             double out[3]) {
    int k = FLOWER_LEN[v];
    M3 M; mat_eye(M);
    for (int t = 0; t < k; t++) {
        int e = FLOWER_E[v][t];
        M3 step;
        movemat(alpha, bend[e], step);
        matmul(M, step, M);
    }
    out[0] = M[0][1];
    out[1] = M[0][2];
    out[2] = M[1][2];
}

/* Build the full residual r[3*N_INT] = stacked vertex residuals. */
static void holonomy_residual(double alpha, const double bend[], double r[]) {
    for (int i = 0; i < N_INT; i++) {
        vertex_residual(INT_VS[i], alpha, bend, &r[3*i]);
    }
}

/* Build analytical Jacobian J[3*N_INT][NVAR] (row-major, stride NVAR). */
static double Jbuf[3*MAXV * MAXN];

static void analytical_jacobian(double alpha, const double bend[]) {
    int rows = 3 * N_INT;
    memset(Jbuf, 0, rows * NVAR * sizeof(double));

    for (int i = 0; i < N_INT; i++) {
        int v = INT_VS[i];
        int k = FLOWER_LEN[v];
        /* per-step movemats */
        M3 Ms[MAXRING];
        for (int t = 0; t < k; t++) {
            movemat(alpha, bend[FLOWER_E[v][t]], Ms[t]);
        }
        /* prefix products P[0]=I, P[t] = M_0 ... M_{t-1} */
        M3 P[MAXRING+1];
        mat_eye(P[0]);
        for (int t = 0; t < k; t++) matmul(P[t], Ms[t], P[t+1]);
        /* suffix products S[k]=I, S[t] = M_t ... M_{k-1} */
        M3 S[MAXRING+1];
        mat_eye(S[k]);
        for (int t = k-1; t >= 0; t--) matmul(Ms[t], S[t+1], S[t]);

        for (int t = 0; t < k; t++) {
            int e = FLOWER_E[v][t];
            int col = VAR_OF_EDGE[e];
            if (col < 0) continue;  /* base edge — skip */
            M3 dM, tmp, contrib;
            dmovemat_dbeta(alpha, bend[e], dM);
            matmul(P[t], dM, tmp);
            matmul(tmp, S[t+1], contrib);
            int row = 3 * i;
            Jbuf[(row+0)*NVAR + col] += contrib[0][1];
            Jbuf[(row+1)*NVAR + col] += contrib[0][2];
            Jbuf[(row+2)*NVAR + col] += contrib[1][2];
        }
    }
}

/* ---------- Newton with dent check ---------------------------------------- */
/* dent at vertex v: sum of bends on incident edges (in flower order) < 0. */
static int has_dent(const double bend[]) {
    /* for each vertex sum bends on its flower edges */
    for (int v = 1; v <= NV; v++) {
        int k = FLOWER_LEN[v];
        double s = 0.0;
        for (int t = 0; t < k; t++) s += bend[FLOWER_E[v][t]];
        if (s < 0.0) return 1;
    }
    return 0;
}

/* Try to solve at α from bend_in. Updates bend_out (with var-edges replaced).
   Returns: 0=ok, 1=dented, 2=lstsq_fail, 3=max_iter. */
static double rvec[3*MAXV];
static double dxvec[MAXN];

static int try_full_newton(double alpha, const double bend_in[],
                            double bend_out[], int max_iter, double tol,
                            int *iters_used) {
    /* Copy bend dict; only var entries change. */
    memcpy(bend_out, bend_in, NE*sizeof(double));

    for (int it = 0; it < max_iter; it++) {
        holonomy_residual(alpha, bend_out, rvec);
        double norm2 = 0.0;
        for (int i = 0; i < 3*N_INT; i++) norm2 += rvec[i]*rvec[i];
        if (sqrt(norm2) <= tol) { *iters_used = it; return 0; }

        analytical_jacobian(alpha, bend_out);
        for (int i = 0; i < 3*N_INT; i++) dxvec[i] = -rvec[i];
        if (lu_solve(Jbuf, dxvec, NVAR, NVAR) < 0) {
            *iters_used = it + 1;
            return 2;
        }
        /* full step */
        for (int i = 0; i < NE; i++) {
            if (VAR_OF_EDGE[i] >= 0) {
                bend_out[i] += dxvec[VAR_OF_EDGE[i]];
            }
        }
        /* dent check on new x */
        if (has_dent(bend_out)) { *iters_used = it + 1; return 1; }
    }
    *iters_used = max_iter;
    return 3;
}

/* ---------- complete base bends (9→3 lstsq via normal equations) ----------
 *
 * After homotopy, the 3 base-edge bends are still at their α=0 values; they
 * satisfy the dent check (sums) but not the base-vertex holonomy. Solve for
 * them by FD-Newton on the 9 base-vertex residuals. */
/* Closed-form base bend imputation via PYP factorization.
 *
 * For each base vertex BASE_VS[bi], the puffup holonomy product is
 *   M_v = Y(α) P(b_first) · K · Y(α) P(b_last) = I,
 * where K is the product of the d-2 non-base movemats in v's flower.
 * Letting X = K · Y(α), algebra gives
 *   X^T = P(b_last) · Y(α) · P(b_first).
 * The user's PYP-factor formula then reads off
 *   b_last = atan2(-X^T[2,0], X^T[1,0]) = atan2(-X[0,2], X[0,1]),
 *   b_first = atan2(-X^T[0,2], -X^T[0,1]) = atan2(-X[2,0], -X[1,0])  (debug check).
 * b_last is the bend on the edge from BASE_VS[bi] to BASE_VS[(bi+1)%3].
 * Compute b_last freshly at each base vertex; b_first is the same value
 * recovered by the previous vertex and serves as a consistency check. */
static int complete_base_bends(double alpha, double bend[]) {
    M3 Yalpha;
    {
        double ca=cos(alpha), sa=sin(alpha);
        Yalpha[0][0]= ca; Yalpha[0][1]=-sa; Yalpha[0][2]=0;
        Yalpha[1][0]= sa; Yalpha[1][1]= ca; Yalpha[1][2]=0;
        Yalpha[2][0]=  0; Yalpha[2][1]=  0; Yalpha[2][2]=1;
    }
    double recovered_A[3], recovered_B[3];
    for (int bi = 0; bi < 3; bi++) {
        int v = BASE_VS[bi];
        int k = FLOWER_LEN[v];
        if (k < 3) return -1;     /* deg too low to have non-base flower */
        M3 K; mat_eye(K);
        for (int t = 1; t < k - 1; t++) {
            int e = FLOWER_E[v][t];
            M3 step;
            movemat(alpha, bend[e], step);
            matmul(K, step, K);
        }
        M3 X; matmul(K, Yalpha, X);
        recovered_A[bi] = atan2(-X[0][2],  X[0][1]);   /* = b on outgoing base edge */
        recovered_B[bi] = atan2(-X[2][0], -X[1][0]);   /* = b on incoming base edge */
    }
    /* Write recovered_A[bi] to bend on edge (BASE_VS[bi], BASE_VS[(bi+1)%3]). */
    for (int bi = 0; bi < 3; bi++) {
        int v_next = BASE_VS[(bi + 1) % 3];
        int e_idx = EDGE_IDX[BASE_VS[bi]][v_next];
        bend[e_idx] = recovered_A[bi];
    }
    return 0;
}

/* ---------- halfway homotopy --------------------------------------------- */
typedef struct {
    int   status;       /* 0=ok, 1=stuck */
    int   n_steps, n_halve, newton_iters;
    double final_alpha;
} HomotopyResult;

static double bends_curr[MAXE];
static double bends_trial[MAXE];

/* Default homotopy params (override via CLI flags from main). */
static double CFG_LOOSE_TOL           = 1e-3;
static int    CFG_MAX_NEWTON          = 8;
static double CFG_ALPHA_STEP_INIT_DEG = 1.0;

/* Retry params if first pass fails. Set via --retry-tol etc. */
static double CFG_RETRY_TOL           = 1e-8;
static int    CFG_RETRY_MAX_NEWTON    = 50;
static double CFG_RETRY_INIT_DEG      = 0.25;

static HomotopyResult homotopy_with(const double bends_init[],
                                     double init_step_deg, int max_newton_iters,
                                     double loose_tol) {
    const double GOAL = M_PI / 3.0;
    const double ALPHA_FINAL = GOAL * (1.0 - 1e-12);
    const double ALPHA_STEP_INIT = init_step_deg * M_PI / 180.0;
    const double MIN_ALPHA_STEP = 1e-12;
    const int    MAX_ATTEMPTS = 8000;             /* large enough for retry params */
    const int    MAX_NEWTON_ITERS = max_newton_iters;
    const double LOOSE_TOL = loose_tol;
    const double TIGHT_TOL = 1e-12;

    memcpy(bends_curr, bends_init, NE*sizeof(double));
    double alpha_curr = 0.0;
    double alpha_step = ALPHA_STEP_INIT;
    int n_steps=0, n_halve=0, newton_total=0;

    for (int attempt = 0; attempt < MAX_ATTEMPTS; attempt++) {
        if (alpha_curr >= ALPHA_FINAL) break;
        double cap = (GOAL - alpha_curr) * 0.5;
        double used = alpha_step < cap ? alpha_step : cap;
        double alpha_try = alpha_curr + used;

        int iters;
        int rc = try_full_newton(alpha_try, bends_curr, bends_trial,
                                  MAX_NEWTON_ITERS, LOOSE_TOL, &iters);
        newton_total += iters;
        if (rc == 0) {
            memcpy(bends_curr, bends_trial, NE*sizeof(double));
            alpha_curr = alpha_try;
            n_steps++;
        } else {
            alpha_step *= 0.5;
            n_halve++;
            if (alpha_step < MIN_ALPHA_STEP) {
                HomotopyResult R = {1, n_steps, n_halve, newton_total, alpha_curr};
                return R;
            }
        }
    }
    /* Reached α_final at loose tol — final tight cleanup at α_curr */
    if (alpha_curr >= ALPHA_FINAL) {
        int iters;
        int rc = try_full_newton(alpha_curr, bends_curr, bends_trial,
                                  30, TIGHT_TOL, &iters);
        newton_total += iters;
        if (rc == 0) {
            memcpy(bends_curr, bends_trial, NE*sizeof(double));
            HomotopyResult R = {0, n_steps, n_halve, newton_total, alpha_curr};
            return R;
        }
    }
    HomotopyResult R = {1, n_steps, n_halve, newton_total, alpha_curr};
    return R;
}

/* Top-level homotopy: try with first-pass params; on failure, retry once
 * with retry params. Sets *retry_used = 1 if pass 2 was the one that
 * succeeded (or failed). The standalone use case wants this to "just
 * work" without external orchestration. */
static HomotopyResult homotopy(const double bends_init[], int *retry_used) {
    HomotopyResult R = homotopy_with(bends_init,
                                       CFG_ALPHA_STEP_INIT_DEG,
                                       CFG_MAX_NEWTON,
                                       CFG_LOOSE_TOL);
    *retry_used = 0;
    if (R.status == 0 || CFG_RETRY_TOL <= 0) return R;
    /* first pass stuck — retry with tighter params */
    *retry_used = 1;
    return homotopy_with(bends_init,
                          CFG_RETRY_INIT_DEG,
                          CFG_RETRY_MAX_NEWTON,
                          CFG_RETRY_TOL);
}

/* ---------- reconstruction (3D coords) -----------------------------------
 *
 * Place base_face vertices at canonical gauge:
 *   base[0] → (0, 0,  1/2)
 *   base[1] → (0, 0, -1/2)
 *   base[2] → (√3/2, 0, 0)
 * BFS dual graph; for each new face with shared edge (a,b) and previous
 * face's third vertex p, place c by
 *   c = m + (√3/2)·(−û·cos θ + v̂·sin θ),
 * where m=(a+b)/2, ê=b−a, û=(p−m)/|p−m|, v̂=û×ê, θ=bend on (a,b). */
static double V_OUT[MAXV+1][3];   /* V_OUT[v] = position of vertex v */
static int    V_PLACED[MAXV+1];

static void vsub(const double a[3], const double b[3], double r[3]) {
    r[0]=a[0]-b[0]; r[1]=a[1]-b[1]; r[2]=a[2]-b[2];
}
static void vadd(const double a[3], const double b[3], double r[3]) {
    r[0]=a[0]+b[0]; r[1]=a[1]+b[1]; r[2]=a[2]+b[2];
}
static void vscale(double s, const double a[3], double r[3]) {
    r[0]=s*a[0]; r[1]=s*a[1]; r[2]=s*a[2];
}
static double vnorm(const double a[3]) {
    return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}
static void vcross(const double a[3], const double b[3], double r[3]) {
    double x = a[1]*b[2] - a[2]*b[1];
    double y = a[2]*b[0] - a[0]*b[2];
    double z = a[0]*b[1] - a[1]*b[0];
    r[0]=x; r[1]=y; r[2]=z;
}

/* place vertex c given placed a, b and previous-face third vertex p, with bend θ */
static void place_third(int a, int b, int p, double theta, int c_out) {
    double m[3], e_hat[3], p_perp[3], u_hat[3], v_hat[3], c_perp[3];
    double tmp[3], tmp2[3];
    /* m = (a+b)/2 */
    vadd(V_OUT[a], V_OUT[b], m); vscale(0.5, m, m);
    /* ê = b - a */
    vsub(V_OUT[b], V_OUT[a], e_hat);
    /* û = (p - m) / |.| */
    vsub(V_OUT[p], m, p_perp);
    double pn = vnorm(p_perp);
    vscale(1.0/pn, p_perp, u_hat);
    /* v̂ = û × ê (then normalize for safety) */
    vcross(u_hat, e_hat, v_hat);
    double vn = vnorm(v_hat);
    vscale(1.0/vn, v_hat, v_hat);
    /* c_perp = -û·cos θ + v̂·sin θ */
    double ct = cos(theta), st = sin(theta);
    vscale(-ct, u_hat, tmp);
    vscale( st, v_hat, tmp2);
    vadd(tmp, tmp2, c_perp);
    /* c = m + (√3/2) · c_perp */
    vscale(sqrt(3.0)/2.0, c_perp, tmp);
    vadd(m, tmp, V_OUT[c_out]);
}

/* BFS reconstruct. Returns 0 on success, -1 if dual graph not connected. */
static int reconstruct(const double bend[]) {
    memset(V_PLACED, 0, (NV+2)*sizeof(int));
    /* base face placement */
    int b0=BASE_VS[0], b1=BASE_VS[1], b2=BASE_VS[2];
    V_OUT[b0][0]=0; V_OUT[b0][1]=0; V_OUT[b0][2]= 0.5; V_PLACED[b0]=1;
    V_OUT[b1][0]=0; V_OUT[b1][1]=0; V_OUT[b1][2]=-0.5; V_PLACED[b1]=1;
    V_OUT[b2][0]=sqrt(3.0)/2.0; V_OUT[b2][1]=0; V_OUT[b2][2]=0; V_PLACED[b2]=1;
    /* find base face index */
    int base_idx = -1;
    for (int i = 0; i < NF; i++) {
        if ((FACES[i].a==b0 || FACES[i].a==b1 || FACES[i].a==b2) &&
            (FACES[i].b==b0 || FACES[i].b==b1 || FACES[i].b==b2) &&
            (FACES[i].c==b0 || FACES[i].c==b1 || FACES[i].c==b2)) {
            base_idx = i; break;
        }
    }
    if (base_idx < 0) return -1;

    static int q[MAXF], placed[MAXF];
    memset(placed, 0, NF*sizeof(int));
    placed[base_idx] = 1;
    int qh = 0, qt = 0;
    q[qt++] = base_idx;
    while (qh < qt) {
        int fi = q[qh++];
        int fa=FACES[fi].a, fb=FACES[fi].b, fc=FACES[fi].c;
        int verts[3] = {fa, fb, fc};
        for (int i = 0; i < 3; i++) {
            int a = verts[i], b = verts[(i+1)%3];
            int other_fi = EM_F[b][a];   /* face with directed edge (b,a) */
            if (other_fi == fi || placed[other_fi]) continue;
            int oa=FACES[other_fi].a, ob=FACES[other_fi].b, oc=FACES[other_fi].c;
            int c = (oa!=a && oa!=b) ? oa : ((ob!=a && ob!=b) ? ob : oc);
            int p = (fa!=a && fa!=b) ? fa : ((fb!=a && fb!=b) ? fb : fc);
            int e = EDGE_IDX[a][b];
            place_third(a, b, p, bend[e], c);
            V_PLACED[c] = 1;
            placed[other_fi] = 1;
            q[qt++] = other_fi;
        }
    }
    for (int v = 1; v <= NV; v++) if (!V_PLACED[v]) return -1;
    return 0;
}

static void write_obj(FILE *fh) {
    for (int v = 1; v <= NV; v++) {
        fprintf(fh, "v %.10f %.10f %.10f\n",
                V_OUT[v][0], V_OUT[v][1], V_OUT[v][2]);
    }
    for (int i = 0; i < NF; i++) {
        fprintf(fh, "f %d %d %d\n", FACES[i].a, FACES[i].b, FACES[i].c);
    }
}

/* ---------- main ---------------------------------------------------------- */
int main(int argc, char **argv) {
    static char line[MAXLINE];
    static double u[MAXV];
    static double bends_init[MAXE];

    /* CLI flags. Two-pass: if pass 1 fails, retry with --retry-* params.
     *   --tol N              pass-1 tolerance (default 1e-3)
     *   --max-newton N       pass-1 Newton max iters (default 8)
     *   --init-step N        pass-1 initial α step in degrees (default 1.0)
     *   --retry-tol N        retry tolerance (default 1e-8)
     *   --retry-max-newton N retry Newton max iters (default 50)
     *   --retry-init-step N  retry initial α step in degrees (default 0.25)
     *   --no-retry           skip retry pass
     *   --obj-dir DIR        write OBJ per successful case to DIR/<n>.obj
     */
    const char *obj_dir = NULL;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--no-retry")) {
            CFG_RETRY_TOL = -1.0;   /* sentinel: skip retry */
            continue;
        }
        if (i >= argc - 1) break;
        if      (!strcmp(argv[i], "--obj-dir"))           obj_dir = argv[i+1];
        else if (!strcmp(argv[i], "--tol"))               CFG_LOOSE_TOL = atof(argv[i+1]);
        else if (!strcmp(argv[i], "--max-newton"))        CFG_MAX_NEWTON = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "--init-step"))         CFG_ALPHA_STEP_INIT_DEG = atof(argv[i+1]);
        else if (!strcmp(argv[i], "--retry-tol"))         CFG_RETRY_TOL = atof(argv[i+1]);
        else if (!strcmp(argv[i], "--retry-max-newton"))  CFG_RETRY_MAX_NEWTON = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "--retry-init-step"))   CFG_RETRY_INIT_DEG = atof(argv[i+1]);
    }
    fprintf(stderr, "puffup: pass1 tol=%g max_newton=%d init=%g°  ",
            CFG_LOOSE_TOL, CFG_MAX_NEWTON, CFG_ALPHA_STEP_INIT_DEG);
    if (CFG_RETRY_TOL > 0)
        fprintf(stderr, "pass2 tol=%g max_newton=%d init=%g°\n",
                CFG_RETRY_TOL, CFG_RETRY_MAX_NEWTON, CFG_RETRY_INIT_DEG);
    else
        fprintf(stderr, "(no retry)\n");

    long n_in = 0, n_ok = 0;
    while (fgets(line, sizeof(line), stdin)) {
        int ll = strlen(line);
        while (ll > 0 && (line[ll-1]=='\n' || line[ll-1]=='\r')) line[--ll] = '\0';
        if (!ll) continue;
        n_in++;

        int pr = parse_facelist(line);
        if (pr < 0) {
            fprintf(stderr, "ERROR: input exceeds MAXF=%d or MAXV=%d "
                            "(rebuild puffup with larger limits)\n", MAXF, MAXV);
            printf("size_too_big %d 0 0 0 0.0 0\n", NV);
            continue;
        }
        if (!pr) {
            printf("parse 0 0 0 0 0.0 0\n");
            continue;
        }
        if (NV > MAXV) {
            fprintf(stderr, "ERROR: NV=%d exceeds MAXV=%d "
                            "(rebuild puffup with larger MAXV)\n", NV, MAXV);
            printf("size_too_big %d 0 0 0 0.0 0\n", NV);
            continue;
        }
        if (build() < 0) {
            printf("deg_too_high %d 0 0 0 0.0 0\n", NV);
            build_clear();
            continue;
        }
        if (horou(u) < 0) {
            printf("horou %d 0 0 0 0.0 0\n", NV);
            build_clear();
            continue;
        }
        compute_bends_at_zero(u, bends_init);
        choose_base_face();

        int retry_used = 0;
        HomotopyResult r = homotopy(bends_init, &retry_used);
        if (r.status == 0) {
            /* solve base bends, then optionally reconstruct + write OBJ */
            if (complete_base_bends(r.final_alpha, bends_curr) < 0) {
                printf("base_bend_fail %d %d %d %d %.10f %d\n",
                       NV, r.n_steps, r.n_halve, r.newton_iters,
                       r.final_alpha * 180.0 / M_PI, retry_used);
                build_clear();
                continue;
            }
            if (obj_dir) {
                if (reconstruct(bends_curr) == 0) {
                    char path[1024];
                    snprintf(path, sizeof(path), "%s/%ld.obj", obj_dir, n_in);
                    FILE *fh = fopen(path, "w");
                    if (fh) { write_obj(fh); fclose(fh); }
                }
            }
        }
        const char *status = (r.status == 0) ? "ok"
                            : (r.status == 1) ? "stuck"
                            : "fail";
        printf("%s %d %d %d %d %.10f %d\n",
               status, NV, r.n_steps, r.n_halve, r.newton_iters,
               r.final_alpha * 180.0 / M_PI, retry_used);
        if (r.status == 0) n_ok++;

        build_clear();
    }
    fflush(stdout);
    fprintf(stderr, "puffup: %ld in, %ld ok\n", n_in, n_ok);
    return 0;
}
