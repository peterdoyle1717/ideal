/*
 * proof_c.c — sub/supersolution proof checker for ideal horoball packings.
 *
 * For each net, computes:
 *   umin = horou(-eps)   — subsolution (angle sum = 2π+eps > 2π at every interior vertex)
 *   umax = horou(+eps)   — supersolution (angle sum = 2π-eps < 2π at every interior vertex)
 * then runs 5 checks that together prove a true solution exists between them.
 *
 * Input:  face lists from stdin, one per line: "a,b,c;d,e,f;..."
 *         (use clers_decode.py to convert CLERS strings)
 * Output: one float32 per net to stdout (binary):
 *           value > 0  → proved (slack/eps; order-1 means healthy margin)
 *           value ≤ 0  → some check failed
 *           NaN        → solver failed
 *
 * Compile:  cc -O3 -o proof_c proof_c.c -lm
 * Run:      python3 ../clers/python/clers_decode.py < prime/60.txt | ./proof_c > proof_60.bin
 *
 * Optional arg: fixed eps (e.g. "0.000125" = 1/8000).
 * In fixed mode: output 1.0 (proved) or NaN (not proved).
 * Used by find_eps.sh. Default: adaptive halving from 1/500.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXV    400
#define MAXF    (2*MAXV + 4)
#define MAXRING 10
#define MAXN    MAXV
#define MAXLINE 16384

/* ── graph state ─────────────────────────────────────────────────────────── */
typedef struct { int a, b, c; } Face;
static int  NV, NF;
static Face F[MAXF];
static int  EM[MAXV+1][MAXV+1];
static int  DEG[MAXV+1];
static int  NBR[MAXV+1][MAXRING];
static int  NNBR[MAXV+1];
static short DU[MAXF*3], DW[MAXF*3]; static int ND;

static void nbr_add(int u, int w) {
    for (int i = 0; i < NNBR[u]; i++) if (NBR[u][i] == w) return;
    NBR[u][NNBR[u]++] = w;
}
static void build_clear(void) {
    for (int i = 0; i < ND; i++) EM[DU[i]][DW[i]] = 0;
    ND = 0;
}
static void build(void) {
    memset(DEG,  0, (NV+2)*sizeof(int));
    memset(NNBR, 0, (NV+2)*sizeof(int));
    ND = 0;
    for (int i = 0; i < NF; i++) {
        int a=F[i].a, b=F[i].b, c=F[i].c;
        EM[a][b]=c; DU[ND]=a; DW[ND]=b; ND++;
        EM[b][c]=a; DU[ND]=b; DW[ND]=c; ND++;
        EM[c][a]=b; DU[ND]=c; DW[ND]=a; ND++;
        DEG[a]++; DEG[b]++; DEG[c]++;
        nbr_add(a,b); nbr_add(b,a);
        nbr_add(b,c); nbr_add(c,b);
        nbr_add(a,c); nbr_add(c,a);
    }
}

/* ── face list parser ────────────────────────────────────────────────────── */
static int parse_facelist(const char *line) {
    NF = 0; NV = 0;
    const char *p = line;
    while (*p) {
        int a, b, c;
        if (sscanf(p, "%d,%d,%d", &a, &b, &c) != 3) break;
        F[NF++] = (Face){a, b, c};
        if (a > NV) NV = a;
        if (b > NV) NV = b;
        if (c > NV) NV = c;
        while (*p && *p != ';') p++;
        if (*p == ';') p++;
    }
    return NF;
}

/* ── cyclic neighbors ────────────────────────────────────────────────────── */
static int cyclic_nbrs(int v, int ring[]) {
    if (!NNBR[v]) return 0;
    int start = NBR[v][0];
    ring[0] = start;
    int k = 1, cur = start;
    for (;;) {
        int nxt = EM[v][cur];
        if (nxt == start) break;
        ring[k++] = nxt;
        cur = nxt;
    }
    return k;
}

/* ── horoball geometry ───────────────────────────────────────────────────── */
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

/* ── Newton solver ───────────────────────────────────────────────────────── */
static double Jmat[MAXN][MAXN];
static double Fvec[MAXN];
static double dxvec[MAXN];

static int lu_solve(double J[][MAXN], double b[], int n) {
    for (int col = 0; col < n; col++) {
        int pivot = col;
        double best = fabs(J[col][col]);
        for (int row = col+1; row < n; row++) {
            if (fabs(J[row][col]) > best) { best = fabs(J[row][col]); pivot = row; }
        }
        if (best < 1e-14) return -1;
        if (pivot != col) {
            for (int k = col; k < n; k++) { double t=J[col][k]; J[col][k]=J[pivot][k]; J[pivot][k]=t; }
            { double t=b[col]; b[col]=b[pivot]; b[pivot]=t; }
        }
        double inv = 1.0 / J[col][col];
        for (int row = col+1; row < n; row++) {
            double fac = J[row][col] * inv;
            for (int k = col; k < n; k++) J[row][k] -= fac * J[col][k];
            b[row] -= fac * b[col];
        }
    }
    for (int i = n-1; i >= 0; i--) {
        double s = b[i];
        for (int j = i+1; j < n; j++) s -= Jmat[i][j] * b[j];
        b[i] = s / J[i][i];
    }
    return 0;
}

static int horou(double defect, double u_out[]) {
    static int    bndry[MAXV+1];
    static int    int_idx[MAXV+1];
    static int    interior[MAXN];
    static int    ring[MAXV+1][MAXRING];
    static int    ringlen[MAXV+1];
    static int    ff_a[MAXF], ff_b[MAXF], ff_c[MAXF];
    static double xvec[MAXN];

    const double TAU = 2.0 * M_PI;
    const double target = TAU - defect;

    memset(bndry,   0, (NV+2)*sizeof(int));
    memset(int_idx,-1, (NV+2)*sizeof(int));
    int bndry_ring[MAXN], n_bndry = cyclic_nbrs(1, bndry_ring);
    for (int i = 0; i < n_bndry; i++) bndry[bndry_ring[i]] = 1;

    int n_int = 0;
    for (int v = 2; v <= NV; v++)
        if (!bndry[v]) { int_idx[v] = n_int; interior[n_int++] = v; }

    #define U(v) (bndry[v] ? 1.0 : xvec[int_idx[v]])

    for (int i = 0; i < n_int; i++) {
        int v = interior[i];
        ringlen[v] = cyclic_nbrs(v, ring[v]);
    }

    int nff = 0;
    for (int i = 0; i < NF; i++) {
        int a=F[i].a, b=F[i].b, c=F[i].c;
        if (a!=1 && b!=1 && c!=1) { ff_a[nff]=a; ff_b[nff]=b; ff_c[nff]=c; nff++; }
    }

    u_out[0] = NAN;
    for (int v = 2; v <= NV; v++)
        u_out[v-1] = bndry[v] ? 1.0 : 0.0;

    if (n_int == 0) goto done;

    for (int i = 0; i < n_int; i++) xvec[i] = 1.0;

    for (int iter = 0; iter < 200; iter++) {
        double res = 0.0;
        for (int i = 0; i < n_int; i++) {
            int v = interior[i]; int k = ringlen[v]; double ui = xvec[i];
            double sv = 0.0;
            for (int j = 0; j < k; j++)
                sv += petal(ui, U(ring[v][j]), U(ring[v][(j+1)%k]));
            Fvec[i] = sv - target;
            double af = fabs(Fvec[i]); if (af > res) res = af;
        }
        if (res < 1e-10) break;

        memset(Jmat, 0, sizeof(double)*n_int*MAXN);
        for (int i = 0; i < n_int; i++) {
            int v = interior[i]; int k = ringlen[v]; double ui = xvec[i];
            for (int j = 0; j < k; j++) {
                int vj = ring[v][j], vk = ring[v][(j+1)%k];
                double uj = U(vj), uk = U(vk);
                double dui, duj, duk;
                petal_grad(ui, uj, uk, &dui, &duj, &duk);
                Jmat[i][i] += dui;
                if (int_idx[vj] >= 0) Jmat[i][int_idx[vj]] += duj;
                if (int_idx[vk] >= 0) Jmat[i][int_idx[vk]] += duk;
            }
        }

        for (int i = 0; i < n_int; i++) dxvec[i] = -Fvec[i];
        if (lu_solve(Jmat, dxvec, n_int) < 0) goto fail;

        double step = 1.0;
        int accepted = 0;
        for (int bt = 0; bt < 60; bt++, step *= 0.5) {
            int ok = 1;
            for (int i = 0; i < n_int; i++)
                if (xvec[i] + step*dxvec[i] <= 0.0) { ok=0; break; }
            if (!ok) continue;
            for (int fi = 0; fi < nff && ok; fi++) {
                double ua = U(ff_a[fi]), ub = U(ff_b[fi]), uc = U(ff_c[fi]);
                if (int_idx[ff_a[fi]] >= 0) ua = xvec[int_idx[ff_a[fi]]] + step*dxvec[int_idx[ff_a[fi]]];
                if (int_idx[ff_b[fi]] >= 0) ub = xvec[int_idx[ff_b[fi]]] + step*dxvec[int_idx[ff_b[fi]]];
                if (int_idx[ff_c[fi]] >= 0) uc = xvec[int_idx[ff_c[fi]]] + step*dxvec[int_idx[ff_c[fi]]];
                double p=ua*ub, q=ua*uc, r=ub*uc;
                if (p+q<=r || p+r<=q || q+r<=p) { ok=0; }
            }
            if (!ok) continue;
            double res2 = 0.0;
            for (int i = 0; i < n_int; i++) {
                int v=interior[i]; int k=ringlen[v];
                double ui = xvec[i] + step*dxvec[i];
                double sv = 0.0;
                for (int j = 0; j < k; j++) {
                    double uj=U(ring[v][j]), uk2=U(ring[v][(j+1)%k]);
                    if (int_idx[ring[v][j]]       >= 0) uj  = xvec[int_idx[ring[v][j]]]       + step*dxvec[int_idx[ring[v][j]]];
                    if (int_idx[ring[v][(j+1)%k]] >= 0) uk2 = xvec[int_idx[ring[v][(j+1)%k]]] + step*dxvec[int_idx[ring[v][(j+1)%k]]];
                    sv += petal(ui, uj, uk2);
                }
                double af = fabs(sv - target); if (af > res2) res2 = af;
            }
            if (res2 < res) { accepted=1; break; }
        }
        if (!accepted) goto fail;

        for (int i = 0; i < n_int; i++) xvec[i] += step * dxvec[i];
    }

done:
    for (int v = 2; v <= NV; v++)
        u_out[v-1] = bndry[v] ? 1.0 : xvec[int_idx[v]];
    #undef U
    return 1;

fail:
    #undef U
    return 0;
}

/* ── proof checks ─────────────────────────────────────────────────────────── */
/*
 * Five checks that together prove a solution u* exists with umin < u* < umax.
 * Returns minimum slack across all checks; positive = proved.
 */
static double proof(double umin[], double umax[],
                    int bndry[], int interior[], int n_int,
                    int ring[MAXV+1][MAXRING], int ringlen[]) {
    const double TAU = 2.0 * M_PI;
    double slack = 1e30;

    /* 1. monocheck: umax[v] > umin[v] for all interior vertices */
    for (int i = 0; i < n_int; i++) {
        int v = interior[i];
        double s = umax[v-1] - umin[v-1];
        if (s < slack) slack = s;
    }
    if (slack <= 0.0) return slack;

    /* 2. excesscheck: umin is sub-solution, umax is super-solution */
    for (int i = 0; i < n_int; i++) {
        int v = interior[i]; int k = ringlen[v];
        double fn_min = 0.0, fn_max = 0.0;
        for (int j = 0; j < k; j++) {
            int vj = ring[v][j], vk = ring[v][(j+1)%k];
            fn_min += petal(umin[v-1], umin[vj-1], umin[vk-1]);
            fn_max += petal(umax[v-1], umax[vj-1], umax[vk-1]);
        }
        double s1 = fn_min - TAU;   /* > 0: angle sum exceeds 2π — subsolution */
        double s2 = TAU - fn_max;   /* > 0: angle sum below 2π — supersolution */
        if (s1 < slack) slack = s1;
        if (s2 < slack) slack = s2;
    }
    if (slack <= 0.0) return slack;

    /* 3. tricheck: triangle inequality holds for any u in [umin, umax] */
    for (int i = 0; i < NF; i++) {
        int a=F[i].a, b=F[i].b, c=F[i].c;
        if (a==1 || b==1 || c==1) continue;
        double uma=umin[a-1], umb=umin[b-1], umc=umin[c-1];
        double xma=umax[a-1], xmb=umax[b-1], xmc=umax[c-1];
        double s1 = uma*(xmb+xmc)/(xmb*xmc) - 1.0;
        double s2 = umb*(xma+xmc)/(xma*xmc) - 1.0;
        double s3 = umc*(xma+xmb)/(xma*xmb) - 1.0;
        if (s1 < slack) slack = s1;
        if (s2 < slack) slack = s2;
        if (s3 < slack) slack = s3;
    }
    if (slack <= 0.0) return slack;

    /* 4. convexcheck: Delaunay / butterfly condition for each interior edge */
    for (int v = 2; v <= NV; v++) {
        for (int w = v+1; w <= NV; w++) {
            if (!EM[v][w]) continue;
            int b = EM[v][w], d = EM[w][v];
            if (b==1 || d==1) continue;
            double ua2=umin[v-1]*umin[v-1], uc2=umin[w-1]*umin[w-1];
            double ub2=umin[b-1]*umin[b-1], ud2=umin[d-1]*umin[d-1];
            double Va2=umax[v-1]*umax[v-1], Vc2=umax[w-1]*umax[w-1];
            double Vb2=umax[b-1]*umax[b-1], Vd2=umax[d-1]*umax[d-1];
            double num = 2.0 * ub2 * ud2 * (ua2 + uc2);
            double den = Va2 * Vc2 * (Vb2 + Vd2);
            double s = num / den - 1.0;
            if (s < slack) slack = s;
        }
    }
    if (slack <= 0.0) return slack;

    /* 5. boundarycheck: boundary angles < π */
    {
        int bndry_ring[MAXV+1];
        int n_bndry = cyclic_nbrs(1, bndry_ring);
        for (int bi = 0; bi < n_bndry; bi++) {
            int b = bndry_ring[bi];
            int b_ring[MAXRING];
            int b_len = cyclic_nbrs(b, b_ring);
            int start = -1;
            for (int j = 0; j < b_len; j++) {
                if (b_ring[j] == 1) { start = (j+1) % b_len; break; }
            }
            if (start < 0) continue;
            int arc[MAXRING]; int arc_len = 0;
            for (int j = 0; j < b_len - 1; j++)
                arc[arc_len++] = b_ring[(start + j) % b_len];
            double slide = 0.0;
            for (int j = 0; j < arc_len - 1; j++) {
                int vj = arc[j], vk = arc[j+1];
                slide += petal(umin[b-1], umax[vj-1], umax[vk-1]);
            }
            double s = (M_PI - slide) / (2.0 * M_PI);
            if (s < slack) slack = s;
        }
    }

    return slack;
}

/* ── main ────────────────────────────────────────────────────────────────── */
#define EPS_START    0.002      /* 1/500 */
#define EPS_HALVINGS 20

int main(int argc, char *argv[]) {
    static char   line[MAXLINE];
    static double umin_d[MAXV], umax_d[MAXV];
    static int    bndry[MAXV+1], int_idx[MAXV+1], interior_v[MAXN];
    static int    ring[MAXV+1][MAXRING];
    static int    ringlen[MAXV+1];

    double fixed_eps = (argc > 1) ? atof(argv[1]) : -1.0;

    while (fgets(line, sizeof(line), stdin)) {
        int ll = strlen(line);
        while (ll > 0 && (line[ll-1]=='\n' || line[ll-1]=='\r')) line[--ll] = '\0';
        if (!ll) continue;

        if (!parse_facelist(line)) { fprintf(stderr, "parse failed: %.60s\n", line); continue; }
        build();

        memset(bndry,   0, (NV+2)*sizeof(int));
        memset(int_idx,-1, (NV+2)*sizeof(int));
        int bndry_ring[MAXV+1], n_bndry = cyclic_nbrs(1, bndry_ring);
        for (int i = 0; i < n_bndry; i++) bndry[bndry_ring[i]] = 1;
        int n_int = 0;
        for (int v = 2; v <= NV; v++)
            if (!bndry[v]) { int_idx[v] = n_int; interior_v[n_int++] = v; }
        for (int i = 0; i < n_int; i++)
            ringlen[interior_v[i]] = cyclic_nbrs(interior_v[i], ring[interior_v[i]]);

        float result = NAN;
        if (fixed_eps > 0.0) {
            if (horou(-fixed_eps, umin_d) && horou(+fixed_eps, umax_d)) {
                double s = proof(umin_d, umax_d, bndry, interior_v, n_int, ring, ringlen);
                if (s > 0.0) result = 1.0f;
            }
        } else {
            double eps = EPS_START;
            for (int attempt = 0; attempt <= EPS_HALVINGS; attempt++, eps *= 0.5) {
                if (!horou(-eps, umin_d) || !horou(+eps, umax_d)) continue;
                double s = proof(umin_d, umax_d, bndry, interior_v, n_int, ring, ringlen);
                if (s > 0.0) { result = (float)(s / eps); break; }
            }
        }

        build_clear();
        fwrite(&result, sizeof(float), 1, stdout);
    }
    return 0;
}
