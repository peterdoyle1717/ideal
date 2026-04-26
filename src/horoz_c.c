/*
 * horoz_c.c — ideal horoball weights and vertex positions for triangulated spheres.
 *
 * Input:  face lists from stdin, one per line: "a,b,c;d,e,f;..."
 *         (use clers_decode.py to convert CLERS strings)
 * Output: NV×3 float64 per net to stdout (binary):
 *           record[0..2]         = NaN, NaN, NaN  (vertex 1 = point at infinity)
 *           record[3*(v-1)+0]    = u[v]            (horoball weight)
 *           record[3*(v-1)+1..2] = x[v], y[v]     (flat UHS position)  v=2..NV
 * Record size = NV * 3 * 4 bytes.
 *
 * Compile:  cc -O3 -o horoz_c horoz_c.c -lm
 * Run:      python3 ../clers/python/clers_decode.py < prime/60.txt | ./horoz_c > horoz/60.bin
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXV    400
#define MAXF    (2*MAXV + 4)
#define MAXRING 10
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
#define MAXN MAXV

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

static void horou(double defect, double u_out[]) {
    static int  bndry[MAXV+1];
    static int  int_idx[MAXV+1];
    static int  interior[MAXN];
    static int  ring[MAXV+1][MAXRING];
    static int  ringlen[MAXV+1];
    static int  ff_a[MAXF], ff_b[MAXF], ff_c[MAXF];
    static double xvec[MAXN];

    const double TAU = 2.0 * M_PI;
    const double target = TAU - defect;

    memset(bndry, 0, (NV+2)*sizeof(int));
    memset(int_idx, -1, (NV+2)*sizeof(int));
    int bndry_ring[MAXN], n_bndry = cyclic_nbrs(1, bndry_ring);
    for (int i = 0; i < n_bndry; i++) bndry[bndry_ring[i]] = 1;

    int n_int = 0;
    for (int v = 2; v <= NV; v++) {
        if (!bndry[v]) { int_idx[v] = n_int; interior[n_int++] = v; }
    }

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

    if (n_int == 0) goto done;

    for (int i = 0; i < n_int; i++) xvec[i] = 1.0;

    for (int iter = 0; iter < 200; iter++) {
        double res = 0.0;
        for (int i = 0; i < n_int; i++) {
            int v = interior[i]; int k = ringlen[v]; double ui = xvec[i];
            double s = 0.0;
            for (int j = 0; j < k; j++)
                s += petal(ui, U(ring[v][j]), U(ring[v][(j+1)%k]));
            Fvec[i] = s - target;
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
        if (lu_solve(Jmat, dxvec, n_int) < 0) break;

        double step = 1.0;
        for (int bt = 0; bt < 60; bt++, step *= 0.5) {
            int ok = 1;
            for (int i = 0; i < n_int; i++) {
                if (xvec[i] + step*dxvec[i] <= 0.0) { ok=0; break; }
            }
            if (!ok) continue;
            for (int f = 0; f < nff && ok; f++) {
                double ua = U(ff_a[f]), ub = U(ff_b[f]), uc = U(ff_c[f]);
                if (int_idx[ff_a[f]] >= 0) ua = xvec[int_idx[ff_a[f]]] + step*dxvec[int_idx[ff_a[f]]];
                if (int_idx[ff_b[f]] >= 0) ub = xvec[int_idx[ff_b[f]]] + step*dxvec[int_idx[ff_b[f]]];
                if (int_idx[ff_c[f]] >= 0) uc = xvec[int_idx[ff_c[f]]] + step*dxvec[int_idx[ff_c[f]]];
                double p=ua*ub, q=ua*uc, r=ub*uc;
                if (p+q<=r || p+r<=q || q+r<=p) { ok=0; }
            }
            if (!ok) continue;
            double res2 = 0.0;
            for (int i = 0; i < n_int; i++) {
                int v=interior[i]; int k=ringlen[v];
                double ui = xvec[i] + step*dxvec[i];
                double s = 0.0;
                for (int j = 0; j < k; j++) {
                    double uj=U(ring[v][j]), uk2=U(ring[v][(j+1)%k]);
                    if (int_idx[ring[v][j]]       >= 0) uj  = xvec[int_idx[ring[v][j]]]       + step*dxvec[int_idx[ring[v][j]]];
                    if (int_idx[ring[v][(j+1)%k]] >= 0) uk2 = xvec[int_idx[ring[v][(j+1)%k]]] + step*dxvec[int_idx[ring[v][(j+1)%k]]];
                    s += petal(ui, uj, uk2);
                }
                double af = fabs(s - target); if (af > res2) res2 = af;
            }
            if (res2 < res) break;
        }

        for (int i = 0; i < n_int; i++) xvec[i] += step * dxvec[i];
    }

done:
    u_out[0] = NAN;
    for (int v = 2; v <= NV; v++) {
        u_out[v-1] = bndry[v] ? 1.0 : xvec[int_idx[v]];
    }
    #undef U
}

/* ── thirdpoint ──────────────────────────────────────────────────────────── */
/*
 * Given directed edge a→b with face (a,b,c) CCW, and horoball edge lengths
 * dA=u[a]*u[c], dB=u[b]*u[c], find (xc,yc) in the upper half-plane.
 */
static void thirdpoint(double xa, double ya,
                       double xb, double yb,
                       double dA, double dB,
                       double *xc, double *yc) {
    double dx = xb - xa, dy = yb - ya;
    double L = sqrt(dx*dx + dy*dy);
    double d0 = dA / L, d1 = dB / L;
    double xn = (1.0 + d0*d0 - d1*d1) / 2.0;
    double yn2 = d0*d0 - xn*xn;
    double yn = sqrt(yn2 > 0.0 ? yn2 : 0.0);
    *xc = xa + dx*xn + dy*yn;
    *yc = ya + dy*xn - dx*yn;
}

/* ── horoz BFS ───────────────────────────────────────────────────────────── */
/*
 * Compute flat vertex positions (x,y) in the upper half-plane from u[].
 * Fixes v=2 at (0,0) and v=3 at (1,0); BFS to place remaining vertices.
 */
static void horoz(double u[], double out[]) {
    static double x[MAXV+1], y[MAXV+1];
    static int placed[MAXV+1];
    static int qa[MAXV*2+8], qb[MAXV*2+8];

    memset(placed, 0, (NV+2)*sizeof(int));

    placed[1] = 1;
    placed[2] = 1; x[2] = 0.0; y[2] = 0.0;
    placed[3] = 1; x[3] = 1.0; y[3] = 0.0;

    int qh = 0, qt = 0;
    qa[qt] = 3; qb[qt] = 2; qt++;

    while (qh < qt) {
        int a = qa[qh], b = qb[qh]; qh++;
        int c = EM[a][b];
        if (c == 0 || c == 1 || placed[c]) continue;

        double dA = u[a-1] * u[c-1];
        double dB = u[b-1] * u[c-1];
        thirdpoint(x[a], y[a], x[b], y[b], dA, dB, &x[c], &y[c]);
        placed[c] = 1;

        qa[qt] = c; qb[qt] = b; qt++;
        qa[qt] = a; qb[qt] = c; qt++;
    }

    /* cleanup: place any vertex not reached by BFS */
    {
        int changed = 1;
        while (changed) {
            changed = 0;
            for (int i = 0; i < NF; i++) {
                int a=F[i].a, b=F[i].b, c=F[i].c;
                if (a==1 || b==1 || c==1) continue;
                if (placed[a] && placed[b] && !placed[c]) {
                    thirdpoint(x[a],y[a],x[b],y[b],u[a-1]*u[c-1],u[b-1]*u[c-1],&x[c],&y[c]);
                    placed[c]=1; changed=1;
                } else if (placed[b] && placed[c] && !placed[a]) {
                    thirdpoint(x[b],y[b],x[c],y[c],u[b-1]*u[a-1],u[c-1]*u[a-1],&x[a],&y[a]);
                    placed[a]=1; changed=1;
                } else if (placed[a] && placed[c] && !placed[b]) {
                    thirdpoint(x[c],y[c],x[a],y[a],u[c-1]*u[b-1],u[a-1]*u[b-1],&x[b],&y[b]);
                    placed[b]=1; changed=1;
                }
            }
        }
    }

    out[0] = NAN; out[1] = NAN; out[2] = NAN;
    for (int v = 2; v <= NV; v++) {
        out[3*(v-1)+0] = u[v-1];
        out[3*(v-1)+1] = x[v];
        out[3*(v-1)+2] = y[v];
    }
}

/* ── main ────────────────────────────────────────────────────────────────── */
int main(void) {
    static char  line[MAXLINE];
    static double u[MAXV];
    static double z_out[MAXV*3];

    while (fgets(line, sizeof(line), stdin)) {
        int ll = strlen(line);
        while (ll > 0 && (line[ll-1]=='\n' || line[ll-1]=='\r')) line[--ll] = '\0';
        if (!ll) continue;

        if (!parse_facelist(line)) { fprintf(stderr, "parse failed: %.60s\n", line); continue; }
        build();
        horou(0.0, u);
        horoz(u, z_out);
        build_clear();

        fwrite(z_out, sizeof(double), NV*3, stdout);
    }
    return 0;
}
