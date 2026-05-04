/*
  euclid_prover.c

  First dependency-free C prototype of an Ellison-style verifier.

  Build:
      cc -O3 -std=c11 -Wall -Wextra -pedantic euclid_prover.c -lm -o euclid_prover

  Do not compile with -ffast-math. The verification layer assumes ordinary
  IEEE-754 binary64 operations and nextafter padding.

  This file intentionally keeps both the untrusted candidate generation and the
  trusted witness checks visible in one translation unit.  The untrusted code
  proposes a singular-value lower candidate and a Cholesky witness; the trusted
  code checks the residual bounds using IEEE-backed interval arithmetic.
*/

#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct { double x, y, z; } Vec3;
typedef struct { int a, b, c; } Face;
typedef struct { int a, b; } Edge;
typedef struct { double lo, hi; } Interval;
typedef struct { int n; int v[3]; } Simplex;

static double ddown(double x) { return nextafter(x, -INFINITY); }
static double dup(double x) { return nextafter(x, INFINITY); }

static void die(const char *msg) {
    fprintf(stderr, "error: %s\n", msg);
    exit(1);
}

static void *xcalloc(size_t n, size_t size) {
    void *p = calloc(n, size);
    if (!p) die("out of memory");
    return p;
}

static void *xrealloc(void *p, size_t size) {
    void *q = realloc(p, size);
    if (!q) die("out of memory");
    return q;
}

static Interval ipoint(double x) {
    if (isnan(x)) die("NaN interval point");
    Interval r = {x, x};
    return r;
}

static Interval izero(void) { return ipoint(0.0); }

static Interval iadd(Interval a, Interval b) {
    Interval r = { ddown(a.lo + b.lo), dup(a.hi + b.hi) };
    if (isnan(r.lo) || isnan(r.hi) || r.lo > r.hi) die("bad interval add");
    return r;
}

static Interval isub(Interval a, Interval b) {
    Interval r = { ddown(a.lo - b.hi), dup(a.hi - b.lo) };
    if (isnan(r.lo) || isnan(r.hi) || r.lo > r.hi) die("bad interval sub");
    return r;
}

static Interval imul(Interval a, Interval b) {
    double v0 = a.lo * b.lo;
    double v1 = a.lo * b.hi;
    double v2 = a.hi * b.lo;
    double v3 = a.hi * b.hi;
    double lo = fmin(fmin(v0, v1), fmin(v2, v3));
    double hi = fmax(fmax(v0, v1), fmax(v2, v3));
    Interval r = { ddown(lo), dup(hi) };
    if (isnan(r.lo) || isnan(r.hi) || r.lo > r.hi) die("bad interval mul");
    return r;
}


static Interval idiv(Interval a, Interval b) {
    if (b.lo <= 0.0 && 0.0 <= b.hi) die("interval division by interval containing zero");
    double v0 = a.lo / b.lo;
    double v1 = a.lo / b.hi;
    double v2 = a.hi / b.lo;
    double v3 = a.hi / b.hi;
    double lo = fmin(fmin(v0, v1), fmin(v2, v3));
    double hi = fmax(fmax(v0, v1), fmax(v2, v3));
    Interval r = { ddown(lo), dup(hi) };
    if (isnan(r.lo) || isnan(r.hi) || r.lo > r.hi) die("bad interval div");
    return r;
}

static Interval isquare(Interval a) {
    if (a.lo <= 0.0 && 0.0 <= a.hi) {
        double hi = fmax(a.lo * a.lo, a.hi * a.hi);
        Interval r = { 0.0, dup(hi) };
        return r;
    }
    double v0 = a.lo * a.lo;
    double v1 = a.hi * a.hi;
    /* Squares are non-negative. Outward rounding ddown(min(v0,v1)) can
       step below 0 if the product underflows; clamp to recover. */
    double lo = ddown(fmin(v0, v1));
    if (lo < 0.0) lo = 0.0;
    Interval r = { lo, dup(fmax(v0, v1)) };
    return r;
}

static Interval clamp_nonneg(Interval a) {
    if (a.lo < 0.0) a.lo = 0.0;
    if (a.hi < 0.0) a.hi = 0.0;
    return a;
}

static Interval isqrt_interval(Interval a) {
    if (a.lo < 0.0) die("sqrt interval with negative lower endpoint");
    Interval r = { ddown(sqrt(a.lo)), dup(sqrt(a.hi)) };
    return r;
}

/* Status-returning sqrt for the sigma-certificate path. Failures inside
   certify_factor_witness must propagate up so certify_sigma_lower can
   shrink s and retry, mirroring the Python ArithmeticError handling. */
static bool isqrt_interval_checked(Interval a, Interval *out, char *msg, size_t msglen) {
    if (a.hi < 0.0) {
        snprintf(msg, msglen, "sqrt of strictly-negative interval: [%.17g, %.17g]", a.lo, a.hi);
        return false;
    }
    if (a.lo < 0.0) {
        snprintf(msg, msglen, "sqrt interval has negative lower endpoint: [%.17g, %.17g]", a.lo, a.hi);
        return false;
    }
    out->lo = ddown(sqrt(a.lo));
    out->hi = dup(sqrt(a.hi));
    return true;
}

static double iabs_upper(Interval a) {
    return dup(fmax(fabs(a.lo), fabs(a.hi)));
}

/* Status-returning Frobenius bound, for sigma-certificate callers. */
static bool frob_upper_intervals_checked(const Interval *x, int n, double *out, char *msg, size_t msglen) {
    Interval total = izero();
    for (int i = 0; i < n; i++) {
        Interval a = ipoint(iabs_upper(x[i]));
        total = iadd(total, isquare(a));
    }
    Interval root;
    if (!isqrt_interval_checked(clamp_nonneg(total), &root, msg, msglen)) return false;
    *out = root.hi;
    return true;
}

static double vecdot(Vec3 a, Vec3 b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

static Vec3 vecsub(Vec3 a, Vec3 b) {
    Vec3 r = {a.x-b.x, a.y-b.y, a.z-b.z};
    return r;
}

static Vec3 vecadd(Vec3 a, Vec3 b) {
    Vec3 r = {a.x+b.x, a.y+b.y, a.z+b.z};
    return r;
}

static Vec3 vecscale(Vec3 a, double s) {
    Vec3 r = {a.x*s, a.y*s, a.z*s};
    return r;
}

static double vecnorm(Vec3 a) { return sqrt(vecdot(a, a)); }

static Vec3 veccross(Vec3 a, Vec3 b) {
    Vec3 r = {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
    return r;
}

static Vec3 vecnormalize(Vec3 a) {
    double n = vecnorm(a);
    if (!(n > 0.0) || !isfinite(n)) die("cannot normalize zero/nonfinite vector");
    return vecscale(a, 1.0 / n);
}


static void parse_obj(const char *path, Vec3 **vout, int *nvout, Face **fout, int *nfout) {
    FILE *fp = fopen(path, "r");
    if (!fp) {
        perror(path);
        exit(1);
    }
    int vc = 0, vcap = 64, fc = 0, fcap = 64;
    Vec3 *verts = (Vec3 *)xcalloc((size_t)vcap, sizeof(Vec3));
    Face *faces = (Face *)xcalloc((size_t)fcap, sizeof(Face));
    char line[4096];
    while (fgets(line, sizeof(line), fp)) {
        char *s = line;
        while (isspace((unsigned char)*s)) s++;
        if (*s == '\0' || *s == '#') continue;
        if (s[0] == 'v' && isspace((unsigned char)s[1])) {
            char *p = s + 1;
            errno = 0;
            double x = strtod(p, &p);
            double y = strtod(p, &p);
            double z = strtod(p, &p);
            if (errno || !isfinite(x) || !isfinite(y) || !isfinite(z)) die("bad OBJ vertex");
            if (vc == vcap) {
                vcap *= 2;
                verts = (Vec3 *)xrealloc(verts, (size_t)vcap * sizeof(Vec3));
            }
            verts[vc++] = (Vec3){x,y,z};
        } else if (s[0] == 'f' && isspace((unsigned char)s[1])) {
            int ids[3];
            char *p = s + 1;
            for (int k = 0; k < 3; k++) {
                while (isspace((unsigned char)*p)) p++;
                if (*p == '\0') die("bad OBJ face");
                errno = 0;
                long idx = strtol(p, &p, 10);
                if (errno || idx == 0) die("bad OBJ face index");
                if (idx < 0) idx = vc + idx + 1;
                ids[k] = (int)idx - 1;
                while (*p && !isspace((unsigned char)*p)) p++;
            }
            while (isspace((unsigned char)*p)) p++;
            if (*p && *p != '#') die("only triangular OBJ faces are supported");
            if (ids[0] == ids[1] || ids[0] == ids[2] || ids[1] == ids[2]) die("degenerate OBJ face");
            if (fc == fcap) {
                fcap *= 2;
                faces = (Face *)xrealloc(faces, (size_t)fcap * sizeof(Face));
            }
            faces[fc++] = (Face){ids[0], ids[1], ids[2]};
        }
    }
    fclose(fp);
    if (vc == 0 || fc == 0) die("OBJ needs vertices and triangular faces");
    for (int i = 0; i < fc; i++) {
        if (faces[i].a < 0 || faces[i].a >= vc || faces[i].b < 0 || faces[i].b >= vc || faces[i].c < 0 || faces[i].c >= vc) die("OBJ face index out of range");
    }
    *vout = verts; *nvout = vc; *fout = faces; *nfout = fc;
}

static int edge_cmp(const void *pa, const void *pb) {
    const Edge *a = (const Edge *)pa, *b = (const Edge *)pb;
    if (a->a != b->a) return a->a - b->a;
    return a->b - b->b;
}

static Edge make_edge(int a, int b) {
    if (a > b) { int t = a; a = b; b = t; }
    return (Edge){a,b};
}

static Edge *edges_from_faces(const Face *faces, int nf, int *neout) {
    int nraw = 3 * nf;
    Edge *raw = (Edge *)xcalloc((size_t)nraw, sizeof(Edge));
    for (int i = 0; i < nf; i++) {
        raw[3*i+0] = make_edge(faces[i].a, faces[i].b);
        raw[3*i+1] = make_edge(faces[i].b, faces[i].c);
        raw[3*i+2] = make_edge(faces[i].c, faces[i].a);
    }
    qsort(raw, (size_t)nraw, sizeof(Edge), edge_cmp);
    int ne = 0;
    for (int i = 0; i < nraw; i++) {
        if (ne == 0 || raw[i].a != raw[ne-1].a || raw[i].b != raw[ne-1].b) raw[ne++] = raw[i];
    }
    *neout = ne;
    return raw;
}

#define MAT(A,n,i,j) ((A)[(size_t)(i)*(size_t)(n)+(size_t)(j)])

static double *build_jacobian(const Vec3 *v, int nv, const Edge *edges, int ne, int *colsout) {
    int cols = 3 * nv;
    double *J = (double *)xcalloc((size_t)ne * (size_t)cols, sizeof(double));
    for (int r = 0; r < ne; r++) {
        int a = edges[r].a, b = edges[r].b;
        Vec3 d = vecscale(vecsub(v[a], v[b]), 2.0);
        J[(size_t)r*cols + 3*a + 0] = d.x;
        J[(size_t)r*cols + 3*a + 1] = d.y;
        J[(size_t)r*cols + 3*a + 2] = d.z;
        J[(size_t)r*cols + 3*b + 0] = -d.x;
        J[(size_t)r*cols + 3*b + 1] = -d.y;
        J[(size_t)r*cols + 3*b + 2] = -d.z;
    }
    *colsout = cols;
    return J;
}

static Interval interval_length_squared(const Vec3 *v, int a, int b) {
    Interval total = izero();
    Interval dx = isub(ipoint(v[a].x), ipoint(v[b].x)); total = iadd(total, isquare(dx));
    Interval dy = isub(ipoint(v[a].y), ipoint(v[b].y)); total = iadd(total, isquare(dy));
    Interval dz = isub(ipoint(v[a].z), ipoint(v[b].z)); total = iadd(total, isquare(dz));
    return total;
}

static double certify_rho_upper(const Vec3 *v, const Edge *edges, int ne) {
    Interval total = izero();
    Interval one = ipoint(1.0);
    for (int i = 0; i < ne; i++) {
        Interval l2 = interval_length_squared(v, edges[i].a, edges[i].b);
        Interval err = isub(l2, one);
        total = iadd(total, isquare(err));
    }
    return isqrt_interval(clamp_nonneg(total)).hi;
}

static double *compute_jjt(const double *J, int m, int cols) {
    double *G = (double *)xcalloc((size_t)m * (size_t)m, sizeof(double));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j <= i; j++) {
            double s = 0.0;
            const double *ri = J + (size_t)i*cols;
            const double *rj = J + (size_t)j*cols;
            for (int k = 0; k < cols; k++) s += ri[k] * rj[k];
            MAT(G,m,i,j) = MAT(G,m,j,i) = s;
        }
    }
    return G;
}

static void jacobi_eigenvalues_symmetric(double *A, int n, double *eig) {
    const int max_sweeps = 100 * n * n + 1000;
    for (int sweep = 0; sweep < max_sweeps; sweep++) {
        int p = 0, q = 1;
        double maxoff = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double a = fabs(MAT(A,n,i,j));
                if (a > maxoff) { maxoff = a; p = i; q = j; }
            }
        }
        double diagmax = 0.0;
        for (int i = 0; i < n; i++) diagmax = fmax(diagmax, fabs(MAT(A,n,i,i)));
        if (maxoff <= 1e-14 * fmax(1.0, diagmax)) break;
        double app = MAT(A,n,p,p), aqq = MAT(A,n,q,q), apq = MAT(A,n,p,q);
        double tau = (aqq - app) / (2.0 * apq);
        double t = copysign(1.0 / (fabs(tau) + sqrt(1.0 + tau*tau)), tau);
        double c = 1.0 / sqrt(1.0 + t*t);
        double s = t * c;
        for (int k = 0; k < n; k++) {
            if (k == p || k == q) continue;
            double akp = MAT(A,n,k,p);
            double akq = MAT(A,n,k,q);
            double nkp = c*akp - s*akq;
            double nkq = s*akp + c*akq;
            MAT(A,n,k,p) = MAT(A,n,p,k) = nkp;
            MAT(A,n,k,q) = MAT(A,n,q,k) = nkq;
        }
        double newpp = c*c*app - 2.0*c*s*apq + s*s*aqq;
        double newqq = s*s*app + 2.0*c*s*apq + c*c*aqq;
        MAT(A,n,p,p) = newpp;
        MAT(A,n,q,q) = newqq;
        MAT(A,n,p,q) = MAT(A,n,q,p) = 0.0;
    }
    for (int i = 0; i < n; i++) eig[i] = MAT(A,n,i,i);
}

static bool cholesky_lower(const double *A, int n, double *C) {
    memset(C, 0, (size_t)n*(size_t)n*sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = MAT(A,n,i,j);
            for (int k = 0; k < j; k++) sum -= MAT(C,n,i,k) * MAT(C,n,j,k);
            if (i == j) {
                if (!(sum > 0.0) || !isfinite(sum)) return false;
                MAT(C,n,i,j) = sqrt(sum);
            } else {
                MAT(C,n,i,j) = sum / MAT(C,n,j,j);
            }
        }
    }
    return true;
}

static bool lower_triangular_inverse(const double *L, int n, double *B) {
    memset(B, 0, (size_t)n*(size_t)n*sizeof(double));
    for (int col = 0; col < n; col++) {
        for (int i = 0; i < n; i++) {
            double rhs = (i == col) ? 1.0 : 0.0;
            for (int k = 0; k < i; k++) rhs -= MAT(L,n,i,k) * MAT(B,n,k,col);
            double diag = MAT(L,n,i,i);
            if (!(fabs(diag) > 0.0) || !isfinite(diag)) return false;
            MAT(B,n,i,col) = rhs / diag;
        }
    }
    return true;
}

static Interval interval_bc_entry(const double *B, const double *C, int n, int i, int j) {
    Interval total = izero();
    for (int k = 0; k < n; k++) total = iadd(total, imul(ipoint(MAT(B,n,i,k)), ipoint(MAT(C,n,k,j))));
    return total;
}

static Interval interval_cct_entry(const double *C, int n, int i, int j) {
    Interval total = izero();
    for (int k = 0; k < n; k++) total = iadd(total, imul(ipoint(MAT(C,n,i,k)), ipoint(MAT(C,n,j,k))));
    return total;
}

static Interval interval_a_entry(const double *J, int m, int cols, double s, int i, int j) {
    Interval total = izero();
    for (int k = 0; k < cols; k++) {
        double x = J[(size_t)i*cols + k], y = J[(size_t)j*cols + k];
        if (x != 0.0 && y != 0.0) total = iadd(total, imul(ipoint(x), ipoint(y)));
    }
    if (i == j) total = isub(total, isquare(ipoint(s)));
    (void)m;
    return total;
}

/* Status-returning Frobenius bound for plain double points. */
static bool frob_upper_points_checked(const double *M, int n, double *out, char *msg, size_t msglen) {
    Interval total = izero();
    for (int i = 0; i < n*n; i++) {
        Interval a = ipoint(dup(fabs(M[i])));
        total = iadd(total, isquare(a));
    }
    Interval root;
    if (!isqrt_interval_checked(clamp_nonneg(total), &root, msg, msglen)) return false;
    *out = root.hi;
    return true;
}

static bool certify_factor_witness(const double *J, int m, int cols, double s, const double *C, char *msg, size_t msglen) {
    double *B = (double *)xcalloc((size_t)m*(size_t)m, sizeof(double));
    if (!lower_triangular_inverse(C, m, B)) {
        snprintf(msg, msglen, "could not invert Cholesky witness");
        free(B); return false;
    }
    Interval *errs = (Interval *)xcalloc((size_t)m*(size_t)m, sizeof(Interval));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            Interval bc = interval_bc_entry(B, C, m, i, j);
            Interval target = ipoint((i == j) ? 1.0 : 0.0);
            errs[(size_t)i*m+j] = isub(target, bc);
        }
    }
    double delta;
    if (!frob_upper_intervals_checked(errs, m*m, &delta, msg, msglen)) {
        free(errs); free(B); return false;
    }
    free(errs);
    if (!(delta < 1.0)) {
        snprintf(msg, msglen, "inverse witness failed: ||I-BC||_F <= %.17g", delta);
        free(B); return false;
    }
    double bnorm;
    if (!frob_upper_points_checked(B, m, &bnorm, msg, msglen)) {
        free(B); return false;
    }
    double factor_sigma_lower = ddown((1.0 - delta) / bnorm);
    if (!(factor_sigma_lower > 0.0)) {
        snprintf(msg, msglen, "factor singular-value bound is not positive");
        free(B); return false;
    }
    Interval *res = (Interval *)xcalloc((size_t)m*(size_t)m, sizeof(Interval));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            Interval aij = interval_a_entry(J, m, cols, s, i, j);
            Interval cij = interval_cct_entry(C, m, i, j);
            res[(size_t)i*m+j] = isub(aij, cij);
        }
    }
    double residual;
    if (!frob_upper_intervals_checked(res, m*m, &residual, msg, msglen)) {
        free(res); free(B); return false;
    }
    free(res);
    double margin = ddown(factor_sigma_lower * factor_sigma_lower - residual);
    if (margin > 0.0) {
        snprintf(msg, msglen, "witness verified: ||I-BC||_F <= %.17g, sigma_min(C) >= %.17g, ||A-CCT||_F <= %.17g, margin >= %.17g", delta, factor_sigma_lower, residual, margin);
        free(B); return true;
    }
    snprintf(msg, msglen, "factor witness margin failed: sigma_min(C)^2 lower %.17g, residual upper %.17g", ddown(factor_sigma_lower*factor_sigma_lower), residual);
    free(B); return false;
}

static bool certify_sigma_lower(const double *J, int m, int cols, double sigma_guess, double shrink, int max_tries, double *sout, char *msg, size_t msglen) {
    double s = ddown(sigma_guess * shrink);
    double *G = compute_jjt(J, m, cols);
    double *A = (double *)xcalloc((size_t)m*(size_t)m, sizeof(double));
    double *C = (double *)xcalloc((size_t)m*(size_t)m, sizeof(double));
    char last[2048] = {0};
    for (int attempt = 0; attempt < max_tries; attempt++) {
        memcpy(A, G, (size_t)m*(size_t)m*sizeof(double));
        for (int i = 0; i < m; i++) MAT(A,m,i,i) -= s*s;
        if (cholesky_lower(A, m, C)) {
            if (certify_factor_witness(J, m, cols, s, C, msg, msglen)) {
                *sout = s;
                free(G); free(A); free(C); return true;
            }
            snprintf(last, sizeof(last), "%s", msg);
        } else {
            snprintf(last, sizeof(last), "ordinary Cholesky failed at s=%.17g", s);
        }
        s = ddown(s * shrink);
    }
    snprintf(msg, msglen, "could not certify sigma lower; last failure: %s", last);
    free(G); free(A); free(C); return false;
}

static bool solve_linear(double *A, double *b, int n, double *x) {
    int N = n;
    double *M = (double *)xcalloc((size_t)N*(size_t)N, sizeof(double));
    double *rhs = (double *)xcalloc((size_t)N, sizeof(double));
    memcpy(M, A, (size_t)N*(size_t)N*sizeof(double));
    memcpy(rhs, b, (size_t)N*sizeof(double));
    for (int k = 0; k < N; k++) {
        int piv = k;
        double best = fabs(MAT(M,N,k,k));
        for (int i = k+1; i < N; i++) {
            double val = fabs(MAT(M,N,i,k));
            if (val > best) { best = val; piv = i; }
        }
        if (!(best > 1e-14)) { free(M); free(rhs); return false; }
        if (piv != k) {
            for (int j = k; j < N; j++) { double tmp = MAT(M,N,k,j); MAT(M,N,k,j)=MAT(M,N,piv,j); MAT(M,N,piv,j)=tmp; }
            double tr = rhs[k]; rhs[k]=rhs[piv]; rhs[piv]=tr;
        }
        for (int i = k+1; i < N; i++) {
            double f = MAT(M,N,i,k) / MAT(M,N,k,k);
            MAT(M,N,i,k) = 0.0;
            for (int j = k+1; j < N; j++) MAT(M,N,i,j) -= f * MAT(M,N,k,j);
            rhs[i] -= f * rhs[k];
        }
    }
    for (int i = N-1; i >= 0; i--) {
        double s = rhs[i];
        for (int j = i+1; j < N; j++) s -= MAT(M,N,i,j) * x[j];
        x[i] = s / MAT(M,N,i,i);
    }
    free(M); free(rhs); return true;
}

static Vec3 affine_combo(const Vec3 *p, const int *idx, const double *w, int n) {
    Vec3 out = {0,0,0};
    for (int i = 0; i < n; i++) out = vecadd(out, vecscale(p[idx[i]], w[i]));
    return out;
}

static void closest_points_active_sets(const Vec3 *verts, const Simplex *A, const Simplex *B, Vec3 *bestp, Vec3 *bestq) {
    double best = INFINITY;
    *bestp = verts[A->v[0]];
    *bestq = verts[B->v[0]];
    int maxmaskA = 1 << A->n, maxmaskB = 1 << B->n;
    for (int ma = 1; ma < maxmaskA; ma++) {
        int ia[3], ca = 0;
        for (int i = 0; i < A->n; i++) if ((ma >> i) & 1) ia[ca++] = A->v[i];
        for (int mb = 1; mb < maxmaskB; mb++) {
            int ib[3], cb = 0;
            for (int i = 0; i < B->n; i++) if ((mb >> i) & 1) ib[cb++] = B->v[i];
            int vars = ca + cb;
            int N = vars + 2;
            double M[64]; double rhs[8]; double sol[8];
            memset(M, 0, sizeof(M)); memset(rhs, 0, sizeof(rhs)); memset(sol, 0, sizeof(sol));
            for (int i = 0; i < ca; i++) for (int j = 0; j < ca; j++) MAT(M,N,i,j) = 2.0 * vecdot(verts[ia[i]], verts[ia[j]]);
            for (int i = 0; i < ca; i++) for (int j = 0; j < cb; j++) MAT(M,N,i,ca+j) = -2.0 * vecdot(verts[ia[i]], verts[ib[j]]);
            for (int i = 0; i < cb; i++) for (int j = 0; j < ca; j++) MAT(M,N,ca+i,j) = -2.0 * vecdot(verts[ib[i]], verts[ia[j]]);
            for (int i = 0; i < cb; i++) for (int j = 0; j < cb; j++) MAT(M,N,ca+i,ca+j) = 2.0 * vecdot(verts[ib[i]], verts[ib[j]]);
            for (int i = 0; i < ca; i++) { MAT(M,N,i,vars) = 1.0; MAT(M,N,vars,i) = 1.0; }
            for (int i = 0; i < cb; i++) { MAT(M,N,ca+i,vars+1) = 1.0; MAT(M,N,vars+1,ca+i) = 1.0; }
            rhs[vars] = 1.0; rhs[vars+1] = 1.0;
            if (!solve_linear(M, rhs, N, sol)) continue;
            bool ok = true;
            for (int i = 0; i < vars; i++) if (sol[i] < -1e-9) ok = false;
            if (!ok) continue;
            double la[3], lb[3], suma = 0.0, sumb = 0.0;
            for (int i = 0; i < ca; i++) { la[i] = fmax(0.0, sol[i]); suma += la[i]; }
            for (int i = 0; i < cb; i++) { lb[i] = fmax(0.0, sol[ca+i]); sumb += lb[i]; }
            if (suma == 0.0 || sumb == 0.0) continue;
            for (int i = 0; i < ca; i++) la[i] /= suma;
            for (int i = 0; i < cb; i++) lb[i] /= sumb;
            Vec3 p = affine_combo(verts, ia, la, ca);
            Vec3 q = affine_combo(verts, ib, lb, cb);
            double d = vecnorm(vecsub(p,q));
            if (d < best) { best = d; *bestp = p; *bestq = q; }
        }
    }
}

static Interval interval_dot_vec(Vec3 n, Vec3 v) {
    Interval total = izero();
    total = iadd(total, imul(ipoint(n.x), ipoint(v.x)));
    total = iadd(total, imul(ipoint(n.y), ipoint(v.y)));
    total = iadd(total, imul(ipoint(n.z), ipoint(v.z)));
    return total;
}

static double norm_upper_vec(Vec3 n) {
    Interval total = izero();
    total = iadd(total, isquare(ipoint(n.x)));
    total = iadd(total, isquare(ipoint(n.y)));
    total = iadd(total, isquare(ipoint(n.z)));
    return isqrt_interval(clamp_nonneg(total)).hi;
}

static double separation_lower_for_normal(const Vec3 *verts, const Simplex *A, const Simplex *B, Vec3 n) {
    double nup = norm_upper_vec(n);
    if (!(nup > 0.0) || !isfinite(nup)) return -INFINITY;
    double minA = INFINITY, maxA = -INFINITY, minB = INFINITY, maxB = -INFINITY;
    for (int i = 0; i < A->n; i++) {
        Interval d = interval_dot_vec(n, verts[A->v[i]]);
        if (d.lo < minA) minA = d.lo;
        if (d.hi > maxA) maxA = d.hi;
    }
    for (int i = 0; i < B->n; i++) {
        Interval d = interval_dot_vec(n, verts[B->v[i]]);
        if (d.lo < minB) minB = d.lo;
        if (d.hi > maxB) maxB = d.hi;
    }
    double gap1 = ddown(minA - maxB);
    double gap2 = ddown(minB - maxA);
    double gap = fmax(gap1, gap2);
    if (gap <= 0.0) return -INFINITY;
    return ddown(gap / nup);
}

static Simplex *build_simplices(int nv, const Edge *edges, int ne, const Face *faces, int nf, int *nsout) {
    int ns = nv + ne + nf;
    Simplex *s = (Simplex *)xcalloc((size_t)ns, sizeof(Simplex));
    int k = 0;
    for (int i = 0; i < nv; i++) { s[k].n = 1; s[k].v[0] = i; k++; }
    for (int i = 0; i < ne; i++) { s[k].n = 2; s[k].v[0] = edges[i].a; s[k].v[1] = edges[i].b; k++; }
    for (int i = 0; i < nf; i++) { s[k].n = 3; s[k].v[0] = faces[i].a; s[k].v[1] = faces[i].b; s[k].v[2] = faces[i].c; k++; }
    *nsout = ns; return s;
}

static bool simplices_disjoint(const Simplex *a, const Simplex *b) {
    for (int i = 0; i < a->n; i++) for (int j = 0; j < b->n; j++) if (a->v[i] == b->v[j]) return false;
    return true;
}

static Vec3 simplex_centroid(const Vec3 *verts, const Simplex *s) {
    Vec3 c = {0,0,0};
    for (int i = 0; i < s->n; i++) c = vecadd(c, verts[s->v[i]]);
    return vecscale(c, 1.0 / (double)s->n);
}

static bool certify_collision_lower(const Vec3 *verts, int nv, const Edge *edges, int ne, const Face *faces, int nf, double *Dout, Simplex *pa, Simplex *pb, int *checkedout) {
    int ns = 0;
    Simplex *simp = build_simplices(nv, edges, ne, faces, nf, &ns);
    double best_lower = INFINITY;
    int checked = 0;
    for (int i = 0; i < ns; i++) {
        for (int j = i+1; j < ns; j++) {
            if (!simplices_disjoint(&simp[i], &simp[j])) continue;
            checked++;
            Vec3 p, q;
            closest_points_active_sets(verts, &simp[i], &simp[j], &p, &q);
            Vec3 candidates[16];
            int nc = 0;
            Vec3 diff = vecsub(p, q);
            if (vecnorm(diff) > 0.0) candidates[nc++] = diff;
            Vec3 cdiff = vecsub(simplex_centroid(verts, &simp[i]), simplex_centroid(verts, &simp[j]));
            if (vecnorm(cdiff) > 0.0 && nc < 16) candidates[nc++] = cdiff;
            for (int a = 0; a < simp[i].n && nc < 16; a++) {
                for (int b = 0; b < simp[j].n && nc < 16; b++) {
                    Vec3 vd = vecsub(verts[simp[i].v[a]], verts[simp[j].v[b]]);
                    if (vecnorm(vd) > 0.0) candidates[nc++] = vd;
                }
            }
            double pair_lower = -INFINITY;
            for (int c = 0; c < nc; c++) {
                double lb = separation_lower_for_normal(verts, &simp[i], &simp[j], candidates[c]);
                if (lb > pair_lower) pair_lower = lb;
            }
            if (!(pair_lower > 0.0)) {
                fprintf(stderr, "failed to certify separation for simplex pair %d,%d\n", i, j);
                free(simp); return false;
            }
            if (pair_lower < best_lower) {
                best_lower = pair_lower;
                *pa = simp[i]; *pb = simp[j];
            }
        }
    }
    *Dout = best_lower; *checkedout = checked;
    free(simp); return true;
}



typedef struct { Interval x, y, z; } IVec3;
typedef struct { int degree; int *nbrs; } Cycle;
typedef struct { int vertex; int degree; double approx; double lower; double upper; bool pass; } TurningRecord;

static IVec3 ivec_point_box(Vec3 p, double r) {
    IVec3 out;
    out.x = (Interval){ ddown(p.x - r), dup(p.x + r) };
    out.y = (Interval){ ddown(p.y - r), dup(p.y + r) };
    out.z = (Interval){ ddown(p.z - r), dup(p.z + r) };
    return out;
}

static IVec3 ivec_sub(IVec3 a, IVec3 b) {
    IVec3 r = { isub(a.x, b.x), isub(a.y, b.y), isub(a.z, b.z) };
    return r;
}

static Interval idot_vec(IVec3 a, IVec3 b) {
    Interval total = izero();
    total = iadd(total, imul(a.x, b.x));
    total = iadd(total, imul(a.y, b.y));
    total = iadd(total, imul(a.z, b.z));
    return total;
}

static IVec3 icross_vec(IVec3 a, IVec3 b) {
    IVec3 r;
    r.x = isub(imul(a.y, b.z), imul(a.z, b.y));
    r.y = isub(imul(a.z, b.x), imul(a.x, b.z));
    r.z = isub(imul(a.x, b.y), imul(a.y, b.x));
    return r;
}

static Interval inorm_vec(IVec3 a) {
    Interval total = izero();
    total = iadd(total, isquare(a.x));
    total = iadd(total, isquare(a.y));
    total = iadd(total, isquare(a.z));
    return clamp_nonneg(isqrt_interval(clamp_nonneg(total)));
}

static bool inormalize_vec(IVec3 a, IVec3 *out, const char *what, char *msg, size_t msglen) {
    Interval n = inorm_vec(a);
    if (!(n.lo > 0.0)) {
        snprintf(msg, msglen, "cannot normalize interval vector for %s: norm interval [%.17g, %.17g]", what, n.lo, n.hi);
        return false;
    }
    out->x = idiv(a.x, n);
    out->y = idiv(a.y, n);
    out->z = idiv(a.z, n);
    return true;
}


static double turning_angle_float(Vec3 uprev, Vec3 ucur, Vec3 unext) {
    Vec3 nprev = vecnormalize(veccross(uprev, ucur));
    Vec3 nnext = vecnormalize(veccross(ucur, unext));
    Vec3 tin = vecnormalize(veccross(nprev, ucur));
    Vec3 tout = vecnormalize(veccross(nnext, ucur));
    double y = vecdot(ucur, veccross(tin, tout));
    double x = vecdot(tin, tout);
    return atan2(y, x);
}

static bool turning_sincos_interval(IVec3 uprev, IVec3 ucur, IVec3 unext, Interval *sout, Interval *cout, char *msg, size_t msglen) {
    IVec3 nprev, nnext, tin, tout;
    if (!inormalize_vec(icross_vec(uprev, ucur), &nprev, "previous great-circle normal", msg, msglen)) return false;
    if (!inormalize_vec(icross_vec(ucur, unext), &nnext, "next great-circle normal", msg, msglen)) return false;
    if (!inormalize_vec(icross_vec(nprev, ucur), &tin, "incoming tangent", msg, msglen)) return false;
    if (!inormalize_vec(icross_vec(nnext, ucur), &tout, "outgoing tangent", msg, msglen)) return false;
    *sout = idot_vec(ucur, icross_vec(tin, tout));
    *cout = idot_vec(tin, tout);
    return true;
}

static void complex_mul_i(Interval ar, Interval ai, Interval br, Interval bi, Interval *rr, Interval *ri) {
    *rr = isub(imul(ar, br), imul(ai, bi));
    *ri = iadd(imul(ar, bi), imul(ai, br));
}

static Cycle *build_vertex_cycles_ccw_outside(int nv, const Face *faces, int nf) {
    int *nxt = (int *)xcalloc((size_t)nv * (size_t)nv, sizeof(int));
    for (int i = 0; i < nv*nv; i++) nxt[i] = -1;
    for (int f = 0; f < nf; f++) {
        int tri[3] = {faces[f].a, faces[f].b, faces[f].c};
        for (int k = 0; k < 3; k++) {
            int v = tri[k];
            int x = tri[(k+1)%3];
            int y = tri[(k+2)%3];
            int tmp = x; x = y; y = tmp; /* fixed convention for ccw-outside faces */
            int *slot = &nxt[(size_t)v*(size_t)nv + (size_t)x];
            if (*slot >= 0 && *slot != y) die("inconsistent oriented link at vertex");
            *slot = y;
        }
    }
    Cycle *cycles = (Cycle *)xcalloc((size_t)nv, sizeof(Cycle));
    for (int v = 0; v < nv; v++) {
        int degree = 0, start = -1;
        for (int j = 0; j < nv; j++) {
            if (nxt[(size_t)v*(size_t)nv + (size_t)j] >= 0) {
                if (start < 0) start = j;
                degree++;
            }
        }
        if (degree == 0) die("isolated vertex in oriented link");
        cycles[v].degree = degree;
        cycles[v].nbrs = (int *)xcalloc((size_t)degree, sizeof(int));
        bool *seen = (bool *)xcalloc((size_t)nv, sizeof(bool));
        int cur = start;
        for (int k = 0; k < degree; k++) {
            if (cur < 0 || cur >= nv) die("broken oriented link");
            if (seen[cur]) die("oriented link is not a single cycle");
            seen[cur] = true;
            cycles[v].nbrs[k] = cur;
            cur = nxt[(size_t)v*(size_t)nv + (size_t)cur];
        }
        if (cur != start) die("oriented link is disconnected");
        free(seen);
    }
    free(nxt);
    return cycles;
}

static void free_cycles(Cycle *cycles, int nv) {
    if (!cycles) return;
    for (int i = 0; i < nv; i++) free(cycles[i].nbrs);
    free(cycles);
}

static double total_turning_float_for_cycle(const Vec3 *verts, int v, const Cycle *cycle) {
    int d = cycle->degree;
    Vec3 *dirs = (Vec3 *)xcalloc((size_t)d, sizeof(Vec3));
    for (int i = 0; i < d; i++) dirs[i] = vecnormalize(vecsub(verts[cycle->nbrs[i]], verts[v]));
    double total = 0.0;
    for (int i = 0; i < d; i++) total += turning_angle_float(dirs[(i+d-1)%d], dirs[i], dirs[(i+1)%d]);
    free(dirs);
    return total;
}

static bool total_turning_sin_half_interval_for_cycle(const IVec3 *boxes, int v, const Cycle *cycle, Interval *out, char *msg, size_t msglen) {
    int d = cycle->degree;
    IVec3 *dirs = (IVec3 *)xcalloc((size_t)d, sizeof(IVec3));
    for (int i = 0; i < d; i++) {
        int a = cycle->nbrs[i];
        if (!inormalize_vec(ivec_sub(boxes[a], boxes[v]), &dirs[i], "vertex direction", msg, msglen)) {
            free(dirs); return false;
        }
    }

    /*
      Each local geodesic turn theta_i gives c_i + I s_i = exp(I theta_i)
      from dot/cross products. Since theta_i is the principal local turn,
      cos(theta_i/2) is positive and equals sqrt((1+c_i)/2). Multiplying
      these half-angle factors gives exp(I T/2), where T is the total
      geodesic turning around the spherical link. Embeddedness is the
      precondition certified externally by scripts/euclid_prover.sh via
      src/embed_check (CGAL exact tri-tri intersection test). With the
      OBJ certified embedded, the spherical link is embedded, so
      -2*pi < T < 2*pi. Hence T > 0 iff sin(T/2) > 0. This proof path
      uses no atan/atan2.
    */
    Interval re = ipoint(1.0);
    Interval im = ipoint(0.0);
    Interval one = ipoint(1.0);
    Interval two = ipoint(2.0);
    for (int i = 0; i < d; i++) {
        Interval s, c;
        if (!turning_sincos_interval(dirs[(i+d-1)%d], dirs[i], dirs[(i+1)%d], &s, &c, msg, msglen)) {
            free(dirs); return false;
        }
        Interval half_cos_sq = idiv(iadd(one, c), two);
        if (!(half_cos_sq.lo > 0.0)) {
            snprintf(msg, msglen, "cannot certify positive half-angle cosine at link corner %d: (1+c)/2 interval=[%.17g, %.17g]", i+1, half_cos_sq.lo, half_cos_sq.hi);
            free(dirs); return false;
        }
        Interval ch = isqrt_interval(half_cos_sq);
        Interval sh = idiv(s, imul(two, ch));
        Interval nr, ni;
        complex_mul_i(re, im, ch, sh, &nr, &ni);
        re = nr; im = ni;
    }
    *out = im;
    free(dirs);
    return true;
}

static bool certify_undented(const Vec3 *verts, int nv, const Face *faces, int nf, double motion_radius, bool print_turning, double *min_lower_out, int *min_vertex_out, char *msg, size_t msglen) {
    Cycle *cycles = build_vertex_cycles_ccw_outside(nv, faces, nf);
    IVec3 *boxes = (IVec3 *)xcalloc((size_t)nv, sizeof(IVec3));
    for (int i = 0; i < nv; i++) boxes[i] = ivec_point_box(verts[i], motion_radius);
    double min_lower = INFINITY;
    int min_vertex = -1;
    bool ok = true;
    char first_failure[1024] = {0};
    for (int v = 0; v < nv; v++) {
        int d = cycles[v].degree;
        Interval sin_half = { -INFINITY, INFINITY };
        char local_msg[512] = {0};
        bool vertex_ok = total_turning_sin_half_interval_for_cycle(boxes, v, &cycles[v], &sin_half, local_msg, sizeof(local_msg));
        double lower = vertex_ok ? sin_half.lo : -INFINITY;
        double upper = vertex_ok ? sin_half.hi : INFINITY;
        bool pass = vertex_ok && lower > 0.0;
        if (lower < min_lower) { min_lower = lower; min_vertex = v + 1; }
        if (print_turning) {
            double approx = NAN;
            /* Diagnostic only. This uses libm atan2, but it is not reached unless
               --print-turning is requested and does not affect certification. */
            approx = total_turning_float_for_cycle(verts, v, &cycles[v]);
            printf("turning_vertex %d: degree=%d approx_total_turning_diagnostic=%.17g sin_half_total_turning_interval=[%.17g, %.17g] status=%s\n", v+1, d, approx, lower, upper, pass ? "PASS" : "FAIL");
        }
        if (!pass && ok) {
            ok = false;
            if (vertex_ok) {
                snprintf(first_failure, sizeof(first_failure), "vertex %d: sin(T/2) lower bound not positive; degree=%d, sin_half_interval=[%.17g, %.17g]", v+1, d, lower, upper);
            } else {
                snprintf(first_failure, sizeof(first_failure), "vertex %d: sin(T/2) interval evaluation failed: %s", v+1, local_msg);
            }
        }
    }
    *min_lower_out = min_lower;
    *min_vertex_out = min_vertex;
    if (ok) snprintf(msg, msglen, "all %d vertex sin(T/2) lower bounds are positive under ccw-outside face convention", nv);
    else snprintf(msg, msglen, "%s", first_failure);
    free(boxes);
    free_cycles(cycles, nv);
    return ok;
}

static double upper_motion_bound(double sigma_lower, double rho_upper, int edge_count) {
    double sqrtE = dup(sqrt((double)edge_count));
    double disc = ddown(sigma_lower * sigma_lower - dup(16.0 * rho_upper * sqrtE));
    if (!(disc > 0.0)) return INFINITY;
    double root = ddown(sqrt(disc));
    double numerator = dup(sigma_lower - root);
    double denominator = ddown(8.0 * sqrtE);
    return dup(numerator / denominator);
}

static void print_simplex(Simplex s) {
    printf("(");
    for (int i = 0; i < s.n; i++) {
        if (i) printf(",");
        printf("%d", s.v[i]+1);
    }
    printf(")");
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s input.obj [--sigma-shrink x] [--max-sigma-tries n] [--sharp-collision] [--print-turning]\n", argv[0]);
        return 1;
    }
    const char *path = argv[1];
    double shrink = 0.98;
    int max_tries = 80;
    bool sharp_collision = false;
    bool print_turning = false;
    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "--sigma-shrink") == 0 && i+1 < argc) shrink = atof(argv[++i]);
        else if (strcmp(argv[i], "--max-sigma-tries") == 0 && i+1 < argc) max_tries = atoi(argv[++i]);
        else if (strcmp(argv[i], "--sharp-collision") == 0) sharp_collision = true;
        else if (strcmp(argv[i], "--print-turning") == 0) print_turning = true;
        else die("unknown command-line option");
    }

    Vec3 *verts = NULL; Face *faces = NULL;
    int nv = 0, nf = 0;
    parse_obj(path, &verts, &nv, &faces, &nf);
    int ne = 0;
    Edge *edges = edges_from_faces(faces, nf, &ne);
    int cols = 0;
    double *J = build_jacobian(verts, nv, edges, ne, &cols);

    double rho_upper = certify_rho_upper(verts, edges, ne);

    double *G = compute_jjt(J, ne, cols);
    double *Gcopy = (double *)xcalloc((size_t)ne*(size_t)ne, sizeof(double));
    memcpy(Gcopy, G, (size_t)ne*(size_t)ne*sizeof(double));
    double *eig = (double *)xcalloc((size_t)ne, sizeof(double));
    jacobi_eigenvalues_symmetric(Gcopy, ne, eig);
    double mineig = INFINITY;
    for (int i = 0; i < ne; i++) if (eig[i] < mineig) mineig = eig[i];
    double sigma_guess = sqrt(fmax(0.0, mineig));
    free(G); free(Gcopy); free(eig);

    double sigma_lower = 0.0;
    char sigma_msg[1024];
    bool sigma_ok = certify_sigma_lower(J, ne, cols, sigma_guess, shrink, max_tries, &sigma_lower, sigma_msg, sizeof(sigma_msg));
    if (!sigma_ok) {
        printf("vertices: %d\n", nv);
        printf("faces: %d\n", nf);
        printf("edges: %d\n", ne);
        printf("rho_upper: %.17g\n", rho_upper);
        printf("sigma_guess_uncertified: %.17g\n", sigma_guess);
        printf("sigma_method: cholesky-witness\n");
        printf("sigma_certificate: FAILED: %s\n", sigma_msg);
        printf("final: REJECT\n");
        free(verts); free(faces); free(edges); free(J);
        return 2;
    }

    double collision_lower = 0.0;
    Simplex ca = {0,{0,0,0}}, cb = {0,{0,0,0}};
    int checked = 0;
    bool collision_ok = certify_collision_lower(verts, nv, edges, ne, faces, nf, &collision_lower, &ca, &cb, &checked);
    if (!collision_ok) {
        printf("vertices: %d\n", nv);
        printf("faces: %d\n", nf);
        printf("edges: %d\n", ne);
        printf("rho_upper: %.17g\n", rho_upper);
        printf("sigma_lower_certified: %.17g\n", sigma_lower);
        printf("collision_certificate: FAILED\n");
        printf("final: REJECT\n");
        free(verts); free(faces); free(edges); free(J);
        return 2;
    }

    double sqrtE = dup(sqrt((double)ne));
    double existence_rhs = ddown((sigma_lower * sigma_lower) / dup(16.0 * sqrtE));
    bool existence_ok = rho_upper < existence_rhs;
    double motion_upper = upper_motion_bound(sigma_lower, rho_upper, ne);

    double turning_min_lower = -INFINITY;
    int turning_min_vertex = -1;
    char turning_msg[1024];
    bool undented_ok = certify_undented(verts, nv, faces, nf, motion_upper, print_turning, &turning_min_lower, &turning_min_vertex, turning_msg, sizeof(turning_msg));

    double two_motion_value = dup(2.0 * motion_upper);
    bool two_motion_ok = two_motion_value < collision_lower;
    double ellison_value = dup(dup(sqrt((double)nv)) * motion_upper);
    bool ellison_ok = ellison_value < collision_lower;
    bool embedding_ok = sharp_collision ? two_motion_ok : ellison_ok;
    bool accepted = existence_ok && embedding_ok && undented_ok;

    printf("vertices: %d\n", nv);
    printf("faces: %d\n", nf);
    printf("edges: %d\n", ne);
    printf("rho_upper: %.17g\n", rho_upper);
    printf("sigma_guess_uncertified: %.17g\n", sigma_guess);
    printf("sigma_method: cholesky-witness\n");
    printf("sigma_lower_certified: %.17g\n", sigma_lower);
    printf("sigma_certificate: %s\n", sigma_msg);
    printf("existence_rhs_sigma2_over_16sqrtE: %.17g\n", existence_rhs);
    printf("existence_test: %s\n", existence_ok ? "PASS" : "FAIL");
    printf("motion_upper: %.17g\n", motion_upper);
    printf("collision_lower: %.17g\n", collision_lower);
    printf("closest_certified_pair: "); print_simplex(ca); printf(" "); print_simplex(cb); printf("\n");
    printf("collision_certificate: checked %d disjoint simplex pairs\n", checked);
    printf("two_motion_test: %s  value=%.17g < %.17g\n", two_motion_ok ? "PASS" : "FAIL", two_motion_value, collision_lower);
    printf("ellison_sqrtV_motion_test: %s  value=%.17g < %.17g\n", ellison_ok ? "PASS" : "FAIL", ellison_value, collision_lower);
    printf("embedding_test_used: %s\n", sharp_collision ? "sharp_two_motion" : "ellison_sqrtV");
    printf("turning_convention: ccw_outside_faces\n");
    printf("turning_certificate: sin_half_total_turning_positive_using_embedding_bound_-2pi_lt_T_lt_2pi\n");
    printf("sin_half_total_turning_min_lower: %.17g\n", turning_min_lower);
    printf("turning_min_vertex: %d\n", turning_min_vertex);
    printf("undented_test: %s  %s\n", undented_ok ? "PASS" : "FAIL", turning_msg);
    printf("final: %s\n", accepted ? "ACCEPT" : "REJECT");

    free(verts); free(faces); free(edges); free(J);
    return accepted ? 0 : 2;
}
