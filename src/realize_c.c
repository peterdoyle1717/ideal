/*
 * realize_c.c — Klein-OBJ realization from a puffup bends file.
 *
 * Reads a single .bends file on stdin (the format puffup_c writes with
 * --bends-out):
 *
 *     puffup-bends 1
 *     NV V NE E alpha_deg A
 *     faces a,b,c;d,e,f;...
 *     bends
 *     u v θ
 *     ...   (NE lines, canonical u<v ordering)
 *
 * For an equilateral hyperbolic triangulation with corner angle α,
 * side length s satisfies cosh(s) = cos α / (1 − cos α). We:
 *   1. Build the same canonical edge/face index as puffup_c (so that
 *      the bend on edge (u,v) lines up with the geometric edge).
 *   2. Place FACES[0] in the hyperboloid model H³ ⊂ R^{3,1} with one
 *      vertex at the apex (0,0,0,1) and the other two boosted out by
 *      hyperbolic distance s.
 *   3. BFS over the dual graph; for each new face we cross edge (a,b)
 *      from a placed face whose third vertex is p, and place the new
 *      third vertex c via the Lorentz analog of puffup's place_third:
 *        c = M cosh(h) + sinh(h) · (−τ_n cos θ + τ_w sin θ)
 *      where M is the hyperbolic midpoint of (a,b), τ_n is the unit
 *      tangent at M in the (a,b,p) face plane pointing toward p,
 *      τ_w completes a right-handed frame {τ_e, τ_n, τ_w} in the
 *      tangent 3-space at M (with τ_e along the edge), θ is the bend
 *      from the bends file (puffup convention: bend = π − dihedral),
 *      and h is the equilateral altitude with cosh(h) = cosh(s)/cosh(s/2).
 *   4. Project to Klein: write `v xi/x4 yi/x4 zi/x4` per vertex.
 *
 * Compile: cc -O3 -o realize_c realize_c.c -lm
 * Usage:   realize_c < net.bends > net.obj
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXV    1200
#define MAXF    (2*MAXV + 4)
#define MAXE    (3*MAXV - 6)
#define MAXLINE 65536

/* ---------- combinatorial state -------------------------------------------- */
typedef struct { int a, b, c; } Face;

static int  NV, NF, NE;
static Face FACES[MAXF];
static int  EDGE_A[MAXE], EDGE_B[MAXE];
static int  EDGE_IDX[MAXV+1][MAXV+1];   /* canonical edge index, -1 if none */
static int  EM_F[MAXV+1][MAXV+1];       /* face containing directed edge (u,v) */
static double BENDS[MAXE];
static double ALPHA_DEG;

/* ---------- placement state ------------------------------------------------ */
static double V[MAXV+1][4];             /* hyperboloid coords (x1,x2,x3,x4) */
static int    V_PLACED[MAXV+1];

/* ---------- Minkowski arithmetic ------------------------------------------- */
/* Signature (+,+,+,−). Hyperboloid: <x,x> = −1, x4 > 0. */

static double mdot(const double a[4], const double b[4]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] - a[3]*b[3];
}
static void mscale(double s, const double a[4], double r[4]) {
    r[0]=s*a[0]; r[1]=s*a[1]; r[2]=s*a[2]; r[3]=s*a[3];
}
static void madd(const double a[4], const double b[4], double r[4]) {
    r[0]=a[0]+b[0]; r[1]=a[1]+b[1]; r[2]=a[2]+b[2]; r[3]=a[3]+b[3];
}
static void msub(const double a[4], const double b[4], double r[4]) {
    r[0]=a[0]-b[0]; r[1]=a[1]-b[1]; r[2]=a[2]-b[2]; r[3]=a[3]-b[3];
}

/* Project v onto the tangent space at M (subtract its M-component).
   Uses <M,M> = −1, so coefficient of M is +<v,M>. */
static void tang_at(const double M[4], const double v[4], double r[4]) {
    double k = mdot(v, M);
    r[0] = v[0] + k*M[0];
    r[1] = v[1] + k*M[1];
    r[2] = v[2] + k*M[2];
    r[3] = v[3] + k*M[3];
}

/* Subtract the spacelike-unit-vector u component from v (Gram-Schmidt). */
static void orth_against(const double u[4], double v[4]) {
    double k = mdot(v, u);
    v[0] -= k*u[0]; v[1] -= k*u[1]; v[2] -= k*u[2]; v[3] -= k*u[3];
}

static int unit_spacelike(double v[4]) {
    double n2 = mdot(v, v);
    if (!(n2 > 1e-24)) return -1;
    double s = 1.0/sqrt(n2);
    v[0]*=s; v[1]*=s; v[2]*=s; v[3]*=s;
    return 0;
}

/* Find the unique unit-spacelike tangent at M orthogonal to (te, tn).
   Tries each Minkowski basis vector in turn, projects, normalizes. */
static int fourth_tangent(const double M[4], const double te[4],
                          const double tn[4], double tw[4]) {
    for (int k = 0; k < 4; k++) {
        double cand[4] = {0,0,0,0};
        cand[k] = 1.0;
        tang_at(M, cand, cand);
        orth_against(te, cand);
        orth_against(tn, cand);
        if (unit_spacelike(cand) == 0) {
            tw[0]=cand[0]; tw[1]=cand[1]; tw[2]=cand[2]; tw[3]=cand[3];
            return 0;
        }
    }
    return -1;
}

/* Set the orientation of tw so that {M, te, tn, tw} is right-handed,
   i.e. det(M | te | tn | tw) > 0 with rows treated as columns. */
static void orient_tw(const double M[4], const double te[4],
                      const double tn[4], double tw[4]) {
    double d =
        + M[0]*(te[1]*(tn[2]*tw[3]-tn[3]*tw[2])
              - te[2]*(tn[1]*tw[3]-tn[3]*tw[1])
              + te[3]*(tn[1]*tw[2]-tn[2]*tw[1]))
        - M[1]*(te[0]*(tn[2]*tw[3]-tn[3]*tw[2])
              - te[2]*(tn[0]*tw[3]-tn[3]*tw[0])
              + te[3]*(tn[0]*tw[2]-tn[2]*tw[0]))
        + M[2]*(te[0]*(tn[1]*tw[3]-tn[3]*tw[1])
              - te[1]*(tn[0]*tw[3]-tn[3]*tw[0])
              + te[3]*(tn[0]*tw[1]-tn[1]*tw[0]))
        - M[3]*(te[0]*(tn[1]*tw[2]-tn[2]*tw[1])
              - te[1]*(tn[0]*tw[2]-tn[2]*tw[0])
              + te[2]*(tn[0]*tw[1]-tn[1]*tw[0]));
    if (d < 0) {
        tw[0]=-tw[0]; tw[1]=-tw[1]; tw[2]=-tw[2]; tw[3]=-tw[3];
    }
}

/* ---------- bends-file parser ---------------------------------------------- */
static int parse_facelist_into_faces(const char *line) {
    NF = 0; NV = 0;
    const char *p = line;
    while (*p) {
        int a, b, c;
        if (sscanf(p, "%d,%d,%d", &a, &b, &c) != 3) break;
        if (NF >= MAXF || a > MAXV || b > MAXV || c > MAXV) return -1;
        FACES[NF++] = (Face){a, b, c};
        if (a > NV) NV = a;
        if (b > NV) NV = b;
        if (c > NV) NV = c;
        while (*p && *p != ';') p++;
        if (*p == ';') p++;
    }
    return NF;
}

static int read_bends_file(FILE *fh) {
    char line[MAXLINE];
    int nv=-1, ne=-1;
    int saw_header=0, saw_dim=0, saw_faces=0, saw_bends_marker=0;
    int n_bends_read = 0;

    while (fgets(line, sizeof(line), fh)) {
        int ll = strlen(line);
        while (ll > 0 && (line[ll-1]=='\n' || line[ll-1]=='\r')) line[--ll] = '\0';
        if (!ll) continue;

        if (!saw_header) {
            int ver;
            if (sscanf(line, "puffup-bends %d", &ver) != 1) {
                fprintf(stderr, "ERROR: expected 'puffup-bends N' header\n");
                return -1;
            }
            saw_header = 1;
            continue;
        }
        if (!saw_dim) {
            if (sscanf(line, "NV %d NE %d alpha_deg %lf", &nv, &ne, &ALPHA_DEG) != 3) {
                fprintf(stderr, "ERROR: expected 'NV V NE E alpha_deg A'\n");
                return -1;
            }
            saw_dim = 1;
            continue;
        }
        if (!saw_faces) {
            const char *q = line;
            while (*q == ' ') q++;
            if (strncmp(q, "faces ", 6) != 0) {
                fprintf(stderr, "ERROR: expected 'faces ...' line\n"); return -1;
            }
            if (parse_facelist_into_faces(q + 6) <= 0) {
                fprintf(stderr, "ERROR: empty/invalid face list\n"); return -1;
            }
            saw_faces = 1;
            continue;
        }
        if (!saw_bends_marker) {
            if (strcmp(line, "bends") != 0) {
                fprintf(stderr, "ERROR: expected 'bends' marker\n"); return -1;
            }
            saw_bends_marker = 1;
            continue;
        }
        int u, v;
        double bend;
        if (sscanf(line, "%d %d %lf", &u, &v, &bend) != 3) {
            fprintf(stderr, "ERROR: bad bend line: %s\n", line); return -1;
        }
        if (n_bends_read >= MAXE || u < 1 || v < 1 || u > MAXV || v > MAXV) {
            fprintf(stderr, "ERROR: bend index out of range\n"); return -1;
        }
        if (u >= v) {
            fprintf(stderr, "ERROR: non-canonical edge (u<v required): %d %d\n", u, v);
            return -1;
        }
        EDGE_A[n_bends_read] = u;
        EDGE_B[n_bends_read] = v;
        BENDS[n_bends_read] = bend;
        n_bends_read++;
    }
    if (!saw_bends_marker || n_bends_read != ne) {
        fprintf(stderr, "ERROR: read %d bends, expected NE=%d\n", n_bends_read, ne);
        return -1;
    }
    NE = ne;
    if (NV != nv) {
        fprintf(stderr, "ERROR: header NV=%d disagrees with face list NV=%d\n", nv, NV);
        return -1;
    }
    return 0;
}

/* ---------- index building (must match puffup_c canonical order) ----------- */
static int build_indices(void) {
    /* Edge index map. Walk faces, register each unordered edge once. */
    for (int u = 0; u <= MAXV; u++)
        for (int v = 0; v <= MAXV; v++) { EDGE_IDX[u][v] = -1; EM_F[u][v] = 0; }

    for (int i = 0; i < NF; i++) {
        int a = FACES[i].a, b = FACES[i].b, c = FACES[i].c;
        EM_F[a][b] = i;
        EM_F[b][c] = i;
        EM_F[c][a] = i;
    }
    /* The bends file already supplies (EDGE_A[i], EDGE_B[i]) per edge in the
       same order puffup wrote. Just install the inverse map. */
    for (int i = 0; i < NE; i++) {
        EDGE_IDX[EDGE_A[i]][EDGE_B[i]] = i;
        EDGE_IDX[EDGE_B[i]][EDGE_A[i]] = i;
    }
    return 0;
}

/* ---------- placement ------------------------------------------------------ */
static double S_LEN, COSH_S, SINH_S, COSH_HALF_S;
static double H_ALT, COSH_H, SINH_H;

static int compute_geom(void) {
    double a_rad = ALPHA_DEG * M_PI / 180.0;
    double ca = cos(a_rad);
    if (!(ca > 0 && ca < 1.0)) {
        fprintf(stderr, "ERROR: alpha must be in (0°, 90°), got %g°\n", ALPHA_DEG);
        return -1;
    }
    COSH_S = ca / (1.0 - ca);
    if (COSH_S <= 1.0) {
        fprintf(stderr, "ERROR: cosh(s)≤1 at α=%g°\n", ALPHA_DEG); return -1;
    }
    S_LEN = acosh(COSH_S);
    SINH_S = sqrt(COSH_S*COSH_S - 1.0);
    COSH_HALF_S = cosh(0.5 * S_LEN);
    COSH_H = COSH_S / COSH_HALF_S;
    if (COSH_H < 1.0) COSH_H = 1.0;
    H_ALT = acosh(COSH_H);
    SINH_H = sqrt(COSH_H*COSH_H - 1.0);
    return 0;
}

static void place_base(void) {
    /* FACES[0] = (a, b, c). Place:
       a at apex (0, 0, 0, 1)
       b at boost(s) along +x:        (sinh s, 0, 0, cosh s)
       c at boost(s) at angle α:      (sinh s · cosα, sinh s · sinα, 0, cosh s)
    */
    int a = FACES[0].a, b = FACES[0].b, c = FACES[0].c;
    V[a][0]=0;             V[a][1]=0;             V[a][2]=0; V[a][3]=1;
    V[b][0]=SINH_S;        V[b][1]=0;             V[b][2]=0; V[b][3]=COSH_S;
    double a_rad = ALPHA_DEG * M_PI / 180.0;
    V[c][0]=SINH_S*cos(a_rad);
    V[c][1]=SINH_S*sin(a_rad);
    V[c][2]=0; V[c][3]=COSH_S;
    V_PLACED[a] = V_PLACED[b] = V_PLACED[c] = 1;
}

/* Lorentz analog of puffup's place_third. Given placed a, b and the
   previous-face third vertex p, place c with bend θ on edge (a,b). */
static int lorentz_place_third(int ai, int bi, int pi, double theta, int ci) {
    double M[4], te[4], tn[4], tw[4], cand[4];
    /* Hyperbolic midpoint M = normalize(a + b) */
    madd(V[ai], V[bi], cand);
    double m2 = -mdot(cand, cand);
    if (!(m2 > 1e-24)) {
        fprintf(stderr, "ERROR: degenerate edge %d-%d\n", ai, bi);
        return -1;
    }
    mscale(1.0/sqrt(m2), cand, M);

    /* τ_e = unit-spacelike along (b − a), tangent at M */
    msub(V[bi], V[ai], te);
    /* (b − a) is automatically tangent at M since <a, M> = <b, M>. Normalize. */
    if (unit_spacelike(te) < 0) {
        fprintf(stderr, "ERROR: edge tangent degenerate\n"); return -1;
    }

    /* τ_n = unit-spacelike at M in face plane, pointing toward p */
    msub(V[pi], M, tn);
    tang_at(M, tn, tn);
    orth_against(te, tn);
    if (unit_spacelike(tn) < 0) {
        fprintf(stderr, "ERROR: face tangent degenerate\n"); return -1;
    }

    /* τ_w = unique unit-spacelike tangent perp to {M, τ_e, τ_n} */
    if (fourth_tangent(M, te, tn, tw) < 0) {
        fprintf(stderr, "ERROR: cannot find 4th tangent\n"); return -1;
    }
    orient_tw(M, te, tn, tw);

    /* c = M · cosh(h) + sinh(h) · (−τ_n cosθ + τ_w sinθ) */
    double ct = cos(theta), st = sin(theta);
    double dir[4];
    dir[0] = -ct*tn[0] + st*tw[0];
    dir[1] = -ct*tn[1] + st*tw[1];
    dir[2] = -ct*tn[2] + st*tw[2];
    dir[3] = -ct*tn[3] + st*tw[3];
    V[ci][0] = COSH_H*M[0] + SINH_H*dir[0];
    V[ci][1] = COSH_H*M[1] + SINH_H*dir[1];
    V[ci][2] = COSH_H*M[2] + SINH_H*dir[2];
    V[ci][3] = COSH_H*M[3] + SINH_H*dir[3];
    V_PLACED[ci] = 1;
    return 0;
}

static int bfs_place(void) {
    static int q[MAXF], face_placed[MAXF];
    memset(face_placed, 0, NF*sizeof(int));
    memset(V_PLACED, 0, (NV+2)*sizeof(int));

    place_base();
    face_placed[0] = 1;
    int qh = 0, qt = 0;
    q[qt++] = 0;

    while (qh < qt) {
        int fi = q[qh++];
        int verts[3] = {FACES[fi].a, FACES[fi].b, FACES[fi].c};
        for (int i = 0; i < 3; i++) {
            int a = verts[i], b = verts[(i+1)%3];
            int other_fi = EM_F[b][a];      /* face containing the reverse edge */
            if (other_fi == fi || face_placed[other_fi]) continue;
            int oa = FACES[other_fi].a, ob = FACES[other_fi].b, oc = FACES[other_fi].c;
            int c = (oa!=a && oa!=b) ? oa : ((ob!=a && ob!=b) ? ob : oc);
            int p = (verts[0]!=a && verts[0]!=b) ? verts[0]
                  : ((verts[1]!=a && verts[1]!=b) ? verts[1] : verts[2]);
            int e = EDGE_IDX[a][b];
            if (e < 0) {
                fprintf(stderr, "ERROR: missing edge index for (%d,%d)\n", a, b);
                return -1;
            }
            if (lorentz_place_third(a, b, p, BENDS[e], c) < 0) return -1;
            face_placed[other_fi] = 1;
            q[qt++] = other_fi;
        }
    }
    for (int v = 1; v <= NV; v++) {
        if (!V_PLACED[v]) {
            fprintf(stderr, "ERROR: vertex %d not placed (dual graph not connected?)\n", v);
            return -1;
        }
    }
    return 0;
}

/* ---------- output (Klein projection) -------------------------------------- */
static void write_klein_obj(FILE *fh) {
    fprintf(fh, "# Klein realization — α=%g°, side s=%.10g\n", ALPHA_DEG, S_LEN);
    for (int v = 1; v <= NV; v++) {
        double w = V[v][3];
        fprintf(fh, "v %.10f %.10f %.10f\n",
                V[v][0]/w, V[v][1]/w, V[v][2]/w);
    }
    for (int i = 0; i < NF; i++) {
        fprintf(fh, "f %d %d %d\n", FACES[i].a, FACES[i].b, FACES[i].c);
    }
}

/* ---------- main ----------------------------------------------------------- */
int main(int argc, char **argv) {
    (void)argc; (void)argv;
    if (read_bends_file(stdin) < 0) return 1;
    if (build_indices() < 0)        return 1;
    if (compute_geom() < 0)         return 1;
    if (bfs_place() < 0)            return 1;
    write_klein_obj(stdout);
    return 0;
}
