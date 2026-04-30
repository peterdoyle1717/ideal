/*
 * realize_c.c — Klein-OBJ realization from a puffup bends file. (v2)
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
 * Each face F carries a 4×4 SO(3,1) frame F_M and a cyclic labeling
 * (V0, V1, V2). Vertices sit at canonical local positions
 *
 *     canon[0] = (0, 0, 0, 1)                                   (apex)
 *     canon[1] = (sinh s, 0, 0, cosh s)                         (boost +x by s)
 *     canon[2] = (sinh s · cos α, sinh s · sin α, 0, cosh s)    (matz(α) · canon[1])
 *
 * with s = arccosh(cos α / (1 − cos α)), and the world position of F's
 * i-th labeled vertex is F_M · canon[i]. FACES[0] is the BFS root with
 * F_M = identity.
 *
 * BFS step: when face F (labels (V0, V1, V2), frame F_M) crosses its
 * canon-edge i (the edge between V_i and V_{(i+1) mod 3}) into neighbor N,
 * with bend β read from the bends file:
 *
 *     F_N_M = F_M · R[i] · boost(+x, s) · matz(π) · matx(−β)
 *
 * where R[i] = R_1^i is the i-th power of the 3-cycle Lorentz isometry
 * that fixes the canonical-triangle centroid and rotates canon[0] →
 * canon[1] → canon[2] → canon[0]. R_1 is precomputed from α as
 *
 *     R_1 = T_G · matz(2π/3) · T_G^{-1}
 *
 * with T_G the boost from the apex to the triangle centroid G. N's
 * labels are (V_{(i+1) mod 3}, V_i, c′) and N's third-vertex world
 * position is F_N_M · canon[2].
 *
 * No Gram-Schmidt, no orientation determinant. Sign and chirality are
 * fixed by the composition. Choosing matx(−β) (rather than matx(+β))
 * makes the per-vertex induced tangent step at the apex equal the
 * solver's movemat(α, β) = matz(α) · matx(−β) — same convention, no
 * D-conjugation.
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

/* per-face frame and labeling (BFS state) */
typedef double M4[4][4];
static M4  FRAMES[MAXF];
static int LABELS[MAXF][3];
static int FACE_PLACED[MAXF];

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
    for (int u = 0; u <= MAXV; u++)
        for (int v = 0; v <= MAXV; v++) { EDGE_IDX[u][v] = -1; EM_F[u][v] = -1; }

    for (int i = 0; i < NF; i++) {
        int a = FACES[i].a, b = FACES[i].b, c = FACES[i].c;
        EM_F[a][b] = i;
        EM_F[b][c] = i;
        EM_F[c][a] = i;
    }
    for (int i = 0; i < NE; i++) {
        EDGE_IDX[EDGE_A[i]][EDGE_B[i]] = i;
        EDGE_IDX[EDGE_B[i]][EDGE_A[i]] = i;
    }
    return 0;
}

/* ---------- 4×4 Lorentz arithmetic ----------------------------------------- */
static void mat4_iden(M4 M) {
    memset(M, 0, sizeof(M4));
    for (int i = 0; i < 4; i++) M[i][i] = 1.0;
}
static void mat4_mul(const M4 A, const M4 B, M4 C) {
    M4 T;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
            double s = 0.0;
            for (int k = 0; k < 4; k++) s += A[i][k] * B[k][j];
            T[i][j] = s;
        }
    memcpy(C, T, sizeof(M4));
}
static void mat4_apply(const M4 A, const double v[4], double r[4]) {
    double t[4];
    for (int i = 0; i < 4; i++) {
        t[i] = A[i][0]*v[0] + A[i][1]*v[1] + A[i][2]*v[2] + A[i][3]*v[3];
    }
    memcpy(r, t, 4*sizeof(double));
}

/* boost(+x, d): rapidity d along +x. Mixes x and t. */
static void mat4_boost_x(double d, M4 M) {
    double cd = cosh(d), sd = sinh(d);
    mat4_iden(M);
    M[0][0] = cd; M[0][3] = sd;
    M[3][0] = sd; M[3][3] = cd;
}
/* matz(theta): rotation around z-axis by theta (mixes x and y). */
static void mat4_rotz(double theta, M4 M) {
    double c = cos(theta), s = sin(theta);
    mat4_iden(M);
    M[0][0] = c;  M[0][1] = -s;
    M[1][0] = s;  M[1][1] =  c;
}
/* matx(theta): rotation around x-axis by theta (mixes y and z). */
static void mat4_rotx(double theta, M4 M) {
    double c = cos(theta), s = sin(theta);
    mat4_iden(M);
    M[1][1] = c;  M[1][2] = -s;
    M[2][1] = s;  M[2][2] =  c;
}
/* boost in arbitrary spatial direction (ux, uy, uz unit) by rapidity d. */
static void mat4_boost_dir(double ux, double uy, double uz, double d, M4 M) {
    double cd = cosh(d), sd = sinh(d);
    mat4_iden(M);
    M[0][0] = 1.0 + (cd - 1.0)*ux*ux;
    M[0][1] =       (cd - 1.0)*ux*uy;
    M[0][2] =       (cd - 1.0)*ux*uz;
    M[1][0] =       (cd - 1.0)*uy*ux;
    M[1][1] = 1.0 + (cd - 1.0)*uy*uy;
    M[1][2] =       (cd - 1.0)*uy*uz;
    M[2][0] =       (cd - 1.0)*uz*ux;
    M[2][1] =       (cd - 1.0)*uz*uy;
    M[2][2] = 1.0 + (cd - 1.0)*uz*uz;
    M[0][3] = sd*ux;  M[3][0] = sd*ux;
    M[1][3] = sd*uy;  M[3][1] = sd*uy;
    M[2][3] = sd*uz;  M[3][2] = sd*uz;
    M[3][3] = cd;
}

/* ---------- geometry constants & precomputed operators --------------------- */
static double S_LEN, COSH_S, SINH_S;
static double ALPHA_RAD;
static double CANON[3][4];        /* canon[0..2] in 4D */
static M4     R_POW[3];           /* R_POW[i] = R_1^i, i = 0,1,2 */
static M4     BOOST_X_S, MATZ_PI; /* precomputed common factors */

static int compute_geom(void) {
    ALPHA_RAD = ALPHA_DEG * M_PI / 180.0;
    double ca = cos(ALPHA_RAD);
    if (!(ca > 0 && ca < 1.0)) {
        fprintf(stderr, "ERROR: alpha must be in (0°, 90°), got %g°\n", ALPHA_DEG);
        return -1;
    }
    COSH_S = ca / (1.0 - ca);
    if (COSH_S <= 1.0) {
        fprintf(stderr, "ERROR: cosh(s) ≤ 1 at α=%g°\n", ALPHA_DEG); return -1;
    }
    S_LEN = acosh(COSH_S);
    SINH_S = sqrt(COSH_S*COSH_S - 1.0);

    CANON[0][0] = 0.0;          CANON[0][1] = 0.0;
    CANON[0][2] = 0.0;          CANON[0][3] = 1.0;
    CANON[1][0] = SINH_S;       CANON[1][1] = 0.0;
    CANON[1][2] = 0.0;          CANON[1][3] = COSH_S;
    CANON[2][0] = SINH_S*ca;    CANON[2][1] = SINH_S*sin(ALPHA_RAD);
    CANON[2][2] = 0.0;          CANON[2][3] = COSH_S;

    /* R_1: 3-cycle around triangle centroid G.
       G = (canon[0]+canon[1]+canon[2]) normalized to hyperboloid;
       R_1 = T_G · matz(2π/3) · T_G^{-1}. */
    double G_un[4] = {
        CANON[0][0]+CANON[1][0]+CANON[2][0],
        CANON[0][1]+CANON[1][1]+CANON[2][1],
        CANON[0][2]+CANON[1][2]+CANON[2][2],
        CANON[0][3]+CANON[1][3]+CANON[2][3]
    };
    double minkowski_norm_sq = -(G_un[0]*G_un[0] + G_un[1]*G_un[1]
                                + G_un[2]*G_un[2]) + G_un[3]*G_un[3];
    double Gnorm = sqrt(minkowski_norm_sq);
    double G[4] = { G_un[0]/Gnorm, G_un[1]/Gnorm, G_un[2]/Gnorm, G_un[3]/Gnorm };

    double G_r = sqrt(G[0]*G[0] + G[1]*G[1] + G[2]*G[2]);
    double dG = (G_r > 1e-15) ? asinh(G_r) : 0.0;
    M4 T_G, T_G_inv, RZ23, T1, R1;
    if (G_r > 1e-15) {
        mat4_boost_dir(G[0]/G_r, G[1]/G_r, G[2]/G_r, dG, T_G);
        mat4_boost_dir(G[0]/G_r, G[1]/G_r, G[2]/G_r, -dG, T_G_inv);
    } else {
        mat4_iden(T_G); mat4_iden(T_G_inv);
    }
    mat4_rotz(2.0*M_PI/3.0, RZ23);
    mat4_mul(T_G, RZ23, T1);
    mat4_mul(T1, T_G_inv, R1);

    mat4_iden(R_POW[0]);
    memcpy(R_POW[1], R1, sizeof(M4));
    mat4_mul(R1, R1, R_POW[2]);

    mat4_boost_x(S_LEN, BOOST_X_S);
    mat4_rotz(M_PI, MATZ_PI);

    return 0;
}

/* T_cross_i(β) = R_POW[i] · BOOST_X_S · MATZ_PI · matx(−β). */
static void compute_T_cross(int i, double beta, M4 T) {
    M4 MX, A, B;
    mat4_rotx(-beta, MX);
    mat4_mul(MATZ_PI, MX, A);          /* MATZ_PI · matx(−β) */
    mat4_mul(BOOST_X_S, A, B);         /* boost · matz(π) · matx(−β) */
    mat4_mul(R_POW[i], B, T);          /* R_pow[i] · ... */
}

/* ---------- BFS placement -------------------------------------------------- */
static int bfs_place(void) {
    static int q[MAXF];
    memset(FACE_PLACED, 0, NF*sizeof(int));
    memset(V_PLACED, 0, (NV+2)*sizeof(int));

    /* Base face: F_M = identity, labels follow stored CCW order. */
    mat4_iden(FRAMES[0]);
    LABELS[0][0] = FACES[0].a;
    LABELS[0][1] = FACES[0].b;
    LABELS[0][2] = FACES[0].c;
    for (int i = 0; i < 3; i++) {
        memcpy(V[LABELS[0][i]], CANON[i], 4*sizeof(double));
        V_PLACED[LABELS[0][i]] = 1;
    }
    FACE_PLACED[0] = 1;

    int qh = 0, qt = 0;
    q[qt++] = 0;

    while (qh < qt) {
        int fi = q[qh++];
        for (int i = 0; i < 3; i++) {
            int a = LABELS[fi][i];
            int b = LABELS[fi][(i+1)%3];
            int other_fi = EM_F[b][a];          /* face holding reverse edge */
            if (other_fi < 0) {
                fprintf(stderr, "ERROR: edge (%d,%d) has no opposite face\n", a, b);
                return -1;
            }
            if (FACE_PLACED[other_fi]) continue;

            int oa = FACES[other_fi].a;
            int ob = FACES[other_fi].b;
            int oc = FACES[other_fi].c;
            int c_id = (oa!=a && oa!=b) ? oa : ((ob!=a && ob!=b) ? ob : oc);

            int e = EDGE_IDX[a][b];
            if (e < 0) {
                fprintf(stderr, "ERROR: missing edge index for (%d,%d)\n", a, b);
                return -1;
            }
            double beta = BENDS[e];

            M4 T_cross_i, F_N;
            compute_T_cross(i, beta, T_cross_i);
            mat4_mul(FRAMES[fi], T_cross_i, F_N);

            memcpy(FRAMES[other_fi], F_N, sizeof(M4));
            LABELS[other_fi][0] = b;            /* N's V0 */
            LABELS[other_fi][1] = a;            /* N's V1 */
            LABELS[other_fi][2] = c_id;         /* N's V2 */
            FACE_PLACED[other_fi] = 1;

            if (!V_PLACED[c_id]) {
                mat4_apply(F_N, CANON[2], V[c_id]);
                V_PLACED[c_id] = 1;
            }
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
    fprintf(fh, "# Klein realization (realize_c v2) — α=%.10g°, side s=%.10g\n",
            ALPHA_DEG, S_LEN);
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
