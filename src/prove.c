/* prove.c — prove existence of unit-edge-length polyhedra
 *
 * C port of prove_float.py. Uses IEEE 754 double precision with
 * rigorous error tracking. Requires LAPACK for SVD (dgesvd).
 *
 * Checks:
 *   0. Undented (all vertex turning sums > 0)
 *   0b. Embedded (no triangle-triangle intersections)
 *   1. 3|V| >= |E|
 *   2. Collision distance CD > 0
 *   3. sigma_min > 0 and rho < sigma_min^2 / (16*sqrt(E))
 *   4. Perturbation bound < CD/sqrt(V)
 *
 * Error bounds (IEEE 754):
 *   sigma_min: |sigma_true - sigma_float| <= 5*n*eps*||J||_F  (Wedin)
 *   rho:       |rho_true - rho_float| <= (2E+10)*eps*rho      (accumulation)
 *   CD:        |CD_true - CD_float| <= 100*eps*CD              (arithmetic)
 *
 * Usage:   prove objdir/ | prove file.obj
 *          Reads OBJ file(s), prints PASS/FAIL per net.
 * Compile: cc -O3 -o prove prove.c -lm -llapack -lblas
 *          (do NOT use -ffast-math)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dirent.h>

#define MAXV   400
#define MAXF   (2*MAXV)
#define MAXE   (3*MAXV)
#define MAXRING 12

#define EPS 2.220446049250313e-16  /* DBL_EPSILON */

/* LAPACK SVD */
extern void dgesvd_(char *jobu, char *jobvt, int *m, int *n, double *a,
                    int *lda, double *s, double *u, int *ldu, double *vt,
                    int *ldvt, double *work, int *lwork, int *info);

/* ── geometry ────────────────────────────────────────────────────────────── */
static double vx[MAXV], vy[MAXV], vz[MAXV];
static int fa[MAXF], fb[MAXF], fc[MAXF];
static int eu[MAXE], ev[MAXE];
static int nv, nf, ne;
static int em[MAXV][MAXV];  /* edge map: em[a][b] = opposite vertex */
static int nbr[MAXV][MAXRING], nnbr[MAXV];

static void nbr_add(int u, int w){
    for(int i=0;i<nnbr[u];i++) if(nbr[u][i]==w) return;
    nbr[u][nnbr[u]++]=w;
}

static int read_obj(const char *path){
    FILE *fp=fopen(path,"r"); if(!fp) return 0;
    char line[1024]; nv=nf=ne=0;
    memset(em,0,sizeof(em)); memset(nnbr,0,sizeof(nnbr));
    while(fgets(line,sizeof(line),fp)){
        if(line[0]=='v'&&line[1]==' ')
            sscanf(line+2,"%lf%lf%lf",&vx[nv],&vy[nv],&vz[nv]), nv++;
        else if(line[0]=='f'){
            int a,b,c; sscanf(line+2,"%d%d%d",&a,&b,&c); a--;b--;c--;
            fa[nf]=a;fb[nf]=b;fc[nf]=c; nf++;
            em[a][b]=c+1; em[b][c]=a+1; em[c][a]=b+1; /* +1 so 0 = no edge */
            nbr_add(a,b);nbr_add(b,a);nbr_add(b,c);nbr_add(c,b);nbr_add(a,c);nbr_add(c,a);
        }
    }
    fclose(fp);
    /* collect edges */
    ne=0;
    for(int i=0;i<nv;i++) for(int j=i+1;j<nv;j++)
        if(em[i][j]){eu[ne]=i;ev[ne]=j;ne++;}
    return nv;
}

static int cyclic_ring(int v, int ring[]){
    if(!nnbr[v]) return 0;
    int start=nbr[v][0]; ring[0]=start; int k=1,cur=start;
    for(int s=0;s<20;s++){
        int nxt=em[v][cur]-1; /* -1 because stored +1 */
        if(nxt<0||nxt==start) break;
        ring[k++]=nxt; cur=nxt;
    }
    return k;
}

/* ── Check 0: undented ───────────────────────────────────────────────────── */
static double signed_sph_angle(double ax,double ay,double az,
                                double bx,double by,double bz,
                                double cx,double cy,double cz){
    double acx=ay*cz-az*cy,acy=az*cx-ax*cz,acz=ax*cy-ay*cx;
    double num=bx*acx+by*acy+bz*acz;
    double den=(ax*bx+ay*by+az*bz)*(bx*cx+by*cy+bz*cz)-(ax*cx+ay*cy+az*cz);
    return atan2(num,den);
}

static double compute_undented(void){
    /* Returns min turning. Undented if min_turning > 100*V*eps */
    double min_t=1e30;
    for(int v=0;v<nv;v++){
        int ring[MAXRING];
        int k=cyclic_ring(v,ring); if(k<3) continue;
        double dx[MAXRING],dy[MAXRING],dz[MAXRING];
        for(int i=0;i<k;i++){
            int w=ring[i];
            double ex=vx[w]-vx[v],ey=vy[w]-vy[v],ez=vz[w]-vz[v];
            double L=sqrt(ex*ex+ey*ey+ez*ez); if(L<1e-30)L=1e-30;
            dx[i]=ex/L;dy[i]=ey/L;dz[i]=ez/L;
        }
        double tv=0;
        for(int i=0;i<k;i++)
            tv+=signed_sph_angle(dx[(i-1+k)%k],dy[(i-1+k)%k],dz[(i-1+k)%k],
                                  dx[i],dy[i],dz[i],
                                  dx[(i+1)%k],dy[(i+1)%k],dz[(i+1)%k]);
        if(tv<min_t) min_t=tv;
    }
    return min_t;
}

/* ── Check 0b: embedded (triangle-triangle intersection) ─────────────── */
static int tri_tri_intersect(int fi, int fj){
    double p0x=vx[fa[fi]],p0y=vy[fa[fi]],p0z=vz[fa[fi]];
    double p1x=vx[fb[fi]],p1y=vy[fb[fi]],p1z=vz[fb[fi]];
    double p2x=vx[fc[fi]],p2y=vy[fc[fi]],p2z=vz[fc[fi]];
    double q0x=vx[fa[fj]],q0y=vy[fa[fj]],q0z=vz[fa[fj]];
    double q1x=vx[fb[fj]],q1y=vy[fb[fj]],q1z=vz[fb[fj]];
    double q2x=vx[fc[fj]],q2y=vy[fc[fj]],q2z=vz[fc[fj]];
    /* plane of tri 1 */
    double e1x=p1x-p0x,e1y=p1y-p0y,e1z=p1z-p0z;
    double e2x=p2x-p0x,e2y=p2y-p0y,e2z=p2z-p0z;
    double n1x=e1y*e2z-e1z*e2y,n1y=e1z*e2x-e1x*e2z,n1z=e1x*e2y-e1y*e2x;
    double d1=-(n1x*p0x+n1y*p0y+n1z*p0z);
    double dq0=n1x*q0x+n1y*q0y+n1z*q0z+d1;
    double dq1=n1x*q1x+n1y*q1y+n1z*q1z+d1;
    double dq2=n1x*q2x+n1y*q2y+n1z*q2z+d1;
    if(dq0*dq1>0&&dq0*dq2>0) return 0;
    /* plane of tri 2 */
    double f1x=q1x-q0x,f1y=q1y-q0y,f1z=q1z-q0z;
    double f2x=q2x-q0x,f2y=q2y-q0y,f2z=q2z-q0z;
    double n2x=f1y*f2z-f1z*f2y,n2y=f1z*f2x-f1x*f2z,n2z=f1x*f2y-f1y*f2x;
    double d2=-(n2x*q0x+n2y*q0y+n2z*q0z);
    double dp0=n2x*p0x+n2y*p0y+n2z*p0z+d2;
    double dp1=n2x*p1x+n2y*p1y+n2z*p1z+d2;
    double dp2=n2x*p2x+n2y*p2y+n2z*p2z+d2;
    if(dp0*dp1>0&&dp0*dp2>0) return 0;
    /* coplanar check */
    double nn1=sqrt(n1x*n1x+n1y*n1y+n1z*n1z);
    double nn2=sqrt(n2x*n2x+n2y*n2y+n2z*n2z);
    if(nn1<1e-30||nn2<1e-30) return 0;
    double Dx=n1y*n2z-n1z*n2y,Dy=n1z*n2x-n1x*n2z,Dz=n1x*n2y-n1y*n2x;
    double Dn=sqrt(Dx*Dx+Dy*Dy+Dz*Dz);
    if(Dn<1e-30) return 0; /* coplanar */
    /* project */
    int ax=0; double am=fabs(Dx);
    if(fabs(Dy)>am){ax=1;am=fabs(Dy);} if(fabs(Dz)>am) ax=2;
    double pp[3],qq[3],dp[3],dq[3];
    double pv[3][3]={{p0x,p0y,p0z},{p1x,p1y,p1z},{p2x,p2y,p2z}};
    double qv[3][3]={{q0x,q0y,q0z},{q1x,q1y,q1z},{q2x,q2y,q2z}};
    for(int i=0;i<3;i++){pp[i]=pv[i][ax];qq[i]=qv[i][ax];}
    dp[0]=dp0;dp[1]=dp1;dp[2]=dp2;
    dq[0]=dq0;dq[1]=dq1;dq[2]=dq2;
    double eps=1e-12*(nn1>nn2?nn1:nn2);
    if(fabs(dp[0])<eps&&fabs(dp[1])<eps&&fabs(dp[2])<eps) return 0;
    if(fabs(dq[0])<eps&&fabs(dq[1])<eps&&fabs(dq[2])<eps) return 0;
    /* sign check */
    double sp[3],sq[3];
    for(int i=0;i<3;i++){sp[i]=(fabs(dp[i])<eps)?0:(dp[i]>0?1:-1);sq[i]=(fabs(dq[i])<eps)?0:(dq[i]>0?1:-1);}
    if((sp[0]>=0&&sp[1]>=0&&sp[2]>=0)||(sp[0]<=0&&sp[1]<=0&&sp[2]<=0)) return 0;
    if((sq[0]>=0&&sq[1]>=0&&sq[2]>=0)||(sq[0]<=0&&sq[1]<=0&&sq[2]<=0)) return 0;
    /* intervals */
    int lone_p,pa0,pa1, lone_q,qa0,qa1;
    if(sp[0]*sp[1]>0||(sp[0]==0&&sp[1]==0)){lone_p=2;pa0=0;pa1=1;}
    else if(sp[0]*sp[2]>0||(sp[0]==0&&sp[2]==0)){lone_p=1;pa0=0;pa1=2;}
    else{lone_p=0;pa0=1;pa1=2;}
    if(sq[0]*sq[1]>0||(sq[0]==0&&sq[1]==0)){lone_q=2;qa0=0;qa1=1;}
    else if(sq[0]*sq[2]>0||(sq[0]==0&&sq[2]==0)){lone_q=1;qa0=0;qa1=2;}
    else{lone_q=0;qa0=1;qa1=2;}
    double t0p=pp[pa0]+(pp[lone_p]-pp[pa0])*dp[pa0]/(dp[pa0]-dp[lone_p]);
    double t1p=pp[pa1]+(pp[lone_p]-pp[pa1])*dp[pa1]/(dp[pa1]-dp[lone_p]);
    if(t0p>t1p){double t=t0p;t0p=t1p;t1p=t;}
    double t0q=qq[qa0]+(qq[lone_q]-qq[qa0])*dq[qa0]/(dq[qa0]-dq[lone_q]);
    double t1q=qq[qa1]+(qq[lone_q]-qq[qa1])*dq[qa1]/(dq[qa1]-dq[lone_q]);
    if(t0q>t1q){double t=t0q;t0q=t1q;t1q=t;}
    double overlap=(t1p<t1q?t1p:t1q)-(t0p>t0q?t0p:t0q);
    return overlap>eps;
}

static int check_embedded(int *n_bad_out){
    int n_bad=0;
    for(int i=0;i<nf;i++) for(int j=i+1;j<nf;j++){
        int a=fa[i],b=fb[i],c=fc[i],d=fa[j],e=fb[j],f=fc[j];
        int shared=0;
        if(a==d||a==e||a==f) shared++;
        if(b==d||b==e||b==f) shared++;
        if(c==d||c==e||c==f) shared++;
        if(shared>=2) continue;
        if(tri_tri_intersect(i,j)) n_bad++;
    }
    *n_bad_out=n_bad;
    return n_bad==0;
}

/* ── Check 2: collision distance ─────────────────────────────────────── */
static double pt_pt_sq(int i, int j){
    double dx=vx[i]-vx[j],dy=vy[i]-vy[j],dz=vz[i]-vz[j];
    return dx*dx+dy*dy+dz*dz;
}
static double pt_seg_sq(int p, int a, int b){
    double abx=vx[b]-vx[a],aby=vy[b]-vy[a],abz=vz[b]-vz[a];
    double apx=vx[p]-vx[a],apy=vy[p]-vy[a],apz=vz[p]-vz[a];
    double ab2=abx*abx+aby*aby+abz*abz;
    double t=(apx*abx+apy*aby+apz*abz)/(ab2>1e-300?ab2:1e-300);
    if(t<0)t=0; if(t>1)t=1;
    double cx=vx[a]+t*abx-vx[p],cy=vy[a]+t*aby-vy[p],cz=vz[a]+t*abz-vz[p];
    return cx*cx+cy*cy+cz*cz;
}
static double seg_seg_sq(int a0,int a1,int b0,int b1){
    double d1x=vx[a1]-vx[a0],d1y=vy[a1]-vy[a0],d1z=vz[a1]-vz[a0];
    double d2x=vx[b1]-vx[b0],d2y=vy[b1]-vy[b0],d2z=vz[b1]-vz[b0];
    double rx=vx[a0]-vx[b0],ry=vy[a0]-vy[b0],rz=vz[a0]-vz[b0];
    double a=d1x*d1x+d1y*d1y+d1z*d1z;
    double e=d2x*d2x+d2y*d2y+d2z*d2z;
    double f=d2x*rx+d2y*ry+d2z*rz;
    if(a<1e-300&&e<1e-300) return rx*rx+ry*ry+rz*rz;
    double c=d1x*rx+d1y*ry+d1z*rz;
    double bv=d1x*d2x+d1y*d2y+d1z*d2z;
    if(a<1e-300) return pt_seg_sq(a0,b0,b1);
    if(e<1e-300) return pt_seg_sq(b0,a0,a1);
    double det=a*e-bv*bv;
    double s=(fabs(det)>1e-300)?(bv*f-c*e)/det:0;
    if(s<0)s=0;if(s>1)s=1;
    double t=(bv*s+f)/e;
    if(t<0)t=0;if(t>1)t=1;
    s=(bv*t-c)/a;
    if(s<0)s=0;if(s>1)s=1;
    double dx=vx[a0]+s*d1x-vx[b0]-t*d2x;
    double dy=vy[a0]+s*d1y-vy[b0]-t*d2y;
    double dz=vz[a0]+s*d1z-vz[b0]-t*d2z;
    return dx*dx+dy*dy+dz*dz;
}
static double pt_tri_sq(int p, int ta, int tb, int tc){
    double abx=vx[tb]-vx[ta],aby=vy[tb]-vy[ta],abz=vz[tb]-vz[ta];
    double acx=vx[tc]-vx[ta],acy=vy[tc]-vy[ta],acz=vz[tc]-vz[ta];
    double apx=vx[p]-vx[ta],apy=vy[p]-vy[ta],apz=vz[p]-vz[ta];
    double d00=abx*abx+aby*aby+abz*abz;
    double d01=abx*acx+aby*acy+abz*acz;
    double d11=acx*acx+acy*acy+acz*acz;
    double d20=apx*abx+apy*aby+apz*abz;
    double d21=apx*acx+apy*acy+apz*acz;
    double det=d00*d11-d01*d01;
    if(fabs(det)<1e-300){
        double r1=pt_seg_sq(p,ta,tb),r2=pt_seg_sq(p,tb,tc),r3=pt_seg_sq(p,ta,tc);
        return r1<r2?(r1<r3?r1:r3):(r2<r3?r2:r3);
    }
    double v=(d11*d20-d01*d21)/det;
    double w=(d00*d21-d01*d20)/det;
    double u=1-v-w;
    if(u>=0&&v>=0&&w>=0){
        double px2=u*vx[ta]+v*vx[tb]+w*vx[tc]-vx[p];
        double py2=u*vy[ta]+v*vy[tb]+w*vy[tc]-vy[p];
        double pz2=u*vz[ta]+v*vz[tb]+w*vz[tc]-vz[p];
        return px2*px2+py2*py2+pz2*pz2;
    }
    double r1=pt_seg_sq(p,ta,tb),r2=pt_seg_sq(p,tb,tc),r3=pt_seg_sq(p,ta,tc);
    return r1<r2?(r1<r3?r1:r3):(r2<r3?r2:r3);
}

static double compute_cd(void){
    /* adjacency: vertices sharing a face */
    static int adj[MAXV][MAXV];
    memset(adj,0,sizeof(adj));
    for(int f=0;f<nf;f++){
        int a=fa[f],b=fb[f],c=fc[f];
        adj[a][b]=adj[b][a]=1;
        adj[a][c]=adj[c][a]=1;
        adj[b][c]=adj[c][b]=1;
    }
    double min_sq=1e30;
    /* vertex-vertex */
    for(int i=0;i<nv;i++) for(int j=i+1;j<nv;j++){
        if(adj[i][j]) continue;
        double d=pt_pt_sq(i,j); if(d<min_sq) min_sq=d;
    }
    /* vertex-edge */
    for(int v=0;v<nv;v++) for(int e=0;e<ne;e++){
        if(v==eu[e]||v==ev[e]||adj[v][eu[e]]||adj[v][ev[e]]) continue;
        double d=pt_seg_sq(v,eu[e],ev[e]); if(d<min_sq) min_sq=d;
    }
    /* vertex-face */
    for(int v=0;v<nv;v++) for(int f=0;f<nf;f++){
        if(v==fa[f]||v==fb[f]||v==fc[f]) continue;
        if(adj[v][fa[f]]||adj[v][fb[f]]||adj[v][fc[f]]) continue;
        double d=pt_tri_sq(v,fa[f],fb[f],fc[f]); if(d<min_sq) min_sq=d;
    }
    /* edge-edge */
    for(int i=0;i<ne;i++) for(int j=i+1;j<ne;j++){
        if(eu[i]==eu[j]||eu[i]==ev[j]||ev[i]==eu[j]||ev[i]==ev[j]) continue;
        if(adj[eu[i]][eu[j]]||adj[eu[i]][ev[j]]||adj[ev[i]][eu[j]]||adj[ev[i]][ev[j]]) continue;
        double d=seg_seg_sq(eu[i],ev[i],eu[j],ev[j]); if(d<min_sq) min_sq=d;
    }
    return sqrt(min_sq);
}

/* ── Check 3: sigma_min via SVD, rho via residuals ───────────────────── */
static double compute_sigma_min(double *sigma_lower_out){
    /* Build Jacobian E x 3V (column-major for LAPACK) */
    int m=ne, n=3*nv;
    static double J[MAXE*3*MAXV]; /* column-major */
    memset(J,0,sizeof(double)*m*n);
    for(int k=0;k<ne;k++){
        int i=eu[k],j=ev[k];
        double dx=vx[i]-vx[j],dy=vy[i]-vy[j],dz=vz[i]-vz[j];
        J[k+m*(3*i  )]=2*dx; J[k+m*(3*i+1)]=2*dy; J[k+m*(3*i+2)]=2*dz;
        J[k+m*(3*j  )]=-2*dx;J[k+m*(3*j+1)]=-2*dy;J[k+m*(3*j+2)]=-2*dz;
    }
    /* SVD */
    int mn=m<n?m:n;
    static double S[MAXE];
    double work_query; int lwork=-1, info;
    char jobu='N', jobvt='N';
    dgesvd_(&jobu,&jobvt,&m,&n,J,&m,S,NULL,&m,NULL,&n,&work_query,&lwork,&info);
    lwork=(int)work_query;
    static double work[1000000];
    if(lwork>1000000) lwork=1000000;
    dgesvd_(&jobu,&jobvt,&m,&n,J,&m,S,NULL,&m,NULL,&n,work,&lwork,&info);
    if(info!=0) return -1;
    double sigma_min=S[mn-1];
    /* Frobenius norm (recompute since dgesvd destroys J) */
    double norm_f=0;
    for(int i=0;i<mn;i++) norm_f+=S[i]*S[i];
    norm_f=sqrt(norm_f);
    /* Wedin bound */
    int nn=m>n?m:n;
    double sv_err=5.0*nn*EPS*norm_f;
    *sigma_lower_out=sigma_min-sv_err;
    return sigma_min;
}

static double compute_rho(double *rho_upper_out){
    double rho_sq=0;
    for(int k=0;k<ne;k++){
        int i=eu[k],j=ev[k];
        double dx=vx[i]-vx[j],dy=vy[i]-vy[j],dz=vz[i]-vz[j];
        double r=dx*dx+dy*dy+dz*dz-1.0;
        rho_sq+=r*r;
    }
    double rho=sqrt(rho_sq);
    double rho_rel_err=(2*ne+10)*EPS;
    *rho_upper_out=rho*(1+rho_rel_err)+ne*5*EPS;
    return rho;
}

/* ── prove ───────────────────────────────────────────────────────────────── */
static int prove(const char *name, char *msg, int msglen){
    /* Check 0: undented */
    double min_turn=compute_undented();
    double turn_bound=100.0*nv*EPS;
    if(min_turn<=turn_bound){
        snprintf(msg,msglen,"Failed: dented (min_turn=%.6e)",min_turn);
        return 0;
    }
    /* Check 0b: embedded */
    int n_bad;
    if(!check_embedded(&n_bad)){
        snprintf(msg,msglen,"Failed: self-intersecting (%d pairs)",n_bad);
        return 0;
    }
    /* Inequality 1: 3V >= E */
    if(3*nv<ne){
        snprintf(msg,msglen,"Failed: 3V=%d < E=%d",3*nv,ne);
        return 0;
    }
    /* Inequality 2: CD > 0 */
    double cd=compute_cd();
    double cd_lower=cd*(1-100*EPS);
    if(cd_lower<=0){
        snprintf(msg,msglen,"Failed: CD not provably positive (%.6e)",cd);
        return 0;
    }
    /* Inequality 3: sigma_min > 0 and rho < sigma_min^2 / (16*sqrt(E)) */
    double sigma_lower;
    double sigma=compute_sigma_min(&sigma_lower);
    if(sigma<0){
        snprintf(msg,msglen,"Failed: SVD error");
        return 0;
    }
    if(sigma_lower<=0){
        snprintf(msg,msglen,"Failed: sigma_min not provably positive (%.6e, bound=%.6e)",sigma,sigma_lower);
        return 0;
    }
    double rho_upper;
    double rho=compute_rho(&rho_upper);
    double sqrtE=sqrt((double)ne);
    double threshold=sigma_lower*sigma_lower/(16.0*sqrtE);
    if(rho_upper>=threshold){
        snprintf(msg,msglen,"Failed: rho=%.6e >= threshold=%.6e",rho_upper,threshold);
        return 0;
    }
    /* Inequality 4: perturbation bound < CD/sqrt(V) */
    double disc=sigma_lower*sigma_lower-16*rho_upper*sqrtE;
    if(disc<=0){
        snprintf(msg,msglen,"Failed: discriminant non-positive");
        return 0;
    }
    double lhs=(sigma_lower-sqrt(disc))/(8*sqrtE);
    double rhs=cd_lower/sqrt((double)nv);
    if(lhs>=rhs){
        snprintf(msg,msglen,"Failed: LHS=%.6e >= CD/sqrt(V)=%.6e",lhs,rhs);
        return 0;
    }
    snprintf(msg,msglen,"Success: existence proven");
    return 1;
}

/* ── main ────────────────────────────────────────────────────────────────── */
int main(int argc, char **argv){
    if(argc<2){fprintf(stderr,"usage: prove file.obj | prove objdir/\n");return 1;}

    /* collect paths */
    static char paths[2000000][256];
    int npaths=0;

    struct dirent *de;
    DIR *dr=opendir(argv[1]);
    if(dr){
        while((de=readdir(dr))!=NULL){
            int len=strlen(de->d_name);
            if(len>4&&!strcmp(de->d_name+len-4,".obj")){
                snprintf(paths[npaths],256,"%s/%s",argv[1],de->d_name);
                npaths++;
            }
        }
        closedir(dr);
    } else {
        strncpy(paths[0],argv[1],255);
        npaths=1;
    }

    int n_pass=0,n_fail=0;
    for(int p=0;p<npaths;p++){
        if(!read_obj(paths[p])){
            fprintf(stderr,"cannot read: %s\n",paths[p]);
            continue;
        }
        /* extract name */
        char *base=strrchr(paths[p],'/');
        base=base?base+1:paths[p];
        char name[256]; strncpy(name,base,255);
        char *dot=strrchr(name,'.'); if(dot) *dot='\0';

        char msg[256];
        int ok=prove(name,msg,sizeof(msg));
        printf("  %-40s %3d %3d %3d  %s  %s\n",name,nv,ne,nf,ok?"PASS":"FAIL",msg);
        if(ok) n_pass++; else n_fail++;
    }
    printf("# pass=%d  fail=%d\n",n_pass,n_fail);
    return 0;
}
