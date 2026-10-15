/* pencils.c -- enumerate 4-pencil line arrangements with collinear centres on H = {z=0}
 * pencils through (1:0:0), (0:1:0), (1:1:0), (p:q:0); lines q_i x - p_i y = v z.
 * Pencils A, B, C have integer offsets in a box; pencil D (direction (p:q)) is derived:
 * its lines pass through >= 3 points of the A,B,C configuration (needed for <= 8 points per line).
 * Output: configurations with b2 = 80 (17 lines with H), every line with 4..8 points.
 * usage: pencils sA sB sC sD K p q
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
static ll gcdll(ll a, ll b){ if(a<0)a=-a; if(b<0)b=-b; while(b){ ll t=a%b; a=b; b=t;} return a; }
static int sA,sB,sC,sD,K; static ll P,Q;
static int dirp[4]={1,0,1,0}, dirq[4]={0,1,1,0};
typedef struct { ll v; int d; } Line;          /* line: q_d x - p_d y = v z   (v scaled by den) */
typedef struct { ll x,y,z; } Pt;
static int cmpPt(const void*a,const void*b){ const Pt*p=a,*q=b; if(p->z!=q->z) return p->z<q->z?-1:1; if(p->x!=q->x) return p->x<q->x?-1:1; if(p->y!=q->y) return p->y<q->y?-1:1; return 0; }
/* intersection of lines (i,v) and (j,w): coefficient vectors (q_i,-p_i,-v),(q_j,-p_j,-w) */
static Pt meet(int i, ll v, int j, ll w){
  ll X = (ll)dirp[i]*w - v*(ll)dirp[j];
  ll Y = -v*(ll)dirq[j] + (ll)dirq[i]*w;
  ll Z = -(ll)dirq[i]*dirp[j] + (ll)dirp[i]*dirq[j];
  if(Z<0){X=-X;Y=-Y;Z=-Z;}
  ll g=gcdll(gcdll(X,Y),Z); if(g>1){X/=g;Y/=g;Z/=g;}
  Pt p={X,Y,Z}; return p;
}
typedef struct { ll num, den; } Frac;
static Frac mkfrac(ll n, ll d){ if(d<0){n=-n;d=-d;} ll g=gcdll(n,d); if(g>1){n/=g;d/=g;} Frac f={n,d}; return f; }
static int cmpFrac(const void*a,const void*b){ const Frac*p=a,*q=b; if(p->den!=q->den) return p->den<q->den?-1:1; if(p->num!=q->num) return p->num<q->num?-1:1; return 0; }
static ll A[8],B[8],C[8]; static Frac D[8];
static long nfound=0, ncheck=0;
/* full lattice check: lines of A,B,C (integer v, den 1) and D (v = num/den). To keep integers, scale
   the whole picture by L = lcm of D denominators: x,y -> L x, L y; then A,B,C offsets become L*v and D offsets integers. */
static int checkAll(ll L, int quiet){
  Line ln[20]; int n=0;
  for(int i=0;i<sA;i++){ln[n].d=0; ln[n].v=A[i]*L; n++;}
  for(int i=0;i<sB;i++){ln[n].d=1; ln[n].v=B[i]*L; n++;}
  for(int i=0;i<sC;i++){ln[n].d=2; ln[n].v=C[i]*L; n++;}
  for(int i=0;i<sD;i++){ln[n].d=3; ln[n].v=D[i].num*(L/D[i].den); n++;}
  /* points */
  Pt pts[200]; int lineOf[200][2]; int np=0;
  for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if(ln[i].d==ln[j].d) continue; pts[np]=meet(ln[i].d,ln[i].v,ln[j].d,ln[j].v); lineOf[np][0]=i; lineOf[np][1]=j; np++; }
  /* group: sort indices by point */
  int idx[200]; for(int i=0;i<np;i++) idx[i]=i;
  /* simple insertion sort on idx by pts */
  for(int i=1;i<np;i++){ int k=idx[i]; int j=i-1; while(j>=0 && cmpPt(&pts[idx[j]],&pts[k])>0){ idx[j+1]=idx[j]; j--; } idx[j+1]=k; }
  int b2=0; int cnt[20]; for(int i=0;i<n;i++) cnt[i]=1; /* the centre on H */
  int mult[20]; int maxm=0;
  int i=0; int onLine[20];
  while(i<np){
    int j=i; while(j<np && cmpPt(&pts[idx[i]],&pts[idx[j]])==0) j++;
    int pairs=j-i; /* pairs = C(m,2) */
    int m=2; while(m*(m-1)/2<pairs) m++;
    if(m*(m-1)/2!=pairs) return 0; /* inconsistent (should not happen) */
    b2 += m-1; if(m>maxm) maxm=m;
    /* lines through this point */
    memset(onLine,0,sizeof(onLine));
    for(int t=i;t<j;t++){ onLine[lineOf[idx[t]][0]]=1; onLine[lineOf[idx[t]][1]]=1; }
    for(int t=0;t<n;t++) if(onLine[t]) cnt[t]++;
    i=j;
  }
  if(b2!=64) return 0;
  if(maxm>5) return 0;
  { int hi = getenv("MAXPTS") ? atoi(getenv("MAXPTS")) : 8; for(int t=0;t<n;t++) if(cnt[t]<4||cnt[t]>hi) return 0; }
  if(quiet) return 1;
  nfound++;
  printf("FOUND sizes %d %d %d %d lambda %lld/%lld A", sA,sB,sC,sD,P,Q);
  for(int t=0;t<sA;t++) printf(" %lld",A[t]); printf(" B"); for(int t=0;t<sB;t++) printf(" %lld",B[t]);
  printf(" C"); for(int t=0;t<sC;t++) printf(" %lld",C[t]); printf(" D"); for(int t=0;t<sD;t++) printf(" %lld/%lld",D[t].num,D[t].den);
  printf(" perline"); for(int t=0;t<n;t++) printf(" %d",cnt[t]); printf("\n"); fflush(stdout);
  return 1;
}
static ll lcmll(ll a, ll b){ return a/gcdll(a,b)*b; }
static void tryD(){
  /* candidate D offsets: v = (q X - p Y)/Z over all A,B,C pairwise points; keep those hit by >= 3 points */
  Frac vals[100]; int nv=0;
  for(int i=0;i<sA;i++) for(int j=0;j<sB;j++){ Pt p=meet(0,A[i],1,B[j]); vals[nv++]=mkfrac(Q*p.x-P*p.y,p.z);}
  for(int i=0;i<sA;i++) for(int j=0;j<sC;j++){ Pt p=meet(0,A[i],2,C[j]); vals[nv++]=mkfrac(Q*p.x-P*p.y,p.z);}
  for(int i=0;i<sB;i++) for(int j=0;j<sC;j++){ Pt p=meet(1,B[i],2,C[j]); vals[nv++]=mkfrac(Q*p.x-P*p.y,p.z);}
  qsort(vals,nv,sizeof(Frac),cmpFrac);
  Frac cand[100]; int nc=0;
  for(int i=0;i<nv;){ int j=i; while(j<nv && cmpFrac(&vals[i],&vals[j])==0) j++; if(j-i>=3) cand[nc++]=vals[i]; i=j; }
  if(nc<sD) return;
  /* choose sD of the candidates */
  int sel[8]; for(int i=0;i<sD;i++) sel[i]=i;
  while(1){
    ll L=1; for(int i=0;i<sD;i++){ D[i]=cand[sel[i]]; L=lcmll(L,D[i].den); }
    ncheck++; checkAll(L,0);
    int i=sD-1; while(i>=0 && sel[i]==nc-sD+i) i--; if(i<0) break; sel[i]++; for(int j=i+1;j<sD;j++) sel[j]=sel[j-1]+1;
  }
}
static void enumC(){
  /* C subset of [-K,K] of size sC, increasing */
  int range=2*K+1; int sel[8]; for(int i=0;i<sC;i++) sel[i]=i;
  while(1){
    for(int i=0;i<sC;i++) C[i]=sel[i]-K;
    tryD();
    int i=sC-1; while(i>=0 && sel[i]==range-sC+i) i--; if(i<0) break; sel[i]++; for(int j=i+1;j<sC;j++) sel[j]=sel[j-1]+1;
  }
}
int main(int argc,char**argv){
  if(argc<8){fprintf(stderr,"usage: pencils sA sB sC sD K p q\n"); return 1;}
  sA=atoi(argv[1]); sB=atoi(argv[2]); sC=atoi(argv[3]); sD=atoi(argv[4]); K=atoi(argv[5]); P=atoll(argv[6]); Q=atoll(argv[7]);
  dirp[3]=(int)P; dirq[3]=(int)Q;
  /* A = {0} u (sA-1)-subset of [1,K]; B = {0} u (sB-1)-subset of [1,K] */
  int selA[8], selB[8];
  for(int i=0;i<sA-1;i++) selA[i]=i+1;
  while(1){
    A[0]=0; for(int i=0;i<sA-1;i++) A[i+1]=selA[i];
    for(int i=0;i<sB-1;i++) selB[i]=i+1;
    while(1){
      B[0]=0; for(int i=0;i<sB-1;i++) B[i+1]=selB[i];
      enumC();
      int i=sB-2; while(i>=0 && selB[i]==K-(sB-1)+i+1) i--; if(i<0) break; selB[i]++; for(int j=i+1;j<sB-1;j++) selB[j]=selB[j-1]+1;
    }
    int i=sA-2; while(i>=0 && selA[i]==K-(sA-1)+i+1) i--; if(i<0) break; selA[i]++; for(int j=i+1;j<sA-1;j++) selA[j]=selA[j-1]+1;
  }
  fprintf(stderr,"sizes %d %d %d %d K=%d lambda=%lld/%lld: checked %ld D-choices, found %ld\n",sA,sB,sC,sD,K,P,Q,ncheck,nfound);
  return 0;
}
