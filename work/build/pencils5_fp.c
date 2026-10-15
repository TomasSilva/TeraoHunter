/* pencils5_fp.c -- census over F_p of 5-pencil arrangements with collinear centres
 * (1:0:0),(0:1:0),(1:1:0),(l1:1:0),(l2:1:0); sizes sA sB sC sE sD; pencils A,B,C,E boxed (all of F_p),
 * pencil D (direction (l2:1)) derived: its lines pass through >= need points of the A,B,C,E configuration.
 * Normalisation a0=0,a1=1 (sA>=2), b0=0.  Output: b2 = 80 with H, every line 4..8 points, mult <= 6.
 * usage: pencils5_fp sA sB sC sE sD p l1 l2
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int sz[5],p,l1,l2; static int dirp[5]={1,0,1,0,0}, dirq[5]={0,1,1,1,1};
static int md(long long a){ a%=p; if(a<0)a+=p; return (int)a; }
static int inv(int a){ int r=1,b=a%p,e=p-2; if(b<0)b+=p; while(e){ if(e&1) r=(int)((long long)r*b%p); b=(int)((long long)b*b%p); e>>=1;} return r; }
static int meetKey(int i,int v,int j,int w){
  int X=md((long long)dirp[i]*w - (long long)v*dirp[j]); int Y=md(-(long long)v*dirq[j] + (long long)dirq[i]*w); int Z=md(-(long long)dirq[i]*dirp[j] + (long long)dirp[i]*dirq[j]); int s;
  if(X){ s=inv(X); X=1; Y=md((long long)Y*s); Z=md((long long)Z*s);} else if(Y){ s=inv(Y); Y=1; Z=md((long long)Z*s);} else Z=1;
  return (X*p+Y)*p+Z; }
static void meetXY(int i,int v,int j,int w,int*xo,int*yo){
  int X=md((long long)dirp[i]*w - (long long)v*dirp[j]); int Y=md(-(long long)v*dirq[j] + (long long)dirq[i]*w); int Z=md(-(long long)dirq[i]*dirp[j] + (long long)dirp[i]*dirq[j]);
  int s=inv(Z); *xo=md((long long)X*s); *yo=md((long long)Y*s); }
static int V[5][8]; static long nfound=0,ncheck=0;
static int cmpint(const void*a,const void*b){ return (*(int*)a)-(*(int*)b); }
static void checkAll(){
  int ld[20],lv[20],n=0; for(int d=0;d<5;d++) for(int i=0;i<sz[d];i++){ ld[n]=d; lv[n]=V[d][i]; n++; }
  int keys[200],l1_[200],l2_[200],np=0;
  for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if(ld[i]==ld[j]) continue; keys[np]=meetKey(ld[i],lv[i],ld[j],lv[j]); l1_[np]=i; l2_[np]=j; np++; }
  int idx[200]; for(int i=0;i<np;i++) idx[i]=i;
  for(int i=1;i<np;i++){ int k=idx[i]; int j=i-1; while(j>=0 && keys[idx[j]]>keys[k]){ idx[j+1]=idx[j]; j--; } idx[j+1]=k; }
  int b2=0,cnt[20],maxm=0,onL[20]; for(int i=0;i<n;i++) cnt[i]=1; int i=0;
  while(i<np){ int j=i; while(j<np && keys[idx[j]]==keys[idx[i]]) j++; int pairs=j-i; int m=2; while(m*(m-1)/2<pairs) m++; if(m*(m-1)/2!=pairs) return;
    b2+=m-1; if(m>maxm) maxm=m; memset(onL,0,sizeof(onL)); for(int t=i;t<j;t++){ onL[l1_[idx[t]]]=1; onL[l2_[idx[t]]]=1; } for(int t=0;t<n;t++) if(onL[t]) cnt[t]++; i=j; }
  if(b2!=64 || maxm>6) return;
  for(int t=0;t<n;t++) if(cnt[t]<4||cnt[t]>8) return;
  nfound++;
  printf("FOUND5 sizes %d %d %d %d %d p %d l1 %d l2 %d",sz[0],sz[1],sz[2],sz[3],sz[4],p,l1,l2);
  const char* nm[5]={"A","B","C","E","D"}; for(int d=0;d<5;d++){ printf(" %s",nm[d]); for(int t=0;t<sz[d];t++) printf(" %d",V[d][t]); }
  printf(" perline"); for(int t=0;t<n;t++) printf(" %d",cnt[t]); printf("\n"); }
static void tryD(){
  int vals[200],nv=0,x,y; int P=dirp[4],Q=dirq[4];
  for(int a=0;a<4;a++) for(int b=a+1;b<4;b++) for(int i=0;i<sz[a];i++) for(int j=0;j<sz[b];j++){ meetXY(a,V[a][i],b,V[b][j],&x,&y); vals[nv++]=md((long long)Q*x-(long long)P*y); }
  qsort(vals,nv,sizeof(int),cmpint);
  int cand[200],nc=0; int need=9-sz[4]; if(need<3) need=3;
  for(int i=0;i<nv;){ int j=i; while(j<nv && vals[j]==vals[i]) j++; if(j-i>=need) cand[nc++]=vals[i]; i=j; }
  int sD=sz[4]; if(nc<sD) return;
  int sel[8]; for(int i=0;i<sD;i++) sel[i]=i;
  while(1){ for(int i=0;i<sD;i++) V[4][i]=cand[sel[i]]; ncheck++; checkAll();
    int i=sD-1; while(i>=0 && sel[i]==nc-sD+i) i--; if(i<0) break; sel[i]++; for(int j=i+1;j<sD;j++) sel[j]=sel[j-1]+1; } }
/* enumerate subsets for pencil d (d=2,3: all of F_p; d=0: {0,1}+rest; d=1: {0}+rest) recursively */
static void enumPencil(int d){
  if(d==4){ tryD(); return; }
  int s=sz[d]; int sel[8]; int lo, fixed;
  if(d==0){ V[0][0]=0; V[0][1]=1; fixed=2; lo=2; } else if(d==1){ V[1][0]=0; fixed=1; lo=1; } else { fixed=0; lo=0; }
  int k=s-fixed; if(k==0){ enumPencil(d+1); return; }
  for(int i=0;i<k;i++) sel[i]=lo+i;
  while(1){ for(int i=0;i<k;i++) V[d][fixed+i]=sel[i]; enumPencil(d+1);
    int i=k-1; while(i>=0 && sel[i]==p-k+i) i--; if(i<0) break; sel[i]++; for(int j=i+1;j<k;j++) sel[j]=sel[j-1]+1; } }
int main(int argc,char**argv){
  if(argc<9){fprintf(stderr,"usage: pencils5_fp sA sB sC sE sD p l1 l2\n");return 1;}
  for(int i=0;i<5;i++) sz[i]=atoi(argv[1+i]); p=atoi(argv[6]); l1=md(atoi(argv[7])); l2=md(atoi(argv[8]));
  dirp[3]=l1; dirp[4]=l2;
  if(sz[0]<2){fprintf(stderr,"sA must be >= 2\n"); return 1;}
  enumPencil(0);
  fprintf(stderr,"sizes %d %d %d %d %d p=%d l=(%d,%d): %ld D-choices, found %ld\n",sz[0],sz[1],sz[2],sz[3],sz[4],p,l1,l2,ncheck,nfound); return 0; }
