/* pencils_fp.c -- complete census over F_p of 4-pencil line arrangements with collinear centres
 * (1:0:0),(0:1:0),(1:1:0),(mu:1:0) on H={z=0}; pencil sizes sA sB sC sD; D derived from the A,B,C points.
 * Normalisation: a0=0, a1=1 (if sA>=2), b0=0.  Output: b2=80 (with H), every line 4..8 points.
 * usage: pencils_fp sA sB sC sD p mu
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int sA,sB,sC,sD,p,mu;
static int dirp[4]={1,0,1,0}, dirq[4]={0,1,1,0};
static int md(long long a){ a%=p; if(a<0)a+=p; return (int)a; }
static int inv(int a){ int r=1,b=a%p,e=p-2; if(b<0)b+=p; while(e){ if(e&1) r=(int)((long long)r*b%p); b=(int)((long long)b*b%p); e>>=1;} return r; }
/* point as normalised triple key = x*p*p + y*p + z with first nonzero coord = 1 */
static int meetKey(int i,int v,int j,int w){
  int X=md((long long)dirp[i]*w - (long long)v*dirp[j]);
  int Y=md(-(long long)v*dirq[j] + (long long)dirq[i]*w);
  int Z=md(-(long long)dirq[i]*dirp[j] + (long long)dirp[i]*dirq[j]);
  int s;
  if(X){ s=inv(X); X=1; Y=md((long long)Y*s); Z=md((long long)Z*s);} else if(Y){ s=inv(Y); Y=1; Z=md((long long)Z*s);} else { Z=1; }
  return (X*p+Y)*p+Z;
}
static int meetXY(int i,int v,int j,int w,int*xo,int*yo){ /* affine coordinates (Z != 0 guaranteed for different directions) */
  int X=md((long long)dirp[i]*w - (long long)v*dirp[j]);
  int Y=md(-(long long)v*dirq[j] + (long long)dirq[i]*w);
  int Z=md(-(long long)dirq[i]*dirp[j] + (long long)dirp[i]*dirq[j]);
  int s=inv(Z); *xo=md((long long)X*s); *yo=md((long long)Y*s); return 1;
}
static int A[8],B[8],C[8],D[8]; static long nfound=0,ncheck=0;
static int cmpint(const void*a,const void*b){ return (*(int*)a)-(*(int*)b); }
static void checkAll(){
  int ld[20],lv[20],n=0;
  for(int i=0;i<sA;i++){ld[n]=0;lv[n]=A[i];n++;} for(int i=0;i<sB;i++){ld[n]=1;lv[n]=B[i];n++;}
  for(int i=0;i<sC;i++){ld[n]=2;lv[n]=C[i];n++;} for(int i=0;i<sD;i++){ld[n]=3;lv[n]=D[i];n++;}
  int keys[200],l1[200],l2[200],np=0;
  for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if(ld[i]==ld[j]) continue; keys[np]=meetKey(ld[i],lv[i],ld[j],lv[j]); l1[np]=i; l2[np]=j; np++; }
  int idx[200]; for(int i=0;i<np;i++) idx[i]=i;
  for(int i=1;i<np;i++){ int k=idx[i]; int j=i-1; while(j>=0 && keys[idx[j]]>keys[k]){ idx[j+1]=idx[j]; j--; } idx[j+1]=k; }
  int b2=0,cnt[20],maxm=0,onL[20]; for(int i=0;i<n;i++) cnt[i]=1;
  int i=0;
  while(i<np){ int j=i; while(j<np && keys[idx[j]]==keys[idx[i]]) j++; int pairs=j-i; int m=2; while(m*(m-1)/2<pairs) m++; if(m*(m-1)/2!=pairs) return;
    b2+=m-1; if(m>maxm) maxm=m; memset(onL,0,sizeof(onL)); for(int t=i;t<j;t++){ onL[l1[idx[t]]]=1; onL[l2[idx[t]]]=1; } for(int t=0;t<n;t++) if(onL[t]) cnt[t]++; i=j; }
  if(b2!=64 || maxm>5) return;
  for(int t=0;t<n;t++) if(cnt[t]<4||cnt[t]>8) return;
  nfound++;
  printf("FOUND sizes %d %d %d %d p %d mu %d A",sA,sB,sC,sD,p,mu); for(int t=0;t<sA;t++) printf(" %d",A[t]); printf(" B"); for(int t=0;t<sB;t++) printf(" %d",B[t]);
  printf(" C"); for(int t=0;t<sC;t++) printf(" %d",C[t]); printf(" D"); for(int t=0;t<sD;t++) printf(" %d",D[t]); printf(" perline"); for(int t=0;t<n;t++) printf(" %d",cnt[t]); printf("\n");
}
static void tryD(){
  int vals[100],nv=0,x,y; int P=dirp[3],Q=dirq[3];
  for(int i=0;i<sA;i++) for(int j=0;j<sB;j++){ meetXY(0,A[i],1,B[j],&x,&y); vals[nv++]=md((long long)Q*x-(long long)P*y);}
  for(int i=0;i<sA;i++) for(int j=0;j<sC;j++){ meetXY(0,A[i],2,C[j],&x,&y); vals[nv++]=md((long long)Q*x-(long long)P*y);}
  for(int i=0;i<sB;i++) for(int j=0;j<sC;j++){ meetXY(1,B[i],2,C[j],&x,&y); vals[nv++]=md((long long)Q*x-(long long)P*y);}
  qsort(vals,nv,sizeof(int),cmpint);
  int cand[100],nc=0; int need = 9-sD; if(need<3) need=3; /* entries >= 2k+t >= 9-sD, and >= 3 */
  for(int i=0;i<nv;){ int j=i; while(j<nv && vals[j]==vals[i]) j++; if(j-i>=need) cand[nc++]=vals[i]; i=j; }
  if(nc<sD) return;
  int sel[8]; for(int i=0;i<sD;i++) sel[i]=i;
  while(1){ for(int i=0;i<sD;i++) D[i]=cand[sel[i]]; ncheck++; checkAll();
    int i=sD-1; while(i>=0 && sel[i]==nc-sD+i) i--; if(i<0) break; sel[i]++; for(int j=i+1;j<sD;j++) sel[j]=sel[j-1]+1; }
}
static void enumC(){ int sel[8]; for(int i=0;i<sC;i++) sel[i]=i;
  while(1){ for(int i=0;i<sC;i++) C[i]=sel[i]; tryD();
    int i=sC-1; while(i>=0 && sel[i]==p-sC+i) i--; if(i<0) break; sel[i]++; for(int j=i+1;j<sC;j++) sel[j]=sel[j-1]+1; } }
int main(int argc,char**argv){
  if(argc<7){fprintf(stderr,"usage: pencils_fp sA sB sC sD p mu\n");return 1;}
  sA=atoi(argv[1]);sB=atoi(argv[2]);sC=atoi(argv[3]);sD=atoi(argv[4]);p=atoi(argv[5]);mu=md(atoi(argv[6]));
  dirp[3]=mu; dirq[3]=1;
  /* A = {0,1} u (sA-2)-subset of [2,p-1]; B = {0} u (sB-1)-subset of [1,p-1] */
  int selA[8],selB[8]; int ka=sA-2, kb=sB-1;
  for(int i=0;i<ka;i++) selA[i]=2+i;
  while(1){ A[0]=0; A[1]=1; for(int i=0;i<ka;i++) A[2+i]=selA[i];
    for(int i=0;i<kb;i++) selB[i]=1+i;
    while(1){ B[0]=0; for(int i=0;i<kb;i++) B[1+i]=selB[i]; enumC();
      int i=kb-1; while(i>=0 && selB[i]==p-kb+i) i--; if(i<0) break; selB[i]++; for(int j=i+1;j<kb;j++) selB[j]=selB[j-1]+1; }
    int i=ka-1; while(i>=0 && selA[i]==p-ka+i) i--; if(i<0) break; selA[i]++; for(int j=i+1;j<ka;j++) selA[j]=selA[j-1]+1; }
  fprintf(stderr,"sizes %d %d %d %d p=%d mu=%d: %ld D-choices, found %ld\n",sA,sB,sC,sD,p,mu,ncheck,nfound); return 0; }
