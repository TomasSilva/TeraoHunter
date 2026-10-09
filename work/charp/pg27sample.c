/* exhaustive scan of n-subsets of the 31 lines of PG(2,5): keep lattices with chibar = (t-d1)(t-d2), integer d1 > 5, 5 not dividing n.
   One representative per invariant key (n, d1, point profile, sorted per-line multiplicity tuples). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define Q 7
#define NP 57
static int pts[NP][3], lines[NP][3];
static int onL[NP][Q+1]; /* points on each line */
static int normalize(int *v){ for(int i=0;i<3;i++) if(v[i]){ int inv=1; for(int k=1;k<Q;k++) if((v[i]*k)%Q==1) inv=k; for(int j=0;j<3;j++) v[j]=(v[j]*inv)%Q; return 1;} return 0; }
typedef struct { uint64_t h; int n; int d1; int subset[NP]; char key[512]; } Entry;
#define HSIZE (1<<22)
static Entry *table;
static long nkeys=0;
static uint64_t fnv(const char *s){ uint64_t h=1469598103934665603ULL; while(*s){ h^=(unsigned char)*s++; h*=1099511628211ULL;} return h; }
static void insert(const char *key, int n, int d1, int *S){
  uint64_t h=fnv(key); long i=h&(HSIZE-1);
  while(table[i].h){ if(table[i].h==h && !strcmp(table[i].key,key)) return; i=(i+1)&(HSIZE-1); }
  table[i].h=h; table[i].n=n; table[i].d1=d1; strcpy(table[i].key,key); memcpy(table[i].subset,S,sizeof(int)*n); nkeys++;
}
int main(int argc, char **argv){
  int np=0;
  for(int a=0;a<Q;a++)for(int b=0;b<Q;b++)for(int c=0;c<Q;c++){ int v[3]={a,b,c}; if(!normalize(v)) continue; int dup=0; for(int i=0;i<np;i++) if(pts[i][0]==v[0]&&pts[i][1]==v[1]&&pts[i][2]==v[2]) dup=1; if(!dup){ memcpy(pts[np],v,sizeof v); np++; } }
  memcpy(lines,pts,sizeof pts);
  for(int l=0;l<NP;l++){ int k=0; for(int p=0;p<NP;p++){ if((pts[p][0]*lines[l][0]+pts[p][1]*lines[l][1]+pts[p][2]*lines[l][2])%Q==0) onL[l][k++]=p; } if(k!=Q+1){fprintf(stderr,"bad incidence\n");return 1;} }
  table=calloc(HSIZE,sizeof(Entry));
  int n=atoi(argv[1]); long nsamp=atol(argv[2]); srand48(12345);
  { long cnt=0, hits=0;
    for(long it=0; it<nsamp; it++){
      int S[NP]; int used[NP]; memset(used,0,sizeof used);
      for(int i=0;i<n;i++){ int l; do { l = (int)(drand48()*NP); } while(used[l]); used[l]=1; S[i]=l; }
      /* sort S */ for(int a=0;a<n;a++)for(int b=a+1;b<n;b++) if(S[b]<S[a]){int t=S[a];S[a]=S[b];S[b]=t;}
      int mult[NP]; memset(mult,0,sizeof mult);
      for(int i=0;i<n;i++){ int l=S[i]; for(int k=0;k<=Q;k++) mult[onL[l][k]]++; }
      int b2=0; for(int p=0;p<NP;p++) if(mult[p]>=2) b2+=mult[p]-1;
      cnt++;
      int disc=(n-1)*(n-1)-4*(b2-n+1);
      if(disc>=0){ int s=0; while((s+1)*(s+1)<=disc) s++;
        if(s*s==disc && ((n-1-s)%2==0)){ int d1=(n-1-s)/2;
          if(d1>Q){ hits++;
            int prof[Q+2]; memset(prof,0,sizeof prof); for(int p=0;p<NP;p++) if(mult[p]>=2) prof[mult[p]]++;
            char tuples[NP][16]; int tl=0;
            for(int i=0;i<n;i++){ int l=S[i]; int ms[Q+1]; for(int k=0;k<=Q;k++) ms[k]=mult[onL[l][k]];
              for(int a=0;a<=Q;a++)for(int b=a+1;b<=Q;b++) if(ms[b]<ms[a]){int t=ms[a];ms[a]=ms[b];ms[b]=t;}
              for(int k=0;k<=Q;k++) tuples[tl][k]='0'+ms[k]; tuples[tl][Q+1]=0; tl++; }
            for(int a=0;a<tl;a++)for(int b=a+1;b<tl;b++) if(strcmp(tuples[b],tuples[a])<0){ char t[16]; strcpy(t,tuples[a]); strcpy(tuples[a],tuples[b]); strcpy(tuples[b],t); }
            char key[512]; int pos=sprintf(key,"n%d d%d p",n,d1); for(int m=2;m<=Q+1;m++) pos+=sprintf(key+pos,"%d.",prof[m]); key[pos++]='|'; key[pos]=0;
            for(int i=0;i<tl;i++){ strcat(key,tuples[i]); strcat(key,","); }
            insert(key,n,d1,S);
          } } }
    }
    fprintf(stderr,"n=%d: %ld samples, %ld hits, %ld keys\n",n,cnt,hits,nkeys);
  }
  for(long i=0;i<HSIZE;i++) if(table[i].h){ printf("%s ;", table[i].key); for(int j=0;j<table[i].n;j++) printf(" %d",table[i].subset[j]); printf("\n"); }
  return 0;
}
