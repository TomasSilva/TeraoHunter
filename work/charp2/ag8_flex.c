/* For every four-pencil hit of ag8_4pencil.c (sub-arrangement of PG(2,8) containing H, d1 >= 9):
   the offsets satisfy linear equations M(mu) s = 0 (one per extra line through an affine point):
     {1,2,3}: c + a + b = 0;  {1,2,4}: d + b + mu a = 0;  {1,3,4}: d + c + (1+mu) a = 0;
     {2,3,4}: d + (1+mu) b + mu c = 0.
   The realisation space of the lattice fibres over mu with fibre ker M(mu) (minus degenerations),
   mod translations (2-dim) and scaling. rank at mu = t (F_8) vs generic mu (primitive elements of
   F_{2^16}): if K - rank_gen >= 3 the lattice has realisations with mu off F_8 (mu-flexible).
   Usage: ag8_flex part nparts   (splits the (T1,T2) loop). Prints the flexible hits and a histogram. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
static int mul8(int a, int b) { int r = 0; for (int i = 0; i < 3; i++) if (b >> i & 1) r ^= a << i;
  for (int i = 4; i >= 3; i--) if (r >> i & 1) r ^= 0xB << (i - 3); return r; }
static int inv8[8];
static uint16_t EXP[2 * 65536], LOG[65536];
static inline uint16_t m16(uint16_t a, uint16_t b) { return (a && b) ? EXP[LOG[a] + LOG[b]] : 0; }
static inline uint16_t i16(uint16_t a) { return EXP[65535 - LOG[a]]; }
/* rank of E x K matrix over F_8 (field=0) or F_2^16 (field=1); entries given as (c0 + c1*mu) */
static int rank_mat(int E, int K, unsigned char co[][32][2], int field, int mu) {
  static int M[160][32];
  for (int i = 0; i < E; i++) for (int j = 0; j < K; j++) {
    int c0 = co[i][j][0], c1 = co[i][j][1];
    M[i][j] = field == 0 ? (c0 ^ (c1 ? mu : 0)) : (c0 ^ (c1 ? mu : 0));
  }
  int r = 0;
  for (int col = 0; col < K && r < E; col++) {
    int piv = -1; for (int i = r; i < E; i++) if (M[i][col]) { piv = i; break; }
    if (piv < 0) continue;
    for (int j = 0; j < K; j++) { int tmp = M[r][j]; M[r][j] = M[piv][j]; M[piv][j] = tmp; }
    int iv = field == 0 ? inv8[M[r][col]] : i16(M[r][col]);
    for (int j = 0; j < K; j++) M[r][j] = field == 0 ? mul8(M[r][j], iv) : m16(M[r][j], iv);
    for (int i = 0; i < E; i++) if (i != r && M[i][col]) {
      int f = M[i][col];
      for (int j = 0; j < K; j++) M[i][j] ^= field == 0 ? mul8(f, M[r][j]) : m16(f, M[r][j]);
    }
    r++;
  }
  return r;
}
int main(int argc, char **argv) {
  int part = atoi(argv[1]), nparts = atoi(argv[2]);
  for (int a = 1; a < 8; a++) for (int b = 1; b < 8; b++) if (mul8(a, b) == 1) inv8[a] = b;
  { uint32_t x = 1; for (int i = 0; i < 65535; i++) { EXP[i] = x; LOG[x] = i; x <<= 1; if (x & 0x10000) x ^= 0x1100B; }
    for (int i = 65535; i < 2 * 65536; i++) EXP[i] = EXP[i - 65535]; }
  int t = 2, mu1 = 2, mu2 = EXP[7];
  static uint64_t m12[65536], m34[65536];
  int v4[64];
  for (int x = 0; x < 8; x++) for (int y = 0; y < 8; y++) v4[8 * x + y] = x ^ mul8(t, y);
  for (int T1 = 0; T1 < 256; T1++) for (int T2 = 0; T2 < 256; T2++) {
    uint64_t m = 0, n = 0;
    for (int x = 0; x < 8; x++) for (int y = 0; y < 8; y++) {
      int P = 8 * x + y;
      if ((T1 >> y & 1) && (T2 >> x & 1)) m |= 1ULL << P;
      if ((T1 >> (x ^ y) & 1) && (T2 >> v4[P] & 1)) n |= 1ULL << P;
    }
    m12[T1 * 256 + T2] = m; m34[T1 * 256 + T2] = n;
  }
  static unsigned char co[160][32][2];
  long hist[40][40] = {{0}}; long hits = 0, flex = 0;
  for (int A = part; A < 65536; A += nparts) {
    int T1 = A >> 8, T2 = A & 255;
    if (T1 == 255 || T2 == 255) continue;
    int tau12 = __builtin_popcount(T1) + __builtin_popcount(T2);
    for (int B = 0; B < 65536; B++) {
      int T3 = B >> 8, T4 = B & 255;
      if (T3 == 255 || T4 == 255) continue;
      int tau = tau12 + __builtin_popcount(T3) + __builtin_popcount(T4);
      int s = 16 - tau; if (s < 2) continue;
      int W = __builtin_popcountll(m12[A] & m34[B]), ok = 0;
      for (int u = 1; 2 * u <= s; u++) if (u * (s - u) == W) ok = 1;
      if (!ok) continue;
      hits++;
      int S[4] = {255 ^ T1, 255 ^ T2, 255 ^ T3, 255 ^ T4}, col[4][8], K = 0;
      for (int j = 0; j < 4; j++) for (int v = 0; v < 8; v++) col[j][v] = (S[j] >> v & 1) ? K++ : -1;
      int E = 0;
      for (int x = 0; x < 8; x++) for (int y = 0; y < 8; y++) {
        int P = 8 * x + y, vv[4] = {y, x, x ^ y, v4[P]}, in[4];
        for (int j = 0; j < 4; j++) in[j] = col[j][vv[j]];
        int cnt = (in[0] >= 0) + (in[1] >= 0) + (in[2] >= 0) + (in[3] >= 0);
        if (cnt < 3) continue;
        #define NEWROW() do { for (int q = 0; q < K; q++) co[E][q][0] = co[E][q][1] = 0; } while (0)
        if (in[0] >= 0 && in[1] >= 0) {
          if (in[2] >= 0) { NEWROW(); co[E][in[2]][0] = 1; co[E][in[0]][0] = 1; co[E][in[1]][0] = 1; E++; }
          if (in[3] >= 0) { NEWROW(); co[E][in[3]][0] = 1; co[E][in[1]][0] = 1; co[E][in[0]][1] = 1; E++; }
        } else if (in[0] >= 0) { NEWROW(); co[E][in[3]][0] = 1; co[E][in[2]][0] = 1; co[E][in[0]][0] = 1; co[E][in[0]][1] = 1; E++; }
        else { NEWROW(); co[E][in[3]][0] = 1; co[E][in[1]][0] = 1; co[E][in[1]][1] = 1; co[E][in[2]][1] = 1; E++; }
      }
      int rt = rank_mat(E, K, co, 0, t);
      int rg = rank_mat(E, K, co, 1, mu1), rg2 = rank_mat(E, K, co, 1, mu2);
      if (rg2 > rg) rg = rg2;
      hist[K - rt][K - rg]++;
      if (K - rg >= 3) { flex++; printf("FLEX %02x %02x %02x %02x K %d W %d dim_t %d dim_gen %d\n", S[0], S[1], S[2], S[3], K, W, K - rt, K - rg); }
    }
  }
  fprintf(stderr, "part %d hits %ld flex %ld\n", part, hits, flex);
  for (int i = 0; i < 40; i++) for (int j = 0; j < 40; j++) if (hist[i][j]) fprintf(stderr, "HIST dim_t %d dim_gen %d : %ld\n", i, j, hist[i][j]);
  return 0;
}
