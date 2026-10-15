/* Orbit representatives of the four-pencil PG(2,8) hits (d1 >= 9) under the 5376 semilinear affine
   maps preserving the four directions, and, for each representative, every mu in F_4, F_16, F_32,
   F_64 \ F_8, F_128, F_256 at which the incidence system M(mu) s = 0 has a non-translation solution
   (rank <= K-3): candidate second realisations of the lattice with cross-ratio off F_8. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
static int mul8(int a, int b) { int r = 0; for (int i = 0; i < 3; i++) if (b >> i & 1) r ^= a << i;
  for (int i = 4; i >= 3; i--) if (r >> i & 1) r ^= 0xB << (i - 3); return r; }
static int sq8(int a) { return mul8(a, a); }
/* GF(2^k) tables */
typedef struct { int k, q; int exp[1024], log[512]; } GF;
static void gf_init(GF *F, int k, int poly) { F->k = k; F->q = 1 << k; int x = 1;
  for (int i = 0; i < F->q - 1; i++) { F->exp[i] = x; F->log[x] = i; x <<= 1; if (x & F->q) x ^= poly; }
  for (int i = F->q - 1; i < 2 * F->q; i++) F->exp[i] = F->exp[i - (F->q - 1)]; }
static inline int gmul(GF *F, int a, int b) { return (a && b) ? F->exp[F->log[a] + F->log[b]] : 0; }
static inline int ginv(GF *F, int a) { return F->exp[F->q - 1 - F->log[a]]; }
static int rankF(GF *F, int E, int K, unsigned char co[][32][2], int mu) {
  static int M[200][32];
  for (int i = 0; i < E; i++) for (int j = 0; j < K; j++) M[i][j] = co[i][j][0] ^ (co[i][j][1] ? mu : 0);
  int r = 0;
  for (int col = 0; col < K && r < E; col++) {
    int piv = -1; for (int i = r; i < E; i++) if (M[i][col]) { piv = i; break; }
    if (piv < 0) continue;
    for (int j = 0; j < K; j++) { int tmp = M[r][j]; M[r][j] = M[piv][j]; M[piv][j] = tmp; }
    int iv = ginv(F, M[r][col]);
    for (int j = 0; j < K; j++) M[r][j] = gmul(F, M[r][j], iv);
    for (int i = 0; i < E; i++) if (i != r && M[i][col]) { int f = M[i][col]; for (int j = 0; j < K; j++) M[i][j] ^= gmul(F, f, M[r][j]); }
    r++;
  }
  return r;
}
/* directions: linear forms l_j(x,y): l1 = y, l2 = x, l3 = x+y, l4 = x + t y (t = 2) */
static int lform(int j, int x, int y) { return j == 0 ? y : j == 1 ? x : j == 2 ? (x ^ y) : (x ^ mul8(2, y)); }
int main(void) {
  /* semilinear maps P -> M sigma^e(P): find those permuting the four line directions */
  int nL = 0; static int perm[200][4], offmap[200][4][8];
  for (int e = 0; e < 3; e++)
  for (int m00 = 0; m00 < 8; m00++) for (int m01 = 0; m01 < 8; m01++) for (int m10 = 0; m10 < 8; m10++) for (int m11 = 0; m11 < 8; m11++) {
    if ((mul8(m00, m11) ^ mul8(m01, m10)) == 0) continue;
    int ok = 1, pm[4], om[4][8];
    for (int j = 0; j < 4 && ok; j++) {
      /* image of the line l_j = v: map its 8 points, check they lie on a line l_j' = v' */
      for (int v = 0; v < 8 && ok; v++) {
        int jj = -1, vv = -1;
        for (int x = 0; x < 8 && ok; x++) for (int y = 0; y < 8 && ok; y++) {
          if (lform(j, x, y) != v) continue;
          int sx = x, sy = y; for (int i = 0; i < e; i++) { sx = sq8(sx); sy = sq8(sy); }
          int X = mul8(m00, sx) ^ mul8(m01, sy), Y = mul8(m10, sx) ^ mul8(m11, sy);
          if (jj < 0) { /* find a direction whose form is constant on the image: test the next point later */
            for (int c = 0; c < 4; c++) { /* tentatively */ }
            /* store first image point */
            static int fx, fy; fx = X; fy = Y; jj = 100; vv = (fx << 3) | fy; continue; }
          if (jj == 100) { int fx = vv >> 3, fy = vv & 7; int c;
            for (c = 0; c < 4; c++) if (lform(c, fx, fy) == lform(c, X, Y)) break;
            if (c == 4) { ok = 0; break; }
            jj = c; vv = lform(c, X, Y); }
          else if (lform(jj, X, Y) != vv) ok = 0;
        }
        if (!ok) break;
        if (v == 0) pm[j] = jj; else if (pm[j] != jj) { ok = 0; break; }
        om[j][v] = vv;
      }
    }
    if (!ok) continue;
    int used = 0; for (int j = 0; j < 4; j++) used |= 1 << pm[j]; if (used != 15) continue;
    memcpy(perm[nL], pm, sizeof pm); memcpy(offmap[nL], om, sizeof om); nL++;
  }
  fprintf(stderr, "semilinear direction-preserving maps: %d (expect 84)\n", nL);
  /* mask tables: image mask of pencil j under map g, then translation by (tx,ty) adds l_{j'}(tx,ty) */
  static unsigned char img[200][4][256], tr[4][64][256];
  for (int g = 0; g < nL; g++) for (int j = 0; j < 4; j++) for (int m = 0; m < 256; m++) {
    int r = 0; for (int v = 0; v < 8; v++) if (m >> v & 1) r |= 1 << offmap[g][j][v]; img[g][j][m] = r; }
  for (int j = 0; j < 4; j++) for (int T = 0; T < 64; T++) for (int m = 0; m < 256; m++) {
    int s = lform(j, T >> 3, T & 7), r = 0; for (int v = 0; v < 8; v++) if (m >> v & 1) r |= 1 << (v ^ s); tr[j][T][m] = r; }
  /* fields */
  static GF F4, F16, F32, F64, F128, F256, F8b;
  gf_init(&F16, 4, 0x13); gf_init(&F32, 5, 0x25); gf_init(&F64, 6, 0x43); gf_init(&F128, 7, 0x89); gf_init(&F256, 8, 0x11D);
  gf_init(&F8b, 3, 0xB);
  GF *fields[5] = {&F16, &F32, &F64, &F128, &F256};
  static uint64_t m12[65536], m34[65536]; int v4[64];
  for (int x = 0; x < 8; x++) for (int y = 0; y < 8; y++) v4[8 * x + y] = x ^ mul8(2, y);
  for (int T1 = 0; T1 < 256; T1++) for (int T2 = 0; T2 < 256; T2++) {
    uint64_t m = 0, n = 0;
    for (int x = 0; x < 8; x++) for (int y = 0; y < 8; y++) { int P = 8 * x + y;
      if ((T1 >> y & 1) && (T2 >> x & 1)) m |= 1ULL << P;
      if ((T1 >> (x ^ y) & 1) && (T2 >> v4[P] & 1)) n |= 1ULL << P; }
    m12[T1 * 256 + T2] = m; m34[T1 * 256 + T2] = n; }
  static unsigned char co[200][32][2];
  long hits = 0, reps = 0, cand = 0;
  for (int A = 0; A < 65536; A++) {
    int T1 = A >> 8, T2 = A & 255; if (T1 == 255 || T2 == 255) continue;
    int tau12 = __builtin_popcount(T1) + __builtin_popcount(T2);
    for (int B = 0; B < 65536; B++) {
      int T3 = B >> 8, T4 = B & 255; if (T3 == 255 || T4 == 255) continue;
      int tau = tau12 + __builtin_popcount(T3) + __builtin_popcount(T4), s = 16 - tau; if (s < 2) continue;
      int W = __builtin_popcountll(m12[A] & m34[B]), ok = 0, uu = 0;
      for (int u = 1; 2 * u <= s; u++) if (u * (s - u) == W) { ok = 1; uu = u; }
      if (!ok) continue;
      hits++;
      int S[4] = {255 ^ T1, 255 ^ T2, 255 ^ T3, 255 ^ T4};
      uint32_t key = (uint32_t)S[0] << 24 | S[1] << 16 | S[2] << 8 | S[3];
      int minimal = 1;
      for (int g = 0; g < nL && minimal; g++) {
        int I[4]; for (int j = 0; j < 4; j++) I[perm[g][j]] = img[g][j][S[j]];
        for (int T = 0; T < 64; T++) {
          uint32_t k2 = (uint32_t)tr[0][T][I[0]] << 24 | tr[1][T][I[1]] << 16 | tr[2][T][I[2]] << 8 | tr[3][T][I[3]];
          if (k2 < key) { minimal = 0; break; }
        }
      }
      if (!minimal) continue;
      reps++;
      /* equations */
      int col[4][8], K = 0;
      for (int j = 0; j < 4; j++) for (int v = 0; v < 8; v++) col[j][v] = (S[j] >> v & 1) ? K++ : -1;
      int E = 0;
      for (int x = 0; x < 8; x++) for (int y = 0; y < 8; y++) {
        int P = 8 * x + y, vv[4] = {y, x, x ^ y, v4[P]}, in[4];
        for (int j = 0; j < 4; j++) in[j] = col[j][vv[j]];
        int cnt = (in[0] >= 0) + (in[1] >= 0) + (in[2] >= 0) + (in[3] >= 0); if (cnt < 3) continue;
        #define NEWROW() do { for (int q = 0; q < K; q++) co[E][q][0] = co[E][q][1] = 0; } while (0)
        if (in[0] >= 0 && in[1] >= 0) {
          if (in[2] >= 0) { NEWROW(); co[E][in[2]][0] = 1; co[E][in[0]][0] = 1; co[E][in[1]][0] = 1; E++; }
          if (in[3] >= 0) { NEWROW(); co[E][in[3]][0] = 1; co[E][in[1]][0] = 1; co[E][in[0]][1] = 1; E++; }
        } else if (in[0] >= 0) { NEWROW(); co[E][in[3]][0] = 1; co[E][in[2]][0] = 1; co[E][in[0]][0] = 1; co[E][in[0]][1] = 1; E++; }
        else { NEWROW(); co[E][in[3]][0] = 1; co[E][in[1]][0] = 1; co[E][in[1]][1] = 1; co[E][in[2]][1] = 1; E++; }
      }
      int rt = rankF(&F8b, E, K, co, 2);
      if (K - rt != 3) fprintf(stderr, "WARN dim_t %d for %08x\n", K - rt, key);
      int found = 0;
      for (int fi = 0; fi < 5; fi++) { GF *F = fields[fi];
        for (int mu = 2; mu < F->q; mu++) {
          /* skip mu in F_4 when F = F_16/F_256 except once, skip F_8 inside F_64: test minimal polynomial degree */
          int d = 1, z = gmul(F, mu, mu); while (z != mu) { z = gmul(F, z, z); d++; }
          if (d == 3) continue;                        /* F_8 elements: the hit's own field */
          if (d == 2 && F->k != 4) continue;           /* F_4 tested once, inside F_16 */
          if (d < F->k && d != 2) continue;            /* proper subfield elements tested in their own field */
          /* one element per Frobenius orbit: keep mu if it is the smallest in its orbit */
          int smallest = 1; z = mu; for (int i = 1; i < d; i++) { z = gmul(F, z, z); if (z < mu) smallest = 0; }
          if (!smallest) continue;
          int r = rankF(F, E, K, co, mu);
          if (K - r >= 3) { printf("CAND %08x K %d exps %d %d field 2^%d mu %d deg %d dim %d\n", key, K, 8 + uu, 8 + s - uu, F->k, mu, d, K - r); found = 1; }
        } }
      if (found) cand++;
    }
  }
  fprintf(stderr, "hits %ld orbit representatives %ld with candidate second realisations %ld\n", hits, reps, cand);
  return 0;
}
