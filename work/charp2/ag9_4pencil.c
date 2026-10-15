/* Characteristic 3 analogue of ag8_orbits.c: four-pencil sub-arrangements of PG(2,9) containing H,
   directions (1:0),(0:1),(1:1),(i:1), F_9 = F_3[i]/(i^2+1) (every generic 4-subset of P^1(F_9) is
   equivalent to this one).  Lines y = a, x = b, x - y = c, x - i y = d.  For A in PG(2,9) with H in A:
   exponents (9+u, 9+v), u+v = K-18, uv = #uncovered affine points; keep u >= 1 (d1 >= 10 > 9 = q, so
   every hit is non-free by the Frobenius bound).  For each orbit representative (semilinear affine maps
   preserving the directions): rank of the incidence system at mu = i, at generic mu (F_{3^7}), and at
   every mu (one per Frobenius orbit) of F_27, F_81 \ F_9, F_243, F_729 \ F_9.
   Usage: ag9_4pencil part nparts */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef unsigned __int128 u128;
static int add9(int u, int v) { return ((u % 3 + v % 3) % 3) + 3 * ((u / 3 + v / 3) % 3); }
static int neg9(int u) { return ((3 - u % 3) % 3) + 3 * ((3 - u / 3) % 3); }
static int mul9(int u, int v) { int a = u % 3, b = u / 3, c = v % 3, d = v / 3;
  int re = ((a * c - b * d) % 3 + 3) % 3, im = (a * d + b * c) % 3; return re + 3 * im; }
static int lform(int j, int x, int y) { /* l1 = y, l2 = x, l3 = x - y, l4 = x - i y (i = 3) */
  return j == 0 ? y : j == 1 ? x : j == 2 ? add9(x, neg9(y)) : add9(x, neg9(mul9(3, y))); }
/* GF(3^k): elements as base-3 integers, exp/log + add tables */
typedef struct { int k, q; int *exp, *log; short *add; int *neg; } GF;
static int polymulx(int e, int k, int *f) { /* multiply element (base-3 digits) by x mod monic f */
  int d[16], top; for (int i = 0; i < k; i++) { d[i] = e % 3; e /= 3; }
  top = d[k - 1]; for (int i = k - 1; i > 0; i--) d[i] = d[i - 1]; d[0] = 0;
  for (int i = 0; i < k; i++) d[i] = ((d[i] - top * f[i]) % 3 + 3) % 3;
  int r = 0; for (int i = k - 1; i >= 0; i--) r = 3 * r + d[i]; return r; }
static int digadd(int u, int v, int k) { int r = 0, p = 1; for (int i = 0; i < k; i++) { r += p * ((u % 3 + v % 3) % 3); u /= 3; v /= 3; p *= 3; } return r; }
static void gf_init(GF *F, int k) {
  F->k = k; F->q = 1; for (int i = 0; i < k; i++) F->q *= 3;
  F->exp = malloc(sizeof(int) * 2 * F->q); F->log = malloc(sizeof(int) * F->q);
  int f[16];
  for (int c = 0; ; c++) { /* search a primitive monic polynomial x^k + f_{k-1} x^{k-1} + ... + f_0 */
    int cc = c; for (int i = 0; i < k; i++) { f[i] = cc % 3; cc /= 3; }
    if (f[0] == 0) continue;
    int x = 1, ord = 0; do { x = polymulx(x, k, f); ord++; } while (x != 1 && ord < F->q);
    if (ord == F->q - 1) break;
  }
  int x = 1; for (int i = 0; i < F->q - 1; i++) { F->exp[i] = x; F->log[x] = i; x = polymulx(x, k, f); }
  for (int i = F->q - 1; i < 2 * F->q; i++) F->exp[i] = F->exp[i - (F->q - 1)];
  F->add = malloc(sizeof(short) * F->q * F->q); F->neg = malloc(sizeof(int) * F->q);
  for (int u = 0; u < F->q; u++) for (int v = 0; v < F->q; v++) F->add[u * F->q + v] = digadd(u, v, k);
  for (int u = 0; u < F->q; u++) F->neg[u] = digadd(u, u, k);   /* -u = 2u */
}
static inline int gm(GF *F, int a, int b) { return (a && b) ? F->exp[F->log[a] + F->log[b]] : 0; }
static inline int ga(GF *F, int a, int b) { return F->add[a * F->q + b]; }
static inline int gi(GF *F, int a) { return F->exp[F->q - 1 - F->log[a]]; }
/* coefficient c0 + c1 mu (c0,c1 in F_3 as 0,1,2) evaluated in F: */
static int coefval(GF *F, int c0, int c1, int mu) {
  int one = 1, r = 0;
  for (int i = 0; i < c0; i++) r = ga(F, r, one);
  for (int i = 0; i < c1; i++) r = ga(F, r, mu);
  return r; }
static int rankF(GF *F, int E, int K, unsigned char co[][40][2], int mu) {
  static int M[300][40];
  for (int i = 0; i < E; i++) for (int j = 0; j < K; j++) M[i][j] = coefval(F, co[i][j][0], co[i][j][1], mu);
  int r = 0;
  for (int col = 0; col < K && r < E; col++) {
    int piv = -1; for (int i = r; i < E; i++) if (M[i][col]) { piv = i; break; }
    if (piv < 0) continue;
    for (int j = 0; j < K; j++) { int tmp = M[r][j]; M[r][j] = M[piv][j]; M[piv][j] = tmp; }
    int iv = gi(F, M[r][col]);
    for (int j = 0; j < K; j++) M[r][j] = gm(F, M[r][j], iv);
    for (int i = 0; i < E; i++) if (i != r && M[i][col]) { int f = F->neg[M[i][col]]; for (int j = 0; j < K; j++) M[i][j] = ga(F, M[i][j], gm(F, f, M[r][j])); }
    r++;
  }
  return r;
}
int main(int argc, char **argv) {
  int part = atoi(argv[1]), nparts = atoi(argv[2]);
  /* F_9 inside the GF machinery: build GF(9) with i <-> element; we need mu = i as an element of a GF(3^2)
     table: map our encoding (a + 3b <-> a + b i) to the table encoding via the minimal polynomial of i */
  /* semilinear maps preserving directions */
  static int perm[400][4], offmap[400][4][9]; int nL = 0;
  for (int e = 0; e < 2; e++)
  for (int m00 = 0; m00 < 9; m00++) for (int m01 = 0; m01 < 9; m01++) for (int m10 = 0; m10 < 9; m10++) for (int m11 = 0; m11 < 9; m11++) {
    if (add9(mul9(m00, m11), neg9(mul9(m01, m10))) == 0) continue;
    int ok = 1, pm[4], om[4][9];
    for (int j = 0; j < 4 && ok; j++) for (int v = 0; v < 9 && ok; v++) {
      int px[9], py[9], np = 0;
      for (int x = 0; x < 9; x++) for (int y = 0; y < 9; y++) if (lform(j, x, y) == v) {
        int sx = x, sy = y; if (e) { sx = mul9(mul9(sx, sx), sx); sy = mul9(mul9(sy, sy), sy); }   /* Frobenius a -> a^3 */
        px[np] = add9(mul9(m00, sx), mul9(m01, sy)); py[np] = add9(mul9(m10, sx), mul9(m11, sy)); np++; }
      int c; for (c = 0; c < 4; c++) { int w = lform(c, px[0], py[0]), good = 1; for (int q = 1; q < np; q++) if (lform(c, px[q], py[q]) != w) good = 0; if (good) break; }
      if (c == 4) { ok = 0; break; }
      if (v == 0) pm[j] = c; else if (pm[j] != c) { ok = 0; break; }
      om[j][v] = lform(c, px[0], py[0]);
    }
    if (!ok) continue;
    int used = 0; for (int j = 0; j < 4; j++) used |= 1 << pm[j]; if (used != 15) continue;
    memcpy(perm[nL], pm, sizeof pm); memcpy(offmap[nL], om, sizeof om); nL++;
  }
  if (part == 0) fprintf(stderr, "semilinear direction-preserving maps: %d\n", nL);
  static unsigned short img[400][4][512], tr[4][81][512];
  for (int g = 0; g < nL; g++) for (int j = 0; j < 4; j++) for (int m = 0; m < 512; m++) {
    int r = 0; for (int v = 0; v < 9; v++) if (m >> v & 1) r |= 1 << offmap[g][j][v]; img[g][j][m] = r; }
  for (int j = 0; j < 4; j++) for (int T = 0; T < 81; T++) for (int m = 0; m < 512; m++) {
    int s = lform(j, T / 9, T % 9), r = 0; for (int v = 0; v < 9; v++) if (m >> v & 1) r |= 1 << add9(v, s); tr[j][T][m] = r; }
  /* fields: F_9 (table) with mu = image of i; generic F_{3^7}; others */
  static GF F9, F27, F81, F243, F729, F2187;
  gf_init(&F9, 2); gf_init(&F27, 3); gf_init(&F81, 4); gf_init(&F243, 5); gf_init(&F729, 6); gf_init(&F2187, 7);
  /* find i in F9 table: an element with i^2 = -1 */
  int iF9 = -1; for (int a = 1; a < 9; a++) if (gm(&F9, a, a) == F9.neg[1]) { iF9 = a; break; }
  GF *others[5] = {&F27, &F81, &F243, &F729, &F2187};
  static u128 m12[262144], m34[262144];
  for (int T1 = 0; T1 < 512; T1++) for (int T2 = 0; T2 < 512; T2++) {
    u128 m = 0, n = 0;
    for (int x = 0; x < 9; x++) for (int y = 0; y < 9; y++) { int P = 9 * x + y;
      if ((T1 >> y & 1) && (T2 >> x & 1)) m |= (u128)1 << P;
      if ((T1 >> lform(2, x, y) & 1) && (T2 >> lform(3, x, y) & 1)) n |= (u128)1 << P; }
    m12[T1 * 512 + T2] = m; m34[T1 * 512 + T2] = n; }
  static unsigned char co[300][40][2];
  long hits = 0, reps = 0, flex = 0, cand = 0;
  for (int A = part; A < 262144; A += nparts) {
    int T1 = A >> 9, T2 = A & 511; if (T1 == 511 || T2 == 511) continue;
    int tau12 = __builtin_popcount(T1) + __builtin_popcount(T2);
    for (int B = 0; B < 262144; B++) {
      int T3 = B >> 9, T4 = B & 511; if (T3 == 511 || T4 == 511) continue;
      int tau = tau12 + __builtin_popcount(T3) + __builtin_popcount(T4), s = 18 - tau; if (s < 2) continue;
      u128 w = m12[A] & m34[B]; int W = __builtin_popcountll((uint64_t)w) + __builtin_popcountll((uint64_t)(w >> 64));
      int ok = 0, uu = 0; for (int u = 1; 2 * u <= s; u++) if (u * (s - u) == W) { ok = 1; uu = u; }
      if (!ok) continue;
      hits++;
      int S[4] = {511 ^ T1, 511 ^ T2, 511 ^ T3, 511 ^ T4};
      uint64_t key = (uint64_t)S[0] << 27 | (uint64_t)S[1] << 18 | S[2] << 9 | S[3];
      int minimal = 1;
      for (int g = 0; g < nL && minimal; g++) {
        int I[4]; for (int j = 0; j < 4; j++) I[perm[g][j]] = img[g][j][S[j]];
        for (int T = 0; T < 81; T++) {
          uint64_t k2 = (uint64_t)tr[0][T][I[0]] << 27 | (uint64_t)tr[1][T][I[1]] << 18 | tr[2][T][I[2]] << 9 | tr[3][T][I[3]];
          if (k2 < key) { minimal = 0; break; } }
      }
      if (!minimal) continue;
      reps++;
      int col[4][9], K = 0;
      for (int j = 0; j < 4; j++) for (int v = 0; v < 9; v++) col[j][v] = (S[j] >> v & 1) ? K++ : -1;
      int E = 0;
      for (int x = 0; x < 9; x++) for (int y = 0; y < 9; y++) {
        int vv[4] = {y, x, lform(2, x, y), lform(3, x, y)}, in[4];
        for (int j = 0; j < 4; j++) in[j] = col[j][vv[j]];
        int cnt = (in[0] >= 0) + (in[1] >= 0) + (in[2] >= 0) + (in[3] >= 0); if (cnt < 3) continue;
        #define NEWROW() do { for (int q = 0; q < K; q++) co[E][q][0] = co[E][q][1] = 0; } while (0)
        /* coefficients mod 3 as (c0, c1): value c0 + c1*mu */
        if (in[0] >= 0 && in[1] >= 0) {
          if (in[2] >= 0) { NEWROW(); co[E][in[2]][0] = 1; co[E][in[1]][0] = 2; co[E][in[0]][0] = 1; E++; }          /* c - b + a */
          if (in[3] >= 0) { NEWROW(); co[E][in[3]][0] = 1; co[E][in[1]][0] = 2; co[E][in[0]][1] = 1; E++; }          /* d - b + mu a */
        } else if (in[0] >= 0) { NEWROW(); co[E][in[3]][0] = 1; co[E][in[2]][0] = 2; co[E][in[0]][0] = 2; co[E][in[0]][1] = 1; E++; } /* d - c - (1-mu) a */
        else { NEWROW(); co[E][in[3]][0] = 1; co[E][in[1]][0] = 2; co[E][in[1]][1] = 1; co[E][in[2]][1] = 2; E++; }        /* d - (1-mu) b - mu c */
      }
      /* the hit's own offsets live in our F_9 encoding; ranks are basis-independent, so evaluate at mu = i in the table field */
      int rt = rankF(&F9, E, K, co, iF9);
      if (K - rt != 3) fprintf(stderr, "WARN dim_t %d for %012llx\n", K - rt, (unsigned long long)key);
      int found = 0, isflex = 0;
      for (int fi = 0; fi < 5; fi++) { GF *F = others[fi];
        for (int mu = 2; mu < F->q; mu++) {
          /* degree over F_3 */
          int d = 1, z = mu; for (;;) { z = gm(F, gm(F, z, z), z); if (z == mu) break; d++; }
          if (d < F->k) continue;                              /* subfield elements are tested in their own field */
          int smallest = 1; z = mu; for (int i2 = 1; i2 < d; i2++) { z = gm(F, gm(F, z, z), z); if (z < mu) smallest = 0; }
          if (!smallest) continue;
          int r = rankF(F, E, K, co, mu);
          if (K - r >= 3) { found = 1; if (F->k == 7) isflex = 1;
            if (F->k < 7) printf("CAND %012llx K %d exps %d %d field 3^%d mu %d dim %d\n", (unsigned long long)key, K, 9 + uu, 9 + s - uu, F->k, mu, K - r); }
          if (F->k == 7) break;                                /* one generic element suffices */
        } }
      if (isflex) { flex++; printf("FLEX %012llx K %d exps %d %d\n", (unsigned long long)key, K, 9 + uu, 9 + s - uu); }
      if (found) cand++;
    }
  }
  fprintf(stderr, "part %d hits %ld orbit reps %ld flexible %ld with other-mu realisations %ld\n", part, hits, reps, flex, cand);
  return 0;
}
