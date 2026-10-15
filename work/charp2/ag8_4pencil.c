/* Four-pencil sub-arrangements of PG(2,8) containing H = {z = 0}.
   Directions (1:0), (0:1), (1:1), (t:1) with F_8 = F_2[t]/(t^3+t+1); every 4-subset of P^1(F_8) is
   PGL(2,8)-equivalent to this one.  Lines: y = a (a in S1), x = b (S2), x+y = c (S3), x+t*y = d (S4).
   For A in PG(2,q) containing H with K = n-1 affine lines covering U of the q^2 affine points:
   chi_0(0) = qK - U, so chi splits iff W = q^2 - U = u*v with u + v = K - 2q, exponents (q+u, q+v).
   Print every (S1..S4) with all four directions used and u >= 1 (d1 >= 9 > 8 = q). */
#include <stdio.h>
#include <stdint.h>
static int mul8(int a, int b) { int r = 0; for (int i = 0; i < 3; i++) if (b >> i & 1) r ^= a << i;
  for (int i = 4; i >= 3; i--) if (r >> i & 1) r ^= 0xB << (i - 3); return r; }   /* t^3 = t + 1 */
int main(void) {
  static uint64_t m12[65536], m34[65536];
  int t = 2;
  for (int T1 = 0; T1 < 256; T1++) for (int T2 = 0; T2 < 256; T2++) {
    uint64_t m = 0, n = 0;
    for (int x = 0; x < 8; x++) for (int y = 0; y < 8; y++) {
      int P = 8 * x + y;
      if ((T1 >> y & 1) && (T2 >> x & 1)) m |= 1ULL << P;
      if ((T1 >> (x ^ y) & 1) && (T2 >> (x ^ mul8(t, y)) & 1)) n |= 1ULL << P;
    }
    m12[T1 * 256 + T2] = m; m34[T1 * 256 + T2] = n;
  }
  long hits = 0;
  for (int A = 0; A < 65536; A++) {
    int T1 = A >> 8, T2 = A & 255;
    if (T1 == 255 || T2 == 255) continue;
    int tau12 = __builtin_popcount(T1) + __builtin_popcount(T2);
    for (int B = 0; B < 65536; B++) {
      int T3 = B >> 8, T4 = B & 255;
      if (T3 == 255 || T4 == 255) continue;
      int tau = tau12 + __builtin_popcount(T3) + __builtin_popcount(T4);
      int s = 16 - tau;                       /* u + v */
      if (s < 2) continue;
      int W = __builtin_popcountll(m12[A] & m34[B]);
      for (int u = 1; 2 * u <= s; u++) if (u * (s - u) == W) {
        printf("%02x %02x %02x %02x K %d W %d exps %d %d\n", 255 ^ T1, 255 ^ T2, 255 ^ T3, 255 ^ T4, 32 - tau, W, 8 + u, 8 + s - u);
        hits++; break;
      }
    }
  }
  fprintf(stderr, "hits %ld\n", hits);
  return 0;
}
