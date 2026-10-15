/* pencils5x.c -- complete census over F_p of 17-line arrangements consisting of a line H = {z=0}
 * with five points and the five pencils through them.
 *   pencil i has direction (dp_i : dq_i : 0) on H and lines  dq_i x - dp_i y = v z  (offset v in F_p)
 *   directions: P0 = (1:0), P1 = (0:1), P2 = (1:1), P3 = (l1:1), P4 = (l2:1)
 *   sizes s0..s4 (the multiplicities m_P - 1 of the five points of H), s0+...+s4 = 16
 * Pencils 0..3 are enumerated (normalised: pencil 0 contains offsets 0 and 1, pencil 1 contains 0);
 * pencil 4 is derived: a line of pencil 4 meets the other 16 - s4 lines, so it has
 * 17 - s4 - ex points where ex = sum over its affine points R of (m_R - 2); all of ex comes from the
 * enumerated pencils. Window: every line has MINPTS..MAXPTS points, b2 = 80.
 * Since the five points of H are fixed to P0,P1,P2 and (l1,l2) runs over all admissible ordered
 * pairs, every configuration over F_p with these sizes is found (up to affine normalisation).
 *
 * usage: pencils5x s0 s1 s2 s3 s4 p l1 [l2]     (no l2: loop over all l2)
 * env: MINPTS (default 4), MAXPTS (default 8)
 * output: FOUND5 lines; stderr: summary.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int S[5], p, MINPTS = 4, MAXPTS = 8;
static int dp[5] = {1, 0, 1, 0, 0}, dq[5] = {0, 1, 1, 1, 1};
static int md(long long a) { a %= p; if (a < 0) a += p; return (int)a; }
static int inv(int a) { int r = 1, b = md(a), e = p - 2; while (e) { if (e & 1) r = (int)((long long)r * b % p); b = (int)((long long)b * b % p); e >>= 1; } return r; }

/* key tables: point of (i,u) and (j,v), i<j, as key = x*p + y */
static int *keytab[5][5];
static int *dvalOfKey;           /* offset of the pencil-4 line through an affine point */
static void buildTables(void) {
    for (int i = 0; i < 5; i++) for (int j = i + 1; j < 5; j++) {
        int det = md(-(long long)dq[i] * dp[j] + (long long)dp[i] * dq[j]);
        int di = inv(det);
        if (!keytab[i][j]) keytab[i][j] = malloc(sizeof(int) * p * p);
        for (int u = 0; u < p; u++) for (int v = 0; v < p; v++) {
            int x = md((long long)di * md(-(long long)dp[j] * u + (long long)dp[i] * v));
            int y = md((long long)di * md((long long)dq[i] * v - (long long)dq[j] * u));
            keytab[i][j][u * p + v] = x * p + y;
        }
    }
    if (!dvalOfKey) dvalOfKey = malloc(sizeof(int) * p * p);
    for (int x = 0; x < p; x++) for (int y = 0; y < p; y++)
        dvalOfKey[x * p + y] = md((long long)dq[4] * x - (long long)dp[4] * y);
}

/* boxed lines */
static int nb;                    /* number of enumerated lines = 16 - S[4] */
static int lpen[16], loff[16];    /* pencil and offset of each enumerated line */
static int lowEx[5], highEx[5];   /* excess bounds per pencil */
static unsigned *stampArr; static unsigned stamp = 0;
static uint32_t *maskArr;
static int touched[256];
static long long nEval = 0, nCombo = 0, nFound = 0;
static int curL1, curL2;

static int dEx[64], dNpts[64], dPts[64][32];
static int lineEx[16];

static void report(int *combo) {
    nFound++;
    printf("FOUND5 sizes %d %d %d %d %d p %d l1 %d l2 %d", S[0], S[1], S[2], S[3], S[4], p, curL1, curL2);
    const char *nm[5] = {"A", "B", "C", "E", "D"};
    int idx = 0;
    for (int i = 0; i < 4; i++) { printf(" %s", nm[i]); for (int t = 0; t < S[i]; t++) printf(" %d", loff[idx++]); }
    printf(" D"); for (int t = 0; t < S[4]; t++) printf(" %d", combo[t]);
    printf("\n"); fflush(stdout);
}

/* choose S[4] candidates with dEx summing to target, then verify the enumerated lines */
static int cand[64], ncand, combo[8];
static void finalCheck(void) {
    nCombo++;
    int fin[16];
    for (int a = 0; a < nb; a++) fin[a] = lineEx[a];
    for (int t = 0; t < S[4]; t++) {
        int d = combo[t];
        for (int q = 0; q < dNpts[d]; q++) {
            uint32_t m = maskArr[dPts[d][q]];
            while (m) { int a = __builtin_ctz(m); m &= m - 1; fin[a]++; }
        }
    }
    for (int a = 0; a < nb; a++) { int i = lpen[a]; if (fin[a] < lowEx[i] || fin[a] > highEx[i]) return; }
    report(combo);
}
static void chooseD(int start, int k, int sumLeft) {
    if (k == S[4]) { if (sumLeft == 0) finalCheck(); return; }
    for (int c = start; c <= ncand - (S[4] - k); c++) {
        int d = cand[c];
        if (dEx[d] > sumLeft) continue;
        combo[k] = d; chooseD(c + 1, k + 1, sumLeft - dEx[d]);
    }
}

static void evaluate(void) {
    nEval++;
    stamp++;
    int nt = 0;
    for (int a = 0; a < nb; a++) for (int b = a + 1; b < nb; b++) {
        int i = lpen[a], j = lpen[b];
        if (i == j) continue;
        int key = keytab[i][j][loff[a] * p + loff[b]];
        if (stampArr[key] != stamp) { stampArr[key] = stamp; maskArr[key] = 0; touched[nt++] = key; }
        maskArr[key] |= (1u << a) | (1u << b);
    }
    int b2box = 0;
    for (int a = 0; a < nb; a++) lineEx[a] = 0;
    for (int d = 0; d < p; d++) { dEx[d] = 0; dNpts[d] = 0; }
    for (int t = 0; t < nt; t++) {
        int key = touched[t]; uint32_t m = maskArr[key]; int k = __builtin_popcount(m);
        b2box += k - 1;
        if (k > 2) { uint32_t mm = m; while (mm) { int a = __builtin_ctz(mm); mm &= mm - 1; lineEx[a] += k - 2; } }
        int d = dvalOfKey[key]; dEx[d] += k - 1; dPts[d][dNpts[d]++] = key;
    }
    /* enumerated lines: excess only grows by at most S[4] */
    for (int a = 0; a < nb; a++) { int i = lpen[a]; if (lineEx[a] > highEx[i] || lineEx[a] + S[4] < lowEx[i]) return; }
    /* b2: final = 16 + b2box + sum over chosen d of (nb - dEx[d]) = 80 */
    int target = 16 + b2box + S[4] * nb - 80;
    ncand = 0;
    for (int d = 0; d < p; d++) if (dEx[d] >= lowEx[4] && dEx[d] <= highEx[4]) cand[ncand++] = d;
    if (ncand < S[4] || target < 0) return;
    chooseD(0, 0, target);
}

/* subset enumeration helpers */
static int boxStart[4];
static void enumPencil(int i, int pos, int lo) {
    /* fill pencil i from position pos with increasing offsets >= lo */
    int base = boxStart[i];
    if (pos == S[i]) { if (i == 3) evaluate(); else {
            int nextFixed = 0;
            enumPencil(i + 1, nextFixed, 0); }
        return; }
    for (int v = lo; v <= p - (S[i] - pos); v++) { lpen[base + pos] = i; loff[base + pos] = v; enumPencil(i, pos + 1, v + 1); }
}
static void enumAll(void) {
    /* pencil 0: offsets 0,1 fixed; pencil 1: offset 0 fixed */
    lpen[boxStart[0]] = 0; loff[boxStart[0]] = 0; lpen[boxStart[0] + 1] = 0; loff[boxStart[0] + 1] = 1;
    lpen[boxStart[1]] = 1; loff[boxStart[1]] = 0;
    /* custom recursion: pencil 0 from position 2 (lo 2), pencil 1 from position 1 (lo 1) */
    int s0 = S[0], s1 = S[1];
    int sel0[8], sel1[8];
    int k0 = s0 - 2, k1 = s1 - 1;
    for (int t = 0; t < k0; t++) sel0[t] = 2 + t;
    while (1) {
        for (int t = 0; t < k0; t++) { lpen[boxStart[0] + 2 + t] = 0; loff[boxStart[0] + 2 + t] = sel0[t]; }
        for (int t = 0; t < k1; t++) sel1[t] = 1 + t;
        while (1) {
            for (int t = 0; t < k1; t++) { lpen[boxStart[1] + 1 + t] = 1; loff[boxStart[1] + 1 + t] = sel1[t]; }
            enumPencil(2, 0, 0);
            int t = k1 - 1; while (t >= 0 && sel1[t] == p - k1 + t) t--; if (t < 0) break;
            sel1[t]++; for (int u = t + 1; u < k1; u++) sel1[u] = sel1[u - 1] + 1;
        }
        int t = k0 - 1; while (t >= 0 && sel0[t] == p - k0 + t) t--; if (t < 0) break;
        sel0[t]++; for (int u = t + 1; u < k0; u++) sel0[u] = sel0[u - 1] + 1;
    }
}

int main(int argc, char **argv) {
    if (argc < 8) { fprintf(stderr, "usage: pencils5x s0 s1 s2 s3 s4 p l1 [l2]\n"); return 1; }
    for (int i = 0; i < 5; i++) S[i] = atoi(argv[1 + i]);
    p = atoi(argv[6]);
    int l1 = md(atoi(argv[7]));
    int l2only = (argc >= 9) ? md(atoi(argv[8])) : -1;
    if (getenv("MINPTS")) MINPTS = atoi(getenv("MINPTS"));
    if (getenv("MAXPTS")) MAXPTS = atoi(getenv("MAXPTS"));
    if (S[0] + S[1] + S[2] + S[3] + S[4] != 16 || S[0] < 2 || S[1] < 1) { fprintf(stderr, "bad sizes\n"); return 1; }
    if (l1 == 0 || l1 == 1) { fprintf(stderr, "bad l1\n"); return 1; }
    nb = 16 - S[4];
    boxStart[0] = 0; boxStart[1] = S[0]; boxStart[2] = S[0] + S[1]; boxStart[3] = S[0] + S[1] + S[2];
    for (int i = 0; i < 5; i++) { lowEx[i] = 17 - S[i] - MAXPTS; highEx[i] = 17 - S[i] - MINPTS; if (lowEx[i] < 0) lowEx[i] = 0; }
    stampArr = calloc(p * p, sizeof(unsigned)); maskArr = calloc(p * p, sizeof(uint32_t));
    curL1 = l1; dp[3] = l1;
    for (int l2 = 2; l2 < p; l2++) {
        if (l2only >= 0 && l2 != l2only) continue;
        if (l2 == l1) continue;
        if (S[3] == S[4] && l2 < l1) continue;   /* pencils 3 and 4 interchangeable */
        curL2 = l2; dp[4] = l2;
        buildTables();
        long long e0 = nEval, f0 = nFound;
        enumAll();
        fprintf(stderr, "sizes %d %d %d %d %d p=%d l=(%d,%d): %lld evaluations, %lld found\n", S[0], S[1], S[2], S[3], S[4], p, l1, l2, nEval - e0, nFound - f0);
    }
    fprintf(stderr, "TOTAL sizes %d %d %d %d %d p=%d l1=%d: %lld evaluations, %lld combos, %lld found\n", S[0], S[1], S[2], S[3], S[4], p, l1, nEval, nCombo, nFound);
    return 0;
}
