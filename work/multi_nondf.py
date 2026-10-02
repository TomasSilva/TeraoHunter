# Multiplicities on the non-DF free families with moduli (13-line, 15-line), over F_p.
# For random m: generic freeness (2 realizations); if free, test all realizations.
import sys, random, time, json, itertools
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice, multi_freeness, _cross3
from hunter.multi import local_exponents, wakamiko
which = sys.argv[1]; p = int(sys.argv[2]); T = float(sys.argv[3]); seed = int(sys.argv[4])
rng = random.Random(seed)
if which == '15':
    pairs = list(itertools.combinations(range(6), 2))
    def arr(d):
        den = (2 * d - 1) % p
        Q = [(1,0,0),(0,1,0),(0,0,1),(1,1,1),(d*d % p, d % p, den),(d*den % p, d*d % p, den)]
        Q = [np.array(x, dtype=np.int64) % p for x in Q]
        return np.array([_cross3(Q[i], Q[j], p) for i, j in pairs], dtype=np.int64) % p
    sig = Lattice(arr(5), p).signature()
    reals = []
    for d in range(2, p):
        try:
            if Lattice(arr(d), p).signature() == sig: reals.append(arr(d))
        except ValueError:
            pass
else:
    reals = list(np.load(f'/Users/agent/Projects/TeraoHunter/work/real13_{p}.npy'))
n = reals[0].shape[0]
L0 = Lattice(reals[0], p); pts = [sorted(S) for S in L0.flats[2]]
print(which, 'realizations', len(reals), 'lines', n, flush=True)
def b2(m, A):
    tot = 0
    for S in pts:
        if len(S) == 2: tot += m[S[0]] * m[S[1]]
        else:
            e = local_exponents(A, m, S, p); tot += e[0] * e[1]
    return tot
def exps_ok(M, b):
    for d1 in range(0, M // 3 + 1):
        for d2 in range(d1, (M - d1) // 2 + 1):
            d3 = M - d1 - d2
            if d1*d2 + d1*d3 + d2*d3 == b: return True
    return False
t0 = time.time(); tested = gfree = 0; seen = set()
while time.time() - t0 < T:
    m = [rng.choice([1, 1, 1, 2, 2, 3]) for _ in range(n)]
    if tuple(m) in seen: continue
    seen.add(tuple(m))
    A1, A2 = reals[rng.randrange(len(reals))], reals[rng.randrange(len(reals))]
    bb = b2(m, A1)
    if bb != b2(m, A2):
        continue   # local exponents vary: local phenomenon, skip
    if not exps_ok(sum(m), bb): continue
    tested += 1
    s1, g1 = multi_freeness(A1, m, p)
    if s1 != 'free': continue
    s2, g2 = multi_freeness(A2, m, p)
    if s2 != 'free': continue
    gfree += 1
    bad = []
    for i, A in enumerate(reals):
        st, g = multi_freeness(A, m, p)
        if st != 'free': bad.append(i)
    print('GENERIC-FREE', which, m, 'exps', g1, 'nonfree realizations', len(bad), 'of', len(reals), flush=True)
print('tested', tested, 'generic free', gfree, time.time() - t0, flush=True)
