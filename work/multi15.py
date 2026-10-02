# Multiplicities on the 15-line non-DF family C_d: generic-free (C_d, m) with d-dependent freeness?
import sys, itertools, time, random, json
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice, multi_freeness, _cross3
from hunter.multi import local_exponents
p = 1009
pairs = list(itertools.combinations(range(6), 2))
def arr(d):
    den = (2 * d - 1) % p
    Q = [(1,0,0),(0,1,0),(0,0,1),(1,1,1),(d*d % p, d % p, den),(d*den % p, d*d % p, den)]
    Q = [np.array(q, dtype=np.int64) % p for q in Q]
    return np.array([_cross3(Q[i], Q[j], p) for i, j in pairs], dtype=np.int64) % p
Lg = Lattice(arr(5), p); sig = Lg.signature()
pts = [sorted(S) for S in Lg.flats[2]]
good_d = [d for d in range(p) if (lambda: True)()]
same_d = []
for d in range(p):
    try:
        if Lattice(arr(d), p).signature() == sig: same_d.append(d)
    except ValueError:
        pass
print('same-lattice d count', len(same_d))
# local exponents at 5-fold points vary with d?  test for multiplicity patterns
def b2(m, d):
    tot = 0
    for S in pts:
        if len(S) == 2: tot += m[S[0]] * m[S[1]]
        else:
            e = local_exponents(arr(d), m, S, p); tot += e[0] * e[1]
    return tot
def exps_for(M, b2v):
    out = []
    for d1 in range(0, M // 3 + 1):
        for d2 in range(d1, (M - d1) // 2 + 1):
            d3 = M - d1 - d2
            if d1*d2 + d1*d3 + d2*d3 == b2v: out.append((d1, d2, d3))
    return out
rng = random.Random(int(sys.argv[1]))
t0 = time.time(); tested = 0
seen = set()
while time.time() - t0 < float(sys.argv[2]):
    k = rng.randint(1, 6)
    J = rng.sample(range(15), k)
    m = [1] * 15
    for j in J: m[j] = rng.choice([2, 2, 3])
    key = tuple(m)
    if key in seen: continue
    seen.add(key)
    d1, d2 = rng.choice(same_d), rng.choice(same_d)
    if b2(m, d1) != b2(m, d2):
        print('LOCAL-EXPONENT VARIATION', m, flush=True); continue
    cand = exps_for(sum(m), b2(m, d1))
    if not cand: continue
    tested += 1
    s1, g1 = multi_freeness(arr(d1), m, p)
    if s1 != 'free': continue
    s2, g2 = multi_freeness(arr(d2), m, p)
    if s2 != 'free': continue
    bad = []
    for d in same_d:
        st, g = multi_freeness(arr(d), m, p)
        if st != 'free': bad.append((d, g))
    print('GENERIC-FREE', m, 'exps', g1, 'nonfree d:', bad, flush=True)
print('tested', tested, 'time', time.time() - t0)
