# P^3 randomized point configurations in special position; A = planes spanned by >= 3 points.
# Filters: chi splits -> locally free -> not DF -> free -> realization enumeration (mixed?).
import sys, random, time, json, itertools
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import rank_mod, nullspace_mod, Lattice, localization_l4
from hunter.ffree import free_rank3_ziegler, free_simple
from hunter.divisional import df_witnesses
from hunter.realize import Realizer
q = int(sys.argv[1]); T = float(sys.argv[2]); seed = int(sys.argv[3]); KLO, KHI = int(sys.argv[4]), int(sys.argv[5])
rng = random.Random(seed)
def norm(v):
    v = [int(x) % q for x in v]
    for c in v:
        if c:
            inv = pow(c, q - 2, q); return tuple(x * inv % q for x in v)
def build(k):
    P = [tuple(int(i == j) for j in range(4)) for i in range(4)] + [tuple([1] * 4)]
    while len(P) < k:
        r = rng.random()
        if r < 0.45:
            a, b = rng.sample(P, 2); t = rng.randrange(1, q)
            v = norm([(a[i] + t * b[i]) % q for i in range(4)])
        elif r < 0.75:
            # intersection of a spanned line with a spanned plane
            a, b = rng.sample(P, 2); c, d, e = rng.sample(P, 3)
            pl = nullspace_mod(np.array([c, d, e], dtype=np.int64) % q, q)
            if pl.shape[0] != 1: continue
            fa, fb = int(np.dot(pl[0], a) % q), int(np.dot(pl[0], b) % q)
            v = norm([(fb * a[i] - fa * b[i]) % q for i in range(4)])
        else:
            a, b, c = rng.sample(P, 3); s, t = rng.randrange(1, q), rng.randrange(1, q)
            v = norm([(a[i] + s * b[i] + t * c[i]) % q for i in range(4)])
        if v and v not in P: P.append(v)
    return P
def spanned(P):
    H = {}
    for tri in itertools.combinations(range(len(P)), 3):
        M = np.array([P[i] for i in tri], dtype=np.int64) % q
        if rank_mod(M, q) < 3: continue
        v = norm(nullspace_mod(M, q)[0]); H.setdefault(v, set()).update(tri)
    return np.array(sorted(H), dtype=np.int64)
def locally_free(N, L):
    for S in L.flats[3]:
        if len(S) > 3 and free_rank3_ziegler(localization_l4(N, S, q), q) != 'free': return False
    return True
st = {'built': 0, 'chi': 0, 'locfree': 0, 'nondf': 0, 'free': 0, 'mixed': 0}
t0 = time.time(); seen = set()
while time.time() - t0 < T:
    P = build(rng.randint(KLO, KHI))
    N = spanned(P)
    if N.shape[0] < 6 or rank_mod(N, q) < 4: continue
    st['built'] += 1
    try:
        L = Lattice(N, q)
    except ValueError:
        continue
    if L.exponents() is None: continue
    st['chi'] += 1
    if not locally_free(N, L): continue
    st['locfree'] += 1
    if df_witnesses(N, q, L): continue
    st['nondf'] += 1
    fr, g = free_simple(N, q, L)
    if fr != 'free': continue
    st['free'] += 1
    key = (len(N), tuple(sorted(len(S) for S in L.flats[2])), tuple(sorted(len(S) for S in L.flats[3])))
    if key in seen: continue
    seen.add(key)
    try:
        sols = Realizer(N, q, max_solutions=3000).solutions()
    except ValueError:
        sols = []
    frs = [free_simple(s, q)[0] for s in sols]
    nf = len(frs) - frs.count('free')
    print('FREE NON-DF', json.dumps({'q': q, 'n': int(len(N)), 'exps': L.exponents(), 'realizations': len(sols), 'nonfree': nf, 'P': [list(x) for x in P]}), flush=True)
    if nf and nf < len(frs):
        st['mixed'] += 1; print('!!! MIXED (TC COUNTEREXAMPLE OVER GF(%d))' % q, flush=True)
print(q, st, time.time() - t0, flush=True)
