# P^4 randomized: points in PG(4,q) built in special position (random / on a spanned line /
# in a spanned plane), A = hyperplanes spanned by >= 4 points spanning a solid.
# Filters: chi splits -> not DF (rank 5) -> free (exact).
import sys, random, time, json, itertools
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import rank_mod, nullspace_mod, Lattice
from hunter.lattice_gen import LatticeGen
from hunter.restrict import ziegler_restriction
from hunter.divisional import df_witnesses, poly_div
from hunter.ffree import ExactModule
q = int(sys.argv[1]); T = float(sys.argv[2]); seed = int(sys.argv[3]); KLO, KHI = int(sys.argv[4]), int(sys.argv[5])
rng = random.Random(seed)
def norm(v):
    v = [int(x) % q for x in v]
    for c in v:
        if c:
            inv = pow(c, q - 2, q); return tuple(x * inv % q for x in v)
def build(k):
    P = [tuple(int(i == j) for j in range(5)) for i in range(5)] + [tuple([1] * 5)]
    while len(P) < k:
        mode = rng.random()
        if mode < 0.4:      # on a line through two points
            a, b = rng.sample(P, 2); t = rng.randrange(1, q)
            v = norm([(a[i] + t * b[i]) % q for i in range(5)])
        elif mode < 0.8:    # in a plane through three points
            a, b, c = rng.sample(P, 3); s, t = rng.randrange(1, q), rng.randrange(1, q)
            v = norm([(a[i] + s * b[i] + t * c[i]) % q for i in range(5)])
        else:
            v = norm([rng.randrange(q) for _ in range(5)])
        if v and v not in P: P.append(v)
    return P
def spanned(P):
    H = {}
    for quad in itertools.combinations(range(len(P)), 4):
        M = np.array([P[i] for i in quad], dtype=np.int64) % q
        if rank_mod(M, q) < 4: continue
        v = norm(nullspace_mod(M, q)[0]); H.setdefault(v, set()).update(quad)
    return np.array(sorted(H), dtype=np.int64)
def df_rank5(N, chi):
    for i in range(N.shape[0]):
        R, mu, g = ziegler_restriction(N, i, q)
        if R.shape[0] < 4 or rank_mod(R, q) < 4: continue
        try:
            LR = Lattice(R, q)
        except ValueError:
            continue
        if poly_div(chi, LR.char_poly()) is None: continue
        if df_witnesses(R, q, LR): return True
    return False
st = {'built': 0, 'chi': 0, 'nondf': 0, 'free': 0}
t0 = time.time()
while time.time() - t0 < T:
    k = rng.randint(KLO, KHI)
    P = build(k)
    N = spanned(P)
    if N.shape[0] < 7 or rank_mod(N, q) < 5: continue
    st['built'] += 1
    try:
        L = LatticeGen(N, q)
    except ValueError:
        continue
    ex = L.exponents()
    if ex is None: continue
    st['chi'] += 1
    if df_rank5(N, L.char_poly()): continue
    st['nondf'] += 1
    fr, g = ExactModule(N, [1] * N.shape[0], q).freeness()
    if fr == 'free':
        st['free'] += 1
        print('FREE NON-DF P^4', json.dumps({'q': q, 'P': [list(x) for x in P], 'n': int(N.shape[0]), 'exps': ex}), flush=True)
print(q, st, time.time() - t0, flush=True)
