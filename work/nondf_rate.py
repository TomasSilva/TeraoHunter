import sys, random, itertools, time
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice
from hunter.ffree import free_rank3_ziegler
from hunter.divisional import rank3_df
q = int(sys.argv[1]); T = float(sys.argv[2])
pts = [v for v in itertools.product(range(q), repeat=3) if any(v) and next(c for c in v if c) == 1]
rng = random.Random(q)
st = {'tested': 0, 'chi': 0, 'free': 0, 'nondf_free': 0}
t0 = time.time()
found = []
while time.time() - t0 < T:
    n = rng.randint(15, min(len(pts) - 2, 45))
    S = rng.sample(range(len(pts)), n)
    N = np.array([pts[i] for i in S], dtype=np.int64)
    try:
        L = Lattice(N, q)
    except ValueError:
        continue
    st['tested'] += 1
    ex = L.exponents()
    if ex is None: continue
    st['chi'] += 1
    if free_rank3_ziegler(N, q, L) != 'free': continue
    st['free'] += 1
    if not rank3_df(N, q, L):
        st['nondf_free'] += 1
        found.append((n, ex, sorted(S)))
print(q, st, [(f[0], f[1]) for f in found[:10]])
