import sys, random, itertools, time, json
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice, freeness, local_freeness_l4
from hunter.divisional import df_witnesses
p = 2147483629
V = []
for v in itertools.product([-1, 0, 1], repeat=4):
    if any(v) and next(c for c in v if c) == 1:
        V.append(list(v))
assert len(V) == 40
seed = int(sys.argv[1]); rng = random.Random(seed)
t0 = time.time(); stats = {'tested': 0, 'chi': 0, 'locfree': 0, 'free': 0, 'nondf': 0}
seen = set()
while time.time() - t0 < float(sys.argv[2]):
    k = rng.randint(10, 26)
    S = sorted(rng.sample(range(40), k))
    N = np.array([V[i] for i in S], dtype=np.int64) % p
    if np.linalg.matrix_rank(np.array([V[i] for i in S])) < 4: continue
    L = Lattice(N, p)
    stats['tested'] += 1
    if L.exponents() is None: continue
    stats['chi'] += 1
    if local_freeness_l4(N, p, L): continue
    stats['locfree'] += 1
    fr = freeness(N, p, L)
    if fr['status'] != 'free': continue
    stats['free'] += 1
    w = df_witnesses(N, p, L)
    if not w:
        stats['nondf'] += 1
        key = (len(S), tuple(sorted(len(F) for F in L.flats[2])), tuple(sorted(len(F) for F in L.flats[3])))
        if key not in seen:
            seen.add(key)
            print('NON-DF FREE', json.dumps({'S': S, 'exps': fr['exps']}), flush=True)
print(stats, time.time() - t0)
