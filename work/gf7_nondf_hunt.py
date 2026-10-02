import sys, random, itertools, time, json
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice
from hunter.ffree import free_rank3_ziegler
from hunter.divisional import rank3_df
from hunter.realize import Realizer
q = int(sys.argv[1]); T = float(sys.argv[2]); seed = int(sys.argv[3])
nlo, nhi = int(sys.argv[4]), int(sys.argv[5])
pts = [v for v in itertools.product(range(q), repeat=3) if any(v) and next(c for c in v if c) == 1]
rng = random.Random(seed)
t0 = time.time(); st = {'tested': 0, 'free': 0, 'nondf_free': 0, 'mixed': 0}
seen = set()
out = open(f'/Users/agent/Projects/TeraoHunter/work/gf{q}_nondf_{seed}.jsonl', 'w')
while time.time() - t0 < T:
    n = rng.randint(nlo, nhi)
    S = rng.sample(range(len(pts)), n)
    N = np.array([pts[i] for i in S], dtype=np.int64)
    try:
        L = Lattice(N, q)
    except ValueError:
        continue
    st['tested'] += 1
    if L.exponents() is None: continue
    if free_rank3_ziegler(N, q, L) != 'free': continue
    st['free'] += 1
    if rank3_df(N, q, L): continue
    st['nondf_free'] += 1
    key = (n, tuple(sorted(len(F) for F in L.flats[2])))
    if key in seen: continue
    seen.add(key)
    try:
        sols = Realizer(N, q, max_solutions=20000).solutions()
    except ValueError as e:
        continue
    fr = [free_rank3_ziegler(s, q) for s in sols]
    nf = len(fr) - fr.count('free')
    rec = {'n': n, 'exps': L.exponents(), 'S': sorted(S), 'realizations': len(sols), 'nonfree': nf}
    out.write(json.dumps(rec) + '\n'); out.flush()
    if nf and nf < len(fr):
        st['mixed'] += 1
        print('!!! MIXED FREENESS (TC COUNTEREXAMPLE OVER GF(%d)):' % q, json.dumps(rec), flush=True)
print(q, st, time.time() - t0, flush=True)
