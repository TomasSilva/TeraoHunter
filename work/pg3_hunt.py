# Rank-4 hunt over GF(q): random plane subsets of PG(3,q); filters
#   chi splits -> locally free (rank-3 localizations, exact Ziegler criterion)
#   -> NOT divisionally free -> free (exact) -> enumerate ALL GF(q)-realizations
#   of the lattice and test each: mixed freeness = Terao counterexample over GF(q).
import sys, random, itertools, time, json
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice, localization_l4
from hunter.ffree import free_rank3_ziegler, free_simple
from hunter.divisional import df_witnesses
from hunter.realize import Realizer

q = int(sys.argv[1]); T = float(sys.argv[2]); seed = int(sys.argv[3])
nlo, nhi = int(sys.argv[4]), int(sys.argv[5])
planes = [v for v in itertools.product(range(q), repeat=4) if any(v) and next(c for c in v if c) == 1]
rng = random.Random(seed)
out = open(f'/Users/agent/Projects/TeraoHunter/work/pg3_{q}_{seed}.jsonl', 'w')
st = {'tested': 0, 'chi': 0, 'locfree': 0, 'nondf': 0, 'free': 0, 'mixed': 0}
t0 = time.time()
seen = set()


def locally_free(N, L):
    for S in L.flats[3]:
        if len(S) <= 3:
            continue
        loc = localization_l4(N, S, q)
        if free_rank3_ziegler(loc, q) != 'free':
            return False
    return True


while time.time() - t0 < T:
    n = rng.randint(nlo, nhi)
    S = rng.sample(range(len(planes)), n)
    N = np.array([planes[i] for i in S], dtype=np.int64)
    if np.linalg.matrix_rank(N.astype(float)) < 4:
        continue
    try:
        L = Lattice(N, q)
    except ValueError:
        continue
    st['tested'] += 1
    if L.exponents() is None:
        continue
    st['chi'] += 1
    if not locally_free(N, L):
        continue
    st['locfree'] += 1
    if df_witnesses(N, q, L):
        continue
    st['nondf'] += 1
    fr, g = free_simple(N, q, L)
    rec = {'n': n, 'exps': L.exponents(), 'S': sorted(S), 'status': fr}
    if fr != 'free':
        out.write(json.dumps(rec) + '\n'); out.flush()
        continue
    st['free'] += 1
    key = (n, tuple(sorted(len(F) for F in L.flats[2])), tuple(sorted(len(F) for F in L.flats[3])))
    if key in seen:
        continue
    seen.add(key)
    try:
        sols = Realizer(N, q, max_solutions=5000).solutions()
    except ValueError:
        sols = []
    frs = [free_simple(s, q)[0] for s in sols]
    nf = len(frs) - frs.count('free')
    rec.update({'realizations': len(sols), 'nonfree': nf})
    out.write(json.dumps(rec) + '\n'); out.flush()
    print('FREE NON-DF', json.dumps(rec), flush=True)
    if nf and nf < len(frs):
        st['mixed'] += 1
        print('!!! MIXED FREENESS (TC COUNTEREXAMPLE OVER GF(%d))' % q, flush=True)
print(q, st, time.time() - t0, flush=True)
