import sys, pickle, random, itertools
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice
from hunter.ffree import free_rank3_ziegler
from hunter.m2 import run_m2
q = 5
levels = pickle.load(open('/Users/agent/Projects/TeraoHunter/work/pg2_5_orbits.pkl', 'rb'))
lines = [v for v in itertools.product(range(q), repeat=3) if any(v) and next(c for c in v if c) == 1]
rng = random.Random(9)
pool = [(k, S) for k, reps in levels.items() if 15 <= k <= 28 for S in reps]
sample = rng.sample(pool, 40)
bad = 0
for k, S in sample:
    N = np.array([lines[i] for i in S], dtype=np.int64)
    L = Lattice(N, q)
    ours = free_rank3_ziegler(N, q, L) == 'free'
    forms = ["+".join(f"({c})*{x}" for c, x in zip(v, 'xyz') if c) for v in N.tolist()]
    out = run_m2(f'needsPackage "HyperplaneArrangements"; R = ZZ/5[x,y,z]; A = arrangement {{{",".join(forms)}}}; print("FREE " | toString isFreeModule prune image der A);')
    theirs = 'FREE true' in out
    bad += ours != theirs
    print(k, L.exponents(), 'ours', ours, 'M2', theirs, '' if ours == theirs else 'MISMATCH', flush=True)
print('mismatches', bad, 'of', len(sample))
