import sys, json, itertools
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice
from hunter.realize import AbstractRealizer
from hunter.ffree import free_rank3_ziegler
rec = [json.loads(l) for l in open('/Users/agent/Projects/TeraoHunter/work/connecting_11.jsonl') if json.loads(l)['n_lines'] == 13][0]
q0 = 11
pts = [v for v in itertools.product(range(q0), repeat=3) if any(v) and next(c for c in v if c) == 1]
P = set(rec['P'])
inc = [[j for j, X in enumerate(pts) if sum(a*b for a, b in zip(L, X)) % q0 == 0] for L in pts]
Ls = [i for i in range(len(pts)) if sum(1 for j in inc[i] if j in P) >= 2]
N = np.array([pts[i] for i in Ls], dtype=np.int64)
L = Lattice(N, q0)
flats = [sorted(S) for S in L.flats[2]]
json.dump({'n': len(Ls), 'flats2': flats}, open('/Users/agent/Projects/TeraoHunter/work/lattice13.json', 'w'))
print('13-line lattice: points', sorted(len(S) for S in flats if len(S) > 2), 'exps', L.exponents())
for p in [101, 1009]:
    sols = AbstractRealizer(len(Ls), flats, p).solutions()
    fr = [free_rank3_ziegler(s, p) for s in sols]
    print('p', p, 'realizations', len(sols), 'free', fr.count('free'))
    np.save(f'/Users/agent/Projects/TeraoHunter/work/real13_{p}.npy', np.array(sols))
