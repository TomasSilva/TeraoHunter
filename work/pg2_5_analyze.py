import sys, pickle, itertools
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
import pynauty
from hunter.core import Lattice
from hunter.ffree import free_rank3_ziegler
from hunter.divisional import rank3_df
q = 5
levels = pickle.load(open('/Users/agent/Projects/TeraoHunter/work/pg2_5_orbits.pkl', 'rb'))
pts = [v for v in itertools.product(range(q), repeat=3) if any(v) and next(c for c in v if c) == 1]
lines = pts
NP = len(pts)
inc = [set(j for j, P in enumerate(pts) if sum(a*b for a, b in zip(L, P)) % q == 0) for L in lines]
def lattice_cert(S):
    S = sorted(S); flats = []
    for j in range(NP):
        th = [i for i in S if j in inc[i]]
        if len(th) >= 2: flats.append(th)
    idx = {L: t for t, L in enumerate(S)}; nv = len(S) + len(flats)
    ad = {len(S) + f: [idx[L] for L in F] for f, F in enumerate(flats)}
    g = pynauty.Graph(nv, directed=False, adjacency_dict=ad, vertex_coloring=[set(range(len(S))), set(range(len(S), nv))])
    return pynauty.certificate(g)
classes = {}
nondf_free = []
for k, reps in levels.items():
    if k < 3: continue
    for S in reps:
        N = np.array([lines[i] for i in S], dtype=np.int64)
        if np.linalg.matrix_rank(N.astype(float)) < 3: continue
        L = Lattice(N % q, q)
        st = free_rank3_ziegler(N, q, L)
        df = rank3_df(N % q, q, L) if L.exponents() else False
        classes.setdefault((k, lattice_cert(S)), []).append((S, st, df, L.exponents()))
        if st == 'free' and not df:
            nondf_free.append((k, S, L.exponents()))
multi = [(k, v) for (k, c), v in classes.items() if len(v) > 1]
print('free non-DF orbits over GF(5):', len(nondf_free), sorted(set((k, e) for k, S, e in nondf_free)))
print('multi-orbit lattice classes:', len(multi))
for k, v in sorted(multi, key=lambda x: x[0]):
    print('  n=%d orbits=%d freeness=%s DF=%s exps=%s' % (k, len(v), [x[1] for x in v], [x[2] for x in v], v[0][3]))
