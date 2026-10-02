# Analyse the free-but-not-DF arrangement over F_3 found by the exhaustive PG(3,3) run.
import sys, json, itertools
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice, rank_mod, integer_roots
from hunter.restrict import ziegler_restriction
from hunter.ffree import ExactModule
from hunter.m2 import m2_free

q = 3
rec = json.loads(open('/Users/agent/.claude/jobs/ef50d01d/tmp/hit.json').read())
N = np.array(rec['normals']) % q
n = len(N)
L = Lattice(N, q)
print('n =', n, 'chi =', L.char_poly(), 'exps =', integer_roots(L.char_poly()))
# restrictions
sizes = []
for i in range(n):
    R, mu, g = ziegler_restriction(N, i, q)
    sizes.append(R.shape[0])
print('|A^H| for H in A:', sizes)
# points of P^3 (rank-3 flats) by multiplicity, lines (rank-2 flats) by multiplicity
from collections import Counter
print('rank-2 flats by size:', sorted(Counter(len(S) for S in L.flats[2]).items()))
print('rank-3 flats by size:', sorted(Counter(len(S) for S in L.flats[3]).items()))
# exact freeness again + M2 over GF(3)
print('ExactModule:', ExactModule(N, [1] * n, q).freeness())
res = m2_free(N.tolist(), field='ZZ/3')
print('M2 over GF(3):', res.get('free'), res.get('gens'), res.get('betti'))
# complement
allpts = [v for v in itertools.product(range(q), repeat=4) if any(v) and next(c for c in v if c) == 1]
S = set(map(tuple, N.tolist()))
comp = [v for v in allpts if v not in S]
print('complement (%d):' % len(comp), comp)
Lc = Lattice(np.array(comp), q)
print('complement chi:', Lc.char_poly())
