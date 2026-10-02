import sys, itertools, random
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice, freeness, local_freeness_l4, _minors3, nullspace_mod
from hunter.divisional import df_witnesses
p = 32003
rng = random.Random(1)
def planes_from_points(P):
    out = []
    for a, b, c in itertools.combinations(range(len(P)), 3):
        K = nullspace_mod(np.array([P[a], P[b], P[c]], dtype=np.int64) % p, p)
        out.append(K[0])
    return np.array(out, dtype=np.int64) % p
for k in [5, 6, 7]:
    for trial in range(2):
        P = [[rng.randrange(p) for _ in range(4)] for _ in range(k)]
        N = planes_from_points(P)
        L = Lattice(N, p)
        loc = local_freeness_l4(N, p, L)
        fr = freeness(N, p, L)
        print(k, 'points ->', len(N), 'planes, chi', L.char_poly(), 'exps', L.exponents(), 'locally non-free', len(loc), fr['status'])
