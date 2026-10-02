import sys, itertools
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice, _cross3, nullspace_mod, rank_mod
from hunter.ffree import free_rank3_ziegler, ExactModule
pairs = list(itertools.combinations(range(6), 2))
def arr(d, p):
    den = (2 * d - 1) % p
    Q = [(1,0,0),(0,1,0),(0,0,1),(1,1,1),(d*d % p, d % p, den),(d*den % p, d*d % p, den)]
    Q = [np.array(x, dtype=np.int64) % p for x in Q]
    return np.array([_cross3(Q[i], Q[j], p) for i, j in pairs], dtype=np.int64) % p, Q
def on_conic(pts, p):
    M = np.array([[x*x, y*y, z*z, x*y, x*z, y*z] for (x, y, z) in pts], dtype=np.int64) % p
    return rank_mod(M, p) < 6
for p in [1009, 1013, 10009, 10037]:
    N0, _ = arr(5, p); sig = Lattice(N0, p).signature()
    hits = []
    for d in range(2, p):
        try:
            N, Q = arr(d, p); L = Lattice(N, p)
        except ValueError:
            continue
        same = L.signature() == sig
        trip = [ [int(c) for c in nullspace_mod(N[sorted(S)], p)[0]] for S in L.flats[2] if len(S) == 3]
        qc = on_conic([tuple(int(c) for c in q) for q in Q], p)
        tc = len(trip) == 6 and on_conic(trip, p)
        if qc or tc:
            hits.append((d, 'Q-on-conic' if qc else '', 'T-on-conic' if tc else '', 'same lattice' if same else 'LATTICE CHANGES',
                         free_rank3_ziegler(N, p, L) if same else '-'))
    print('p', p, hits)
