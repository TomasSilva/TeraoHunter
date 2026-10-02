import sys, itertools, time
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice, _cross3
from hunter.realize import Realizer
from hunter.ffree import free_rank3_ziegler
pairs = list(itertools.combinations(range(6), 2))
def arr(d, p):
    den = (2 * d - 1) % p
    Q = [(1,0,0),(0,1,0),(0,0,1),(1,1,1),(d*d % p, d % p, den),(d*den % p, d*d % p, den)]
    Q = [np.array(x, dtype=np.int64) % p for x in Q]
    return np.array([_cross3(Q[i], Q[j], p) for i, j in pairs], dtype=np.int64) % p
target_t = sorted([3]*6 + [5]*6)
for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61]:
    ref = None
    for d in range(2, p):
        try:
            L = Lattice(arr(d, p), p)
        except ValueError:
            continue
        if sorted(len(S) for S in L.flats[2] if len(S) > 2) == target_t and L.exponents() == (1, 7, 7):
            ref = arr(d, p); break
    if ref is None:
        print('p', p, ': no realization of the generic lattice on the curve'); continue
    t = time.time()
    sols = Realizer(ref, p).solutions()
    fr = [free_rank3_ziegler(s, p) for s in sols]
    print('p %d: %d realizations up to PGL; free: %d, non-free: %d  (%.1fs)' % (p, len(sols), fr.count('free'), len(fr) - fr.count('free'), time.time() - t), flush=True)
