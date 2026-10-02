import sys, itertools, time
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice, freeness, _cross3
p = int(sys.argv[1])
pairs = list(itertools.combinations(range(6), 2))
def arr(d):
    den = (2 * d - 1) % p
    Q = [(1,0,0),(0,1,0),(0,0,1),(1,1,1),(d*d % p, d % p, den),(d*den % p, d*d % p, den)]
    Q = [np.array(q, dtype=np.int64) % p for q in Q]
    return np.array([_cross3(Q[i], Q[j], p) for i, j in pairs], dtype=np.int64) % p
sig = None
for d0 in [5, 7, 11]:
    try:
        L0 = Lattice(arr(d0), p); sig = L0.signature(); break
    except ValueError:
        pass
print('generic t-vector', sorted(len(S) for S in L0.flats[2] if len(S) > 2), 'exps', L0.exponents())
same, bad, other = 0, [], {}
t0 = time.time()
for d in range(p):
    try:
        L = Lattice(arr(d), p)
    except ValueError:
        other['degenerate'] = other.get('degenerate', 0) + 1
        continue
    if L.signature() != sig:
        key = tuple(sorted(len(S) for S in L.flats[2] if len(S) > 2))
        other[key] = other.get(key, 0) + 1
        continue
    same += 1
    r = freeness(arr(d), p, L)
    if r['status'] != 'free':
        bad.append((d, r['status'], r['hilbert']))
print('p', p, 'same-lattice d:', same, 'NON-FREE with same lattice:', bad[:10], len(bad), 'time %.0f' % (time.time() - t0))
print('other lattices:', other)
