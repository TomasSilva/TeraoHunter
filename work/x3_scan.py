import sys, itertools, time
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import multi_freeness, Lattice
p = 32003
def X3(t):
    return np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,t % p]], dtype=np.int64) % p

ts = [1, 2, 3, 5, 7, 11, 12345, (3 * pow(2, p-2, p)) % p, p - 2, (p+1)//2]  # include t=-2, t=1/2
sig0 = Lattice(X3(12345), p).signature()
good_ts = []
for t in ts:
    s = Lattice(X3(t), p).signature()
    print('t', t, 'same lattice' if s == sig0 else 'DIFFERENT lattice', [sorted(f) for f in s[0] if len(f) > 2])
    if s == sig0: good_ts.append(t)
t0 = time.time()
diff = []
for m in itertools.product(range(1, 5), repeat=6):
    res = []
    for t in good_ts:
        st, gens = multi_freeness(X3(t), m, p)
        res.append((st == 'free', tuple(gens)))
    frees = set(r[0] for r in res)
    if len(frees) > 1:
        diff.append((m, list(zip(good_ts, res))))
        print('FREENESS DEPENDS ON t:', m, list(zip(good_ts, res)))
print('done', time.time() - t0, 'n diff', len(diff))
