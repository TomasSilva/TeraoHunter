# Which multiplicities on X3(t) are free for generic t?  Then scan all t in F_p for non-free special t.
import sys, itertools, time
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import multi_freeness, Lattice
p = 1009
def X3(t):
    return np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,t % p]], dtype=np.int64) % p
sig0 = Lattice(X3(123), p).signature()
good_ts = [t for t in range(1, p) if Lattice(X3(t), p).signature() == sig0]
print('good t count', len(good_ts), 'bad t', [t for t in range(1,p) if t not in set(good_ts)])
t0 = time.time()
free_m = []
for m in itertools.product(range(1, 5), repeat=6):
    st, gens = multi_freeness(X3(123), m, p)
    st2, gens2 = multi_freeness(X3(777), m, p)
    if st == 'free' or st2 == 'free':
        free_m.append((m, gens, st, st2))
print('generically free count', len(free_m), time.time() - t0, flush=True)
for m, gens, st, st2 in free_m:
    print(m, gens, st, st2, flush=True)
