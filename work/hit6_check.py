import sys, random
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice, multi_freeness, module_generators
p = 1009
m = [1, 1, 3, 2, 2, 4]
def B(t):
    return np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,1],[(1+t) % p, t % p, t % p],[0,1,1]], dtype=np.int64) % p
sig = Lattice(B(500), p).signature()
for t in [500, 77, 3, 815]:
    st, g = multi_freeness(B(t), m, p, rng=random.Random(t), tries=8)
    gens, _ = module_generators(B(t), m, p, kmax=13, stop_at_rank=False)
    print('t', t, 'same lattice', Lattice(B(t), p).signature() == sig, '| status', st, g, '| all min gens up to deg 13:', gens)
bad = []
for t in range(1, p):
    try:
        L = Lattice(B(t), p)
    except ValueError:
        continue
    if L.signature() != sig: continue
    st, g = multi_freeness(B(t), m, p, rng=random.Random(t), tries=4)
    if st != 'free': bad.append((t, st))
print('nonfree t over F_1009:', bad)
