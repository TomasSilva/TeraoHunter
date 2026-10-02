# Does adding a generic line L (multiplicity mu) to a core multiarrangement give
# free (B, m) for generic L?  If yes, scan positions of L for non-free special ones.
import sys, itertools, random, time
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import multi_freeness, Lattice
p = 32003
rng = random.Random(1)
cores = {
 'bool': [[1,0,0],[0,1,0],[0,0,1]],
 'A3': [[1,0,0],[0,1,0],[0,0,1],[1,-1,0],[1,0,-1],[0,1,-1]],
 'delA3': [[1,0,0],[0,1,0],[0,0,1],[1,-1,0],[1,0,-1]],
 'X3_2': [[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,2]],
 'nearpencil4': [[1,0,0],[0,1,0],[1,1,0],[0,0,1]],
 'B3sub': [[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,-1,0],[1,0,1],[1,0,-1]],
}
found = []
for name, core in cores.items():
    nc = len(core)
    for mcore in itertools.product(range(1, 4), repeat=nc):
        if nc >= 6 and sum(mcore) > 12: continue
        for mu in range(1, 7):
            m = list(mcore) + [mu]
            frees = []
            for trial in range(2):
                L = [rng.randrange(1, p), rng.randrange(1, p), 1]
                N = np.array(core + [L], dtype=np.int64) % p
                st, g = multi_freeness(N, m, p)
                frees.append((st, tuple(g)))
            if any(s == 'free' for s, _ in frees):
                found.append((name, tuple(mcore), mu, frees))
                print(name, mcore, mu, frees, flush=True)
print('done', len(found))
