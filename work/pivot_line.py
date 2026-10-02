# B_s = B' + L0(s), L0(s) a line through a point X1 of B' with slope s, multiplicity m0.
# Look for (B', m', X1, m0) with (B_s, m) free for generic s; then scan s for non-free.
import sys, itertools, random, time
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import multi_freeness, Lattice
p = 32003
rng = random.Random(2)

def cross(u, w):
    return [u[1]*w[2]-u[2]*w[1], u[2]*w[0]-u[0]*w[2], u[0]*w[1]-u[1]*w[0]]

cores = {
 'A3': [[1,0,0],[0,1,0],[0,0,1],[1,-1,0],[1,0,-1],[0,1,-1]],
 'delA3': [[1,0,0],[0,1,0],[0,0,1],[1,-1,0],[1,0,-1]],
 'X3_2': [[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,2]],
 'bool': [[1,0,0],[0,1,0],[0,0,1]],
 'gen4': [[1,0,0],[0,1,0],[0,0,1],[1,1,1]],
 'B3sub7': [[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,-1,0],[1,0,1],[1,0,-1]],
}

def points_of(core):
    L = Lattice(np.array(core) % p, p)
    return [sorted(S) for S in L.flats[2]]

found = []
t0 = time.time()
for name, core in cores.items():
    pts = points_of(core)
    nc = len(core)
    for S in pts:
        u, w = core[S[0]], core[S[1]]
        def B(s):
            n0 = [(u[i] + s * w[i]) % p for i in range(3)]
            return np.array(core + [n0], dtype=np.int64) % p
        # generic lattice
        for mc in itertools.product(range(1, 4), repeat=nc):
            if sum(mc) > 14: continue
            M = sum(mc)
            for m0 in [M - 1, M, M + 2]:
                if m0 < 1: continue
                m = list(mc) + [m0]
                res = []
                for _ in range(2):
                    s = rng.randrange(2, p)
                    st, g = multi_freeness(B(s), m, p)
                    res.append((st, tuple(g)))
                if any(r[0] == 'free' for r in res):
                    found.append((name, S, mc, m0, res))
                    print(name, 'X1=', S, 'm\'=', mc, 'm0=', m0, res, flush=True)
print('done', len(found), time.time() - t0)
