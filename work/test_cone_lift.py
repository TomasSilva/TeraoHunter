# Validate: A = cone(C) + {w} + translates {alpha_L - c w} (generic c) has
#   A free  <=>  (C, m) free,  m = 1 + #translates,  when A is locally free along w=0.
import sys, random, itertools
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice, freeness, local_freeness_l4, multi_freeness
p = 2147483629
rng = random.Random(3)
def lift(C, m):
    planes = [list(a) + [0] for a in C] + [[0, 0, 0, 1]]
    for a, k in zip(C, m):
        for _ in range(k - 1):
            planes.append(list(a) + [rng.randrange(1, 10**6)])
    return np.array(planes, dtype=np.int64) % p
tests = []
A3 = [[1,0,0],[0,1,0],[0,0,1],[1,-1,0],[1,0,-1],[0,1,-1]]
B3 = [[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,-1,0],[1,0,1],[1,0,-1],[0,1,1],[0,1,-1]]
for C in (A3, B3):
    n = len(C)
    L = Lattice(np.array(C) % p, p)
    pts = [S for S in L.flats[2] if len(S) >= 3]
    for trial in range(12):
        m = [1] * n
        # choose an 'independent' set J: no two J-lines share a multiple point
        order = list(range(n)); rng.shuffle(order)
        J = []
        for i in order:
            if all(not (i in S and j in S) for S in pts for j in J) and rng.random() < 0.6:
                J.append(i)
        for i in J: m[i] = rng.randint(2, 3)
        A = lift(C, m)
        LA = Lattice(A, p)
        loc = local_freeness_l4(A, p, LA)
        rA = freeness(A, p, LA)
        sm, g = multi_freeness(np.array(C) % p, m, p)
        print(len(C), 'm', m, '| (C,m):', sm, g, '| A:', rA['status'], rA['exps'], '| locally non-free flats:', len(loc))
