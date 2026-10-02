import sys, random, time, itertools
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice
from hunter.lattice_gen import LatticeGen
from hunter.ffree import ExactModule
p = 1000003
rng = random.Random(1)
# compare with old Lattice for l = 3, 4
for l in [3, 4]:
    for _ in range(5):
        S = set()
        while len(S) < rng.randint(8, 14):
            v = tuple(rng.randint(-1, 1) for _ in range(l))
            if any(v):
                f = next(c for c in v if c); S.add(tuple(c * f for c in v))
        N = np.array(list(S)) % p
        if np.linalg.matrix_rank(np.array(list(S)).astype(float)) < l: continue
        a = Lattice(N, p).char_poly(); b = LatticeGen(N, p).char_poly()
        print(l, len(S), a == b, a, b)
# braid A5 in K^5
braid5 = [[1 if k == i else 0 for k in range(5)] for i in range(5)] + [[1 if k == i else (-1 if k == j else 0) for k in range(5)] for i in range(5) for j in range(i + 1, 5)]
t = time.time(); L = LatticeGen(np.array(braid5) % p, p)
print('A5 braid chi', L.char_poly(), 'exps', L.exponents(), '%.2fs' % (time.time() - t))
t = time.time(); print('A5 freeness:', ExactModule(np.array(braid5) % 101, [1] * 15, 101).freeness(), '%.1fs' % (time.time() - t))
