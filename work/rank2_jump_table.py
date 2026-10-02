# For 4 points (0, inf, 1, lam) with multiplicities (a,b,c,d), find lam values where exponents jump.
import sys, itertools
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import module_generators
p = 1009
res = {}
for ms in itertools.combinations_with_replacement(range(1, 6), 4):
    M = sum(ms)
    if 2 * max(ms) >= M: continue
    # all distinct assignments up to symmetry: fix the multiset; positions matter up to S4 acting on lam
    jumps = {}
    for lam in range(2, p):
        N = np.array([[1,0],[0,1],[1,p-1],[1,(p - lam) % p]], dtype=np.int64)
        gens, _ = module_generators(N, ms, p)
        e = tuple(gens[:2])
        jumps.setdefault(e, []).append(lam)
    gen = max(jumps, key=lambda e: len(jumps[e]))
    others = {e: v for e, v in jumps.items() if e != gen}
    if others:
        print(ms, 'generic', gen, 'jumps', {e: v for e, v in others.items()}, flush=True)
    else:
        print(ms, 'generic', gen, 'no jumps', flush=True)
