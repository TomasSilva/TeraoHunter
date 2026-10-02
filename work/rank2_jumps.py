import sys, itertools
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import module_generators
p = 1009   # p = 1 mod 3 and 1 mod 4 -> 1009 = 1 mod 12? 1009 % 12 = 1
assert p % 12 == 1
for m in [(2,2,2,2),(1,2,2,2),(2,2,2,3),(3,3,3,3),(1,1,2,2),(2,2,3,3),(1,2,2,3),(1,1,1,2),(1,1,2,3),(2,2,2,4)]:
    jumps = {}
    for lam in range(2, p):
        N = np.array([[1,0],[0,1],[1,p-1],[1,(p - lam) % p]], dtype=np.int64)
        gens, _ = module_generators(N, m, p)
        e = tuple(gens)
        jumps.setdefault(e, []).append(lam)
    summary = {e: (len(v) if len(v) > 20 else v) for e, v in jumps.items()}
    print(m, summary)
