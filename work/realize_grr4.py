import sys, itertools, time
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice
from hunter.realize import Realizer
from hunter.ffree import free_simple
from hunter.divisional import df_witnesses
def grr(r, l, p):
    mu = sorted(set(a for a in range(1, p) if pow(a, r, p) == 1))
    assert len(mu) == r
    out = []
    for i in range(l):
        for j in range(i + 1, l):
            for z in mu:
                v = [0] * l; v[i] = 1; v[j] = (-z) % p
                out.append(v)
    return np.array(out, dtype=np.int64)
for (r, p) in [(3, 7), (3, 13), (4, 5), (4, 13), (3, 19), (6, 7), (5, 11), (4, 17)]:
    A = grr(r, 4, p)
    L = Lattice(A, p)
    t = time.time()
    st, g = free_simple(A, p, L)
    dfw = df_witnesses(A, p, L)
    sols = Realizer(A, p, max_solutions=5000).solutions()
    fr = []
    for s in sols:
        fr.append(free_simple(s, p)[0])
    print(f'G({r},{r},4) over F_{p}: n={len(A)} exps={L.exponents()} standard={st} DF-witnesses={len(dfw)} | realizations={len(sols)} free={fr.count("free")} nonfree={len(fr)-fr.count("free")} ({time.time()-t:.0f}s)', flush=True)
