# 3-parameter family: l1 = span(e1,e2) with points e1, e2, e1 + a e2;  l2 = span(e3,e4) with e3, e4,
# e3 + b e4;  transversal t through X1 = e1+e2, X2 = e3+e4 with points X1 + X2, X1 + c X2.
# A = planes spanned by >= 3 non-collinear of the 8 points.
import sys, random, itertools
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice, rank_mod, nullspace_mod, localization_l4
from hunter.ffree import free_rank3_ziegler, free_simple
from hunter.divisional import df_witnesses
def config(a, b, c, p):
    e = np.eye(4, dtype=np.int64)
    X1, X2 = (e[0] + e[1]) % p, (e[2] + e[3]) % p
    return [e[0], e[1], (e[0] + a * e[1]) % p, (X1 + X2) % p, (X1 + c * X2) % p, e[2], e[3], (e[2] + b * e[3]) % p]
def spanned(P, p):
    planes = {}
    for tri in itertools.combinations(range(len(P)), 3):
        M = np.array([P[i] for i in tri]) % p
        if rank_mod(M, p) < 3: continue
        K = nullspace_mod(M, p)
        v = K[0]
        e = next(i for i in range(4) if v[i]); v = v * pow(int(v[e]), p - 2, p) % p
        planes.setdefault(tuple(int(x) for x in v), set()).update(tri)
    keys = sorted(planes)
    return np.array(keys, dtype=np.int64), [sorted(planes[k]) for k in keys]
def summary(N, p):
    L = Lattice(N, p)
    return L, (len(N), tuple(sorted(len(S) for S in L.flats[2] if len(S) > 2)), tuple(sorted(len(S) for S in L.flats[3] if len(S) > 3)))
if __name__ == "__main__":
    p = int(sys.argv[1]) if (__name__ == "__main__" and len(sys.argv) > 1) else 1009
    rng = random.Random(1)
    ref7 = None
    types = {}
    for trial in range(30):
        a, b, c = (rng.randrange(2, p) for _ in range(3))
        N, sets = spanned(config(a, b, c, p), p)
        L, key = summary(N, p)
        types.setdefault(key, []).append((a, b, c))
    print('generic types over F_%d:' % p)
    for key, v in types.items():
        print('  planes', key[0], 'rank2', {m: key[1].count(m) for m in set(key[1])}, 'rank3', {m: key[2].count(m) for m in set(key[2])}, 'count', len(v))
    a, b, c = types[max(types, key=lambda k: len(types[k]))][0]
    N, sets = spanned(config(a, b, c, p), p); L = Lattice(N, p)
    print('generic: exps', L.exponents(), 'chi', L.char_poly())
    if L.exponents():
        print('  free?', free_simple(N, p, L)[0], 'DF witnesses', len(df_witnesses(N, p, L)))
    