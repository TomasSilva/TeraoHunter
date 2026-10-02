import sys, json, itertools
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice, rank_mod, localization_l4
from hunter.ffree import free_rank3_ziegler, free_simple
from hunter.divisional import df_witnesses
q = 7
pts = [v for v in itertools.product(range(q), repeat=4) if any(v) and next(c for c in v if c) == 1]
P = [0, 1, 2, 8, 18, 57, 129, 224]
PP = np.array([pts[i] for i in P])
print('points:', PP.tolist())
# collinearities / coplanarities among the 8 points
for r, name in [(2, 'collinear triples'), (3, 'coplanar quadruples')]:
    cnt = []
    for c in itertools.combinations(range(8), r + 1):
        if rank_mod(PP[list(c)] % q, q) <= r:
            cnt.append(c)
    print(name, len(cnt), cnt[:20])
P_arr = np.array(pts, dtype=np.int64)
inc = (P_arr @ PP.T) % q == 0
Ls = [i for i in range(len(pts)) if inc[i].sum() >= 3 and rank_mod(PP[np.nonzero(inc[i])[0]] % q, q) >= 3]
N = P_arr[Ls] % q
L = Lattice(N, q)
print('planes', len(Ls), 'chi', L.char_poly(), 'exps', L.exponents())
print('rank-2 flat sizes', sorted(len(S) for S in L.flats[2] if len(S) > 2))
print('rank-3 flat sizes', sorted(len(S) for S in L.flats[3] if len(S) > 3))
print('points of P on each plane:', sorted(int(inc[i].sum()) for i in Ls))
print('DF witnesses:', df_witnesses(N, q, L))
from hunter.moduli import moduli_dimension
json.dump({'N': N.tolist(), 'q': q}, open('/Users/agent/Projects/TeraoHunter/work/arr26.json', 'w'))
