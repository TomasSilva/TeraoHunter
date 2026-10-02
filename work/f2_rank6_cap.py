# F_2, rank 6: A = hyperplanes whose normals lie outside a complete cap K (|K| = 13) of PG(5,2).
# Predicted: free with exponents (1,2,4,8,16,19) (covering lemma) and not DF (no restriction has
# a characteristic polynomial dividing chi(A)).
import sys, json, time
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.lattice_gen import LatticeGen
from hunter.restrict import ziegler_restriction
from hunter.core import rank_mod, integer_roots
from hunter.divisional import poly_div

q, n = 2, 6
caps = json.load(open('/Users/agent/Projects/TeraoHunter/work/caps_pg52.json'))
K = caps['13']
vec = lambda p: [(p >> (n - 1 - i)) & 1 for i in range(n)]
A = [vec(p) for p in range(1, 2 ** n) if p not in K]
N = np.array(A, dtype=np.int64)
print('|K| =', len(K), '|A| =', len(A), 'rank', rank_mod(N, q))
# cap / completeness sanity
Ks = set(K)
assert all((a ^ b) not in Ks for a in K for b in K if a != b), 'not a cap'
sec = {a ^ b for a in K for b in K if a != b}
assert all(p in sec or p in Ks for p in range(1, 2 ** n)), 'not complete'
t0 = time.time()
chi = LatticeGen(N, q).char_poly()
print('chi(A) =', chi, 'roots', integer_roots(chi), '(%.0fs)' % (time.time() - t0))
divides = []
sizes = []
for i in range(len(A)):
    R, mu, g = ziegler_restriction(N, i, q)
    sizes.append(R.shape[0])
    c = LatticeGen(R, q).char_poly()
    if poly_div(chi, c) is not None:
        divides.append((i, c))
print('|A^H| over H:', sorted(set(sizes)), ' (full restriction would be 31)')
print('restrictions with chi(A^H) | chi(A):', divides, '(%.0fs)' % (time.time() - t0))
json.dump({'K': K, 'normals': A, 'chi': chi}, open('/Users/agent/Projects/TeraoHunter/work/f2_rank6_cap.json', 'w'))
