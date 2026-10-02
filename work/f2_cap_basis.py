# Put the 13-point complete cap of PG(5,2) in a basis made of 6 of its points; re-check cap,
# completeness, blocking + minimality of the complement, and the root-sum non-DF argument.
import sys, json, itertools
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter/work')
from covering_check import is_covering, private_lines

n = 6
K = json.load(open('/Users/agent/Projects/TeraoHunter/work/caps_pg52.json'))['13']
vec = lambda p: [(p >> (n - 1 - i)) & 1 for i in range(n)]
# choose 6 independent points of K greedily
basis, span = [], {0}
for p in K:
    if p not in span:
        basis.append(p)
        span |= {s ^ p for s in span}
    if len(basis) == n:
        break
assert len(basis) == n
# coordinates of every point in this basis (solve over F_2 by brute force over subsets)
coord = {}
for mask in range(1, 2 ** n):
    v = 0
    for i in range(n):
        if mask >> i & 1:
            v ^= basis[i]
    coord[v] = tuple((mask >> i) & 1 for i in range(n))
Kc = sorted((coord[p] for p in K), key=lambda c: (sum(c), c))
print('cap in the new basis (13 vectors of F_2^6):')
for c in Kc:
    print('  ', ''.join(map(str, c)))
Ki = set(int(''.join(map(str, c)), 2) for c in Kc)
assert all((a ^ b) not in Ki for a in Ki for b in Ki if a != b), 'not a cap'
sums = {a ^ b for a in Ki for b in Ki if a != b}
print('cap:', True, ' complete:', all(p in sums or p in Ki for p in range(1, 64)), ' |K + K| =', len(sums))
A = [[(p >> (n - 1 - i)) & 1 for i in range(n)] for p in range(1, 64) if p not in Ki]
t = private_lines(A, 2)
print('|A| =', len(A), ' A* blocking:', is_covering(A, 2), ' min tangents per point:', min(t.values()),
      ' d =', len(A) - 31, '>= 16:', len(A) - 31 >= 16)
json.dump({'K': [list(c) for c in Kc], 'normals': A}, open('/Users/agent/Projects/TeraoHunter/work/f2_rank6_cap_nice.json', 'w'))
