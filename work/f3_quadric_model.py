# Clean model of the 22-plane example over F_3:
#   Q(y) = y1^2 + y2^2 + y3*y4 (elliptic quadric in the dual PG(3,3)), p = (0,0,1,0), T_p = {y4 = 0},
#   L = {y2 = y4 = 0} (a tangent line at p).   A* = L  u  {y : y4 != 0, Q(y) != 0}.
# Planes: x1, x3, x1 + x3, x1 - x3  and  a x1 + b x2 + c x3 + x4 with c != -(a^2 + b^2).
import sys, json, itertools
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter/work')
import numpy as np
import pynauty
from hunter.core import Lattice, integer_roots
from hunter.restrict import ziegler_restriction
from hunter.ffree import ExactModule
from covering_check import is_covering, private_lines, proj_points, normalize

q = 3
Qf = lambda y: (y[0] ** 2 + y[1] ** 2 + y[2] * y[3]) % q
pts = proj_points(q, 4)
print('points on Q:', sum(1 for y in pts if Qf(y) == 0), '(10 = elliptic)')
model = [(1, 0, 0, 0), (0, 0, 1, 0), (1, 0, 1, 0), (1, 0, 2, 0)]
model += [(a, b, c, 1) for a in range(3) for b in range(3) for c in range(3) if (c + a * a + b * b) % 3 != 0]
N = np.array(model)
L = Lattice(N, q)
print('|A| =', len(model), 'chi', L.char_poly(), 'exps', integer_roots(L.char_poly()))
print('free:', ExactModule(N, [1] * len(model), q).freeness())
print('|A^H|:', sorted(set(ziegler_restriction(N, i, q)[0].shape[0] for i in range(len(model)))))
print('covering:', is_covering(model, q), ' min tangents:', min(private_lines(model, q).values()))

# which PG(3,3)-class? compare nauty certificates with the 7 hits
NP = len(pts)
P = np.array(pts)
inc = (P @ P.T) % q == 0
adj = {NP + j: list(np.nonzero(inc[:, j])[0]) for j in range(NP)}
idx = {v: i for i, v in enumerate(pts)}


def cert(normals):
    S = set(idx[normalize([int(x) for x in v], q)] for v in normals)
    g = pynauty.Graph(2 * NP, directed=False, adjacency_dict=adj,
                      vertex_coloring=[set(range(NP)), set(NP + i for i in S), set(NP + i for i in range(NP) if i not in S)])
    return pynauty.certificate(g)


hits = json.load(open('/Users/agent/Projects/TeraoHunter/work/f3_rank4_free_nondf.json'))
cm = cert(model)
for h in hits:
    print('hit chunk %d (%d planes, exps %s): covering %s, min tangents %d, same class as model: %s' % (
        h['chunk'], h['k'], h['exps'], is_covering(h['normals'], q), min(private_lines(h['normals'], q).values()),
        cert(h['normals']) == cm))
