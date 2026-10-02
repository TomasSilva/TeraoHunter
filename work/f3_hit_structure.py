# Symmetry and geometry of the 22-plane free non-DF arrangement over F_3.
import sys, json, itertools
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
import pynauty
from collections import Counter

q = 3
rec = json.loads(open('/Users/agent/.claude/jobs/ef50d01d/tmp/hit.json').read())
pts = [v for v in itertools.product(range(q), repeat=4) if any(v) and next(c for c in v if c) == 1]
NP = len(pts)
idx = {v: i for i, v in enumerate(pts)}
A = sorted(idx[tuple(v)] for v in rec['normals'])
P = np.array(pts)
inc = (P @ P.T) % q == 0          # point i (primal) on plane j (dual normal) iff dot = 0
adj = {NP + j: list(np.nonzero(inc[:, j])[0]) for j in range(NP)}
g = pynauty.Graph(2 * NP, directed=False, adjacency_dict=adj,
                  vertex_coloring=[set(range(NP)), set(NP + i for i in A), set(NP + i for i in range(NP) if i not in A)])
gens, grpsize1, grpsize2, orbits, numorb = pynauty.autgrp(g)
order = grpsize1 * 10 ** grpsize2
print('|Stab(A)| in PGL(4,3) =', order)
orb = {}
for v, o in enumerate(orbits):
    orb.setdefault(o, []).append(v)
plane_orbits = [sorted(pts[v - NP] for v in vs) for vs in orb.values() if vs[0] >= NP and vs[0] - NP in A]
print('orbits on A:', [len(o) for o in plane_orbits])
for o in plane_orbits:
    print('  ', o)
pt_orbits = [vs for vs in orb.values() if vs[0] < NP]
print('orbits on points:', sorted(len(o) for o in pt_orbits))
# multiplicity of each primal point = number of planes of A through it
mult = {i: sum(1 for j in A if inc[i, j]) for i in range(NP)}
for vs in pt_orbits:
    print('  point orbit size %d, planes of A through each: %d, e.g. %s' % (len(vs), mult[vs[0]], pts[vs[0]]))
