import sys, json, glob, random
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter'); sys.path.insert(0, '/Users/agent/Projects/TeraoHunter/work')
import numpy as np
from hunter.core import Lattice, freeness, local_freeness_l4
from hunter.lift import LiftProblem
from lift_hit_lib import replay, p
hits = []
for fn in sorted(glob.glob('/Users/agent/Projects/TeraoHunter/work/family_search_*.out')):
    for line in open(fn):
        if line.startswith('HIT'):
            hits.append(json.loads(line[4:]))
k = int(sys.argv[1])
rec = hits[k]; recipe, m = rec['recipe'], rec['m']
ts = rec['nonfree_t'][0][0]
rng = random.Random(5)
P0 = LiftProblem(replay(recipe, 500), m, p)
pat = None
while pat is None:
    pat = P0.greedy(rng)
def lat(t):
    P = LiftProblem(replay(recipe, t), m, p)
    c = P.solve(pat, random.Random(9))
    A = P.arrangement(c)
    return P, A, Lattice(A, p)
Pg, Ag, Lg = lat(500)
Ps, As, Ls = lat(ts)
names = ['H'] + [f'L{L}.{i}' for (L, i) in Pg.planes]
for r in (2, 3):
    g = set(Lg.flats[r]); s = set(Ls.flats[r])
    extra = [sorted(S) for S in s - g]
    lost = [sorted(S) for S in g - s]
    print(f'rank {r}: flats only at special t:', [[names[i] for i in S] for S in extra])
    print(f'rank {r}: flats only at generic t:', [[names[i] for i in S] for S in lost][:10])
print('B lines (t generic):', [list(map(int, v)) for v in replay(recipe, 500)])
print('multiple points of B:', Pg.points, 'R:', Pg.R)
print('pattern:', pat)
bad = local_freeness_l4(As, p, Ls)
print('nonfree localizations at special t:', [[names[i] for i in sorted(S)] for S, st in bad])
