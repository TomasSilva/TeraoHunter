import sys, json, glob, random
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter'); sys.path.insert(0, '/Users/agent/Projects/TeraoHunter/work')
import numpy as np
from hunter.core import Lattice, freeness, multi_freeness, rank_mod
from hunter.lift import LiftProblem
from hunter.restrict import ziegler_restriction, heavy_lines
from lift_hit_lib import replay, p
hits = []
for fn in sorted(glob.glob('/Users/agent/Projects/TeraoHunter/work/family_search_*.out')):
    for line in open(fn):
        if line.startswith('HIT'):
            hits.append(json.loads(line[4:]))
k = int(sys.argv[1])
rec = hits[k]; recipe, m = rec['recipe'], rec['m']
rng = random.Random(3)
P0 = LiftProblem(replay(recipe, 500), m, p)
while True:
    pat = P0.greedy(rng)
    if pat is None: continue
    c = P0.solve(pat, random.Random(1)); A = P0.arrangement(c)
    if freeness(A, p)['status'] == 'free': break
names = ['H'] + [f'L{L}.{i}' for (L, i) in P0.planes]
print('hit', k, 'A has', A.shape[0], 'planes; free exps', freeness(A, p)['exps'])
for i in range(A.shape[0]):
    N, mu, groups = ziegler_restriction(A, i, p)
    st, g = multi_freeness(N, mu, p)
    hv = heavy_lines(N, mu, p)
    # combinatorially forced if a locally heavy line has Euler-Ziegler with <= 3 points or a dominant point
    forced = [h for h in hv if len(h[1]) <= 3 or 2 * max(h[1]) >= sum(h[1])]
    print(f'{names[i]:>6}: lines={len(mu)} mults={sorted(mu, reverse=True)} free={st} {g} heavy={[(h[0], h[1]) for h in hv]} FORCED={bool(forced)}')
from hunter.moduli import moduli_dimension
print('--- moduli dimension of each restriction (line arrangement, tangent space) ---')
for i in range(A.shape[0]):
    N, mu, groups = ziegler_restriction(A, i, p)
    d, _ = moduli_dimension(N, p)
    print(f'{names[i]:>6}: moduli_dim={d}')
