import sys, json, glob, random
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter'); sys.path.insert(0, '/Users/agent/Projects/TeraoHunter/work')
import numpy as np
from hunter.core import Lattice, freeness
from hunter.lift import LiftProblem
from hunter.divisional import df_witnesses
from hunter.restrict import ziegler_restriction
from lift_hit_lib import replay, p
hits = []
for fn in sorted(glob.glob('/Users/agent/Projects/TeraoHunter/work/family_search_*.out')):
    for line in open(fn):
        if line.startswith('HIT'):
            hits.append(json.loads(line[4:]))
for k in [5, 9, 10, 11, 13]:
    rec = hits[k]; recipe, m = rec['recipe'], rec['m']
    P0 = LiftProblem(replay(recipe, 500), m, p)
    rng = random.Random(7)
    seen = 0
    witness_classes = {}
    for attempt in range(120):
        if seen >= 8: break
        pat = P0.greedy(rng)
        if pat is None: continue
        c = P0.solve(pat, random.Random(1)); A = P0.arrangement(c); LA = Lattice(A, p)
        if freeness(A, p, LA)['status'] != 'free': continue
        seen += 1
        w = df_witnesses(A, p, LA)
        for i in w:
            cls = 'H' if i == 0 else P0.planes[i - 1][0]
            witness_classes[cls] = witness_classes.get(cls, 0) + 1
        print('  hit', k, 'lift', seen, 'witness planes', [('H' if i == 0 else P0.planes[i-1]) for i in w], flush=True)
    kk = rec['exps']
    print('hit', k, 'm', m, 'exps', kk, '| witness classes (class: count over 8 lifts):', witness_classes,
          '| mult of witness classes:', {c: (m[c] if c != 'H' else '-') for c in witness_classes})
