import sys, json, glob, random
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter'); sys.path.insert(0, '/Users/agent/Projects/TeraoHunter/work')
import numpy as np
from hunter.core import Lattice, freeness
from hunter.lift import LiftProblem
from hunter.divisional import df_witnesses
from lift_hit_lib import replay, p
recs = []
for fn in sorted(glob.glob('/Users/agent/Projects/TeraoHunter/work/family_search_*.jsonl')):
    for line in open(fn):
        recs.append(json.loads(line))
part, nparts = int(sys.argv[1]), int(sys.argv[2])
recs = recs[part::nparts]
for idx, rec in enumerate(recs):
    recipe, m = rec['recipe'], rec['m']
    B = replay(recipe, 500)
    if B is None: continue
    try:
        P0 = LiftProblem(B, m, p)
    except Exception as e:
        continue
    rng = random.Random(idx)
    nvalid = nondf = 0
    seen = set()
    for attempt in range(25):
        pat = P0.greedy(rng)
        if pat is None: continue
        key = tuple(sorted((X, tuple(t)) for X, t in pat))
        if key in seen: continue
        seen.add(key)
        c = P0.solve(pat, random.Random(1)); A = P0.arrangement(c)
        try:
            LA = Lattice(A, p)
        except ValueError:
            continue
        if freeness(A, p, LA)['status'] != 'free': continue
        nvalid += 1
        w = df_witnesses(A, p, LA)
        if not w:
            nondf += 1
            print('NON-DF FREE LIFT', json.dumps({'rec': rec, 'pattern': pat}), flush=True)
    print('rec', idx, 'n', rec['n'], 'm', m, 'special', [x[0] for x in rec['nonfree_t']], 'valid free lifts', nvalid, 'non-DF', nondf, flush=True)
