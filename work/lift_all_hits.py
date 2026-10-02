import sys, json, glob, random
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter'); sys.path.insert(0, '/Users/agent/Projects/TeraoHunter/work')
import numpy as np
from hunter.core import multi_freeness
from lift_hit_lib import replay, p
from work_lift import try_lift
hits = []
for fn in sorted(glob.glob('/Users/agent/Projects/TeraoHunter/work/family_search_*.out')):
    for line in open(fn):
        if line.startswith('HIT'):
            hits.append(json.loads(line[4:]))
k = int(sys.argv[1]); tries = int(sys.argv[2])
rec = hits[k]
recipe, m = rec['recipe'], rec['m']
ts = [500, 77] + [x[0] for x in rec['nonfree_t']]
for t in ts:
    B = replay(recipe, t)
    if B is None: continue
    st, g = multi_freeness(B, m, p, tries=12)
    P, found = try_lift(B, m, p, tries=tries, seed=t)
    print('hit', k, 'n', rec['n'], 'm', m, 't', t, 'multi', st, g, '| R', list(P.R.values()),
          '| lifts', len(found), sorted(set((f[4], len(f[3]), f[2]) for f in found)), flush=True)
