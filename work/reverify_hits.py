import sys, json, glob, random
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter'); sys.path.insert(0, '/Users/agent/Projects/TeraoHunter/work')
import numpy as np
from hunter.core import Lattice, multi_freeness, module_generators
from lift_hit_lib import replay, p
for fn in sorted(glob.glob('/Users/agent/Projects/TeraoHunter/work/family_search_*.out')):
    for line in open(fn):
        if not line.startswith('HIT'): continue
        rec = json.loads(line[4:])
        recipe, m = rec['recipe'], rec['m']
        res = []
        for t, g in rec['nonfree_t']:
            st, gg = multi_freeness(replay(recipe, t), m, p, rng=random.Random(1000 + t), tries=16)
            full, _ = module_generators(replay(recipe, t), m, p, kmax=sum(m), stop_at_rank=False)
            res.append((t, st, full))
        print('n', rec['n'], 'm', m, 'generic', rec['exps'], '->', res, flush=True)
