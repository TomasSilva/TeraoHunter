import sys, json, glob
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter'); sys.path.insert(0, '/Users/agent/Projects/TeraoHunter/work')
import numpy as np
from hunter.core import Lattice, multi_freeness
from hunter.multi import local_exponents
from lift_hit_lib import replay, p
for fn in sorted(glob.glob('/Users/agent/Projects/TeraoHunter/work/family_search_*.out')):
    for line in open(fn):
        if not line.startswith('HIT'): continue
        rec = json.loads(line[4:])
        recipe, m = rec['recipe'], rec['m']
        def locs(t):
            B = replay(recipe, t)
            L = Lattice(B, p)
            return {tuple(sorted(S)): local_exponents(B, m, S, p) for S in L.flats[2] if len(S) >= 3}
        g = locs(500 if replay(recipe, 500) is not None else 77)
        kinds = []
        for t, gens in rec['nonfree_t']:
            s = locs(t)
            kinds.append('GLOBAL' if s == g else 'local')
        print(kinds, 'n', rec['n'], 'm', m, 'exps', rec['exps'], 'special t', [x[0] for x in rec['nonfree_t']], 'pts', rec['mult_points'])
