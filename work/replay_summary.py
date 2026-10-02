import sys, json, glob, random
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter'); sys.path.insert(0, '/Users/agent/Projects/TeraoHunter/work')
import numpy as np
from hunter.core import Lattice, freeness, local_freeness_l4, rank_mod
from hunter.lift import LiftProblem
from lift_hit_lib import replay, p
hits = []
for fn in sorted(glob.glob('/Users/agent/Projects/TeraoHunter/work/family_search_*.out')):
    for line in open(fn):
        if line.startswith('HIT'):
            hits.append(json.loads(line[4:]))
k = int(sys.argv[1]); tries = int(sys.argv[2])
rec = hits[k]; recipe, m = rec['recipe'], rec['m']
specials = [x[0] for x in rec['nonfree_t']]
rng = random.Random(100 + k)
P0 = LiftProblem(replay(recipe, 500), m, p)
seen = set(); stats = {}
survivors = []
for trial in range(tries):
    pat = P0.greedy(rng)
    if pat is None: continue
    key = tuple(sorted((X, tuple(t)) for X, t in pat))
    if key in seen: continue
    seen.add(key)
    M0 = P0.system(pat); dim0 = len(P0.planes) - rank_mod(M0, p)
    c0 = P0.solve(pat, random.Random(1)); A0 = P0.arrangement(c0); L0 = Lattice(A0, p)
    sig0 = L0.signature()
    for ts in specials:
        Ps = LiftProblem(replay(recipe, ts), m, p)
        if Ps.points != P0.points: continue
        Ms = Ps.system(pat); dims = len(Ps.planes) - rank_mod(Ms, p)
        same_any = False
        for s in range(6):   # several random solutions (matters if dims > 4)
            cs = Ps.solve(pat, random.Random(s))
            if cs is None or not Ps._distinct(cs): continue
            try:
                Ls = Lattice(Ps.arrangement(cs), p)
            except ValueError:
                continue
            if Ls.signature() == sig0:
                same_any = True
                survivors.append((pat, ts, cs))
                break
        stats.setdefault((dim0, dims, same_any), 0)
        stats[(dim0, dims, same_any)] += 1
print('hit', k, 'n', rec['n'], 'm', m, 'special', specials, 'patterns', len(seen), 'stats (dim_generic, dim_special, sameLattice):', stats, flush=True)
for pat, ts, cs in survivors[:3]:
    print('SURVIVOR', ts, pat, flush=True)
