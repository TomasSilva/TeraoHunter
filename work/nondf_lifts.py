import sys, json, glob, random, time
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter'); sys.path.insert(0, '/Users/agent/Projects/TeraoHunter/work')
import numpy as np
from hunter.core import Lattice, freeness, rank_mod
from hunter.lift import LiftProblem
from hunter.divisional import df_witnesses
from lift_hit_lib import replay, p
hits = []
for fn in sorted(glob.glob('/Users/agent/Projects/TeraoHunter/work/family_search_*.out')):
    for line in open(fn):
        if line.startswith('HIT'):
            hits.append(json.loads(line[4:]))
k = int(sys.argv[1]); tries = int(sys.argv[2])
rec = hits[k]; recipe, m = rec['recipe'], rec['m']
specials = [x[0] for x in rec['nonfree_t']]
P0 = LiftProblem(replay(recipe, 500), m, p)
rng = random.Random(1000 + k)
seen = set(); nvalid = nondf = 0; t0 = time.time()
for trial in range(tries):
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
    if w: continue
    nondf += 1
    sig0 = LA.signature()
    verdicts = []
    for ts in specials:
        Ps = LiftProblem(replay(recipe, ts), m, p)
        cs = Ps.solve(pat, random.Random(2))
        try:
            Ls = Lattice(Ps.arrangement(cs), p)
            same = Ls.signature() == sig0
            fs = freeness(Ps.arrangement(cs), p, Ls)['status'] if same else '-'
        except ValueError:
            same, fs = False, 'degenerate'
        verdicts.append((ts, same, fs))
    print('NON-DF FREE LIFT hit', k, 'm', m, 'replay at special t:', verdicts, flush=True)
    if any(s and f != 'free' for _, s, f in verdicts):
        print('!!! CANDIDATE COUNTEREXAMPLE', json.dumps({'hit': k, 'pattern': pat}), flush=True)
print('hit', k, 'm', m, 'valid free lifts', nvalid, 'non-DF', nondf, 'time %.0f' % (time.time() - t0), flush=True)
