# For a hit: find lift patterns at generic t, replay each pattern at other t (incl. special),
# compare labelled lattices of A_t, and test freeness of A_t directly in rank 4.
import sys, json, glob, random
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter'); sys.path.insert(0, '/Users/agent/Projects/TeraoHunter/work')
import numpy as np
from hunter.core import Lattice, freeness, local_freeness_l4, multi_freeness, rank_mod
from hunter.lift import LiftProblem
from lift_hit_lib import replay, p
hits = []
for fn in sorted(glob.glob('/Users/agent/Projects/TeraoHunter/work/family_search_*.out')):
    for line in open(fn):
        if line.startswith('HIT'):
            hits.append(json.loads(line[4:]))
k = int(sys.argv[1]); tries = int(sys.argv[2])
rec = hits[k]
recipe, m = rec['recipe'], rec['m']
specials = [x[0] for x in rec['nonfree_t']]
t0 = 500
rng = random.Random(k)
P0 = LiftProblem(replay(recipe, t0), m, p)
def build(t, pattern, seed):
    P = LiftProblem(replay(recipe, t), m, p)
    if P.points != P0.points:
        return None, 'B-lattice changed'
    M = P.system(pattern)
    rk = rank_mod(M, p)
    c = P.solve(pattern, random.Random(seed))
    if c is None or not P._distinct(c):
        return None, 'degenerate'
    A = P.arrangement(c)
    try:
        LA = Lattice(A, p)
    except ValueError:
        return None, 'repeated plane'
    return (A, LA, rk, P.savings(c) == P.R), 'ok'
seen = set()
for trial in range(tries):
    pat = P0.greedy(rng)
    if pat is None: continue
    key = tuple(sorted((X, tuple(t)) for X, t in pat))
    if key in seen: continue
    seen.add(key)
    base, why = build(t0, pat, 1)
    if base is None: continue
    A0, L0, rk0, ok0 = base
    f0 = freeness(A0, p, L0)
    sig0 = L0.signature()
    line = [f'pattern#{len(seen)} rank {rk0} savingsOK {ok0} A_{t0}: {f0["status"]} {f0["exps"]}']
    for t in [77, 123] + specials:
        res, why = build(t, pat, 2)
        if res is None:
            line.append(f't={t}: {why}'); continue
        A, LA, rk, ok = res
        same = LA.signature() == sig0
        fr = freeness(A, p, LA)
        loc = local_freeness_l4(A, p, LA)
        line.append(f't={t}: rank {rk} savingsOK {ok} sameLattice {same} A: {fr["status"]} {fr["exps"]} nonfreeLoc {len(loc)}')
    print(' | '.join(line), flush=True)
print('distinct patterns', len(seen))
