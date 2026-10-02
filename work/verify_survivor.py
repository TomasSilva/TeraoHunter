import sys, json, glob, random, ast
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter'); sys.path.insert(0, '/Users/agent/Projects/TeraoHunter/work')
import numpy as np
from hunter.core import Lattice, freeness, local_freeness_l4, rank_mod, multi_freeness
from hunter.lift import LiftProblem
from lift_hit_lib import replay, p
hits = []
for fn in sorted(glob.glob('/Users/agent/Projects/TeraoHunter/work/family_search_*.out')):
    for line in open(fn):
        if line.startswith('HIT'):
            hits.append(json.loads(line[4:]))
k = int(sys.argv[1]); which = int(sys.argv[2]) if len(sys.argv) > 2 else 0
rec = hits[k]; recipe, m = rec['recipe'], rec['m']
surv = [l for l in open(f'/Users/agent/Projects/TeraoHunter/work/replay_summary_{k}.out') if l.startswith('SURVIVOR')]
line = surv[which]
ts = int(line.split()[1]); pat = ast.literal_eval(line.split(' ', 2)[2])
pat = [(tuple(X), [tuple(x) for x in tup]) for X, tup in pat]
print('hit', k, 'm', m, 'special t', ts)
results = {}
for t in [500, 77, 123, ts]:
    P = LiftProblem(replay(recipe, t), m, p)
    M = P.system(pat); dim = len(P.planes) - rank_mod(M, p)
    for s in range(3):
        c = P.solve(pat, random.Random(1000 * t + s))
        if c is None or not P._distinct(c):
            print('  t', t, 'seed', s, 'degenerate'); continue
        A = P.arrangement(c)
        L = Lattice(A, p)
        fr = freeness(A, p, L, full=True)
        loc = local_freeness_l4(A, p, L)
        zr = multi_freeness(replay(recipe, t), m, p, tries=12)
        results[(t, s)] = L.signature()
        print(f'  t={t} seed={s} soldim={dim} planes={A.shape[0]} chi={L.char_poly()} A: {fr["status"]} exps={fr["exps"]}'
              f' hilbert={[(a,b,c) for a,b,c in fr["hilbert"] if b != c]} nonfreeLoc={len(loc)} | Ziegler restr: {zr}', flush=True)
ref = results.get((500, 0))
print('labelled lattice equal to (500,0):', {key: (sig == ref) for key, sig in results.items()})
