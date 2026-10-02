import sys, json, glob, random
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter'); sys.path.insert(0, '/Users/agent/Projects/TeraoHunter/work')
import numpy as np
from hunter.core import Lattice, freeness, rank_mod
from hunter.lift import LiftProblem
from hunter.restrict import ziegler_restriction
from lift_hit_lib import replay, p
hits = []
for fn in sorted(glob.glob('/Users/agent/Projects/TeraoHunter/work/family_search_*.out')):
    for line in open(fn):
        if line.startswith('HIT'):
            hits.append(json.loads(line[4:]))
def poly_div(a, b):
    a = list(a); out = []
    while len(a) >= len(b):
        q = a[0] // b[0] if b[0] else 0
        if q * b[0] != a[0]: return None
        out.append(q)
        for i in range(len(b)): a[i] -= q * b[i]
        a.pop(0)
    return out if all(x == 0 for x in a) else None
for k in [5, 9, 10, 11, 13]:
    rec = hits[k]; recipe, m = rec['recipe'], rec['m']
    rng = random.Random(3)
    P0 = LiftProblem(replay(recipe, 500), m, p)
    found = 0
    for attempt in range(40):
        pat = P0.greedy(rng)
        if pat is None: continue
        c = P0.solve(pat, random.Random(1)); A = P0.arrangement(c)
        LA = Lattice(A, p)
        fA = freeness(A, p, LA)
        if fA['status'] != 'free': continue
        found += 1
        chiA = LA.char_poly()
        divs = []
        for i in range(A.shape[0]):
            N, mu, groups = ziegler_restriction(A, i, p)
            LR = Lattice(N, p)
            chiR = LR.char_poly()
            q = poly_div(chiA, chiR)
            fr = freeness(N, p, LR)['status']
            divs.append((i, len(mu), fr, q is not None))
        dfree = [d for d in divs if d[2] == 'free' and d[3]]
        print('hit', k, 'lift', found, 'A exps', fA['exps'], '| planes with free restriction & chi division:', [(d[0], d[1]) for d in dfree],
              '| free restrictions:', sum(1 for d in divs if d[2] == 'free'), '/', len(divs), flush=True)
        if found >= 3: break
