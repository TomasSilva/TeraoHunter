# Step 2 prefilter: for a multiarrangement hit (B, m) with exps (d1,d2,d3), a lift's plane of a
# multiplicity-1 class L has |A^P| = (n-1) - sum_{X in L, |B_X|>=3} (concurrency excess at X);
# with only triple concurrencies at 3-line points, the class-L plane is in all R_X triple lines.
# DF-size: 1 + d_i + d_j.  Report hits where no multiplicity-1 class (with only 3-line points)
# is forced to a DF size, and |B| is not a DF size (or B is not DF).
import sys, json, glob
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter'); sys.path.insert(0, '/Users/agent/Projects/TeraoHunter/work')
import numpy as np
from hunter.lift import LiftProblem
from lift_hit_lib import replay, p
hits = []
for fn in sorted(glob.glob('/Users/agent/Projects/TeraoHunter/work/family_search_*.out')):
    for line in open(fn):
        if line.startswith('HIT'):
            hits.append(json.loads(line[4:]))
for k, rec in enumerate(hits):
    m = rec['m']; d = sorted(rec['exps'])
    n = 1 + sum(m)
    dfsizes = {1 + d[0] + d[1], 1 + d[0] + d[2], 1 + d[1] + d[2]}
    P = LiftProblem(replay(rec['recipe'], 500), m, p)
    forced, undetermined = [], []
    for L in range(len(m)):
        if m[L] != 1: continue
        pts = [X for X in P.points if L in X]
        if any(len(X) >= 4 for X in pts):
            undetermined.append(L); continue
        tau = sum(P.R[X] for X in pts)
        size = (n - 1) - tau
        if size in dfsizes: forced.append((L, size))
    hsize = len(m)
    print('hit', k, 'n_lift', n, 'm', m, 'exps', d, 'DF sizes', sorted(dfsizes), '| forced DF classes', forced,
          '| undetermined m=1 classes', undetermined, '| |B|=', hsize, ('(DF size!)' if hsize in dfsizes else ''),
          '  PASS' if not forced and hsize not in dfsizes else '')
