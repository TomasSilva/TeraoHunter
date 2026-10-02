# Which covering classes of PG(3,3) are minimal blocking sets?  Compare with the 7 free non-DF classes.
import sys, pickle, itertools, json
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter/work')
import numpy as np
exec(open('/Users/agent/Projects/TeraoHunter/work/covering_chi_exhaustive.py').read().split("reps = pickle")[0])
reps = pickle.load(open('/Users/agent/Projects/TeraoHunter/work/free_df_par_3_4_reps.pkl', 'rb'))
hits = json.load(open('/Users/agent/Projects/TeraoHunter/work/f3_rank4_free_nondf.json'))
hitsets = {tuple(h['S']) for h in hits}
minimal = {}
for S in reps:
    if len(S) < 13:
        continue
    mask = np.zeros(len(pts), bool)
    mask[list(S)] = True
    cnt = [mask[L].sum() for L in lines]
    if min(cnt) == 0:
        continue
    tangent_pts = set()
    for L, c in zip(lines, cnt):
        if c == 1:
            tangent_pts.add(int(L[mask[L]][0]))
    if len(tangent_pts) == len(S):
        minimal.setdefault(len(S), []).append(tuple(S))
print('minimal blocking sets by size:', sorted((k, len(v)) for k, v in minimal.items()))
big = [S for k, v in minimal.items() if k >= 22 for S in v]
print('minimal of size >= 22:', len(big), ' all among the 7 free non-DF:', all(S in hitsets for S in big),
      ' 7 hits all minimal:', all(any(S == T for T in big) for S in hitsets))
