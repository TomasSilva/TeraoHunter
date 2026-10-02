# Independent check of Lemma 1 with Macaulay2: random covering classes of PG(3,3) (from the exhaustive
# list) and random covering subsets of PG(4,2) must be free with exps (1,q,...,q^(l-2), |A|-[l-1]_q).
import sys, pickle, itertools, random, time
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter/work')
import numpy as np
from hunter.m2 import m2_free
from covering_check import proj_points, lines_of

rng = random.Random(3)


def predicted(q, l, n):
    return sorted([q ** i for i in range(l - 1)] + [n - sum(q ** i for i in range(l - 1))])


def check(q, l, normals):
    res = m2_free([list(map(int, v)) for v in normals], field='ZZ/%d' % q, timeout=3000)
    gens = sorted(int(x) for x in res['gens'].strip('{} ').split(',')) if res.get('gens') else None
    return res.get('free'), gens


# PG(3,3): sample from the exhaustive orbit list
q, l = 3, 4
pts = proj_points(q, l)
P = np.array(pts)
lines = [np.array(sorted(pts.index(p) for p in L)) for L in lines_of(q, l)]
reps = pickle.load(open('/Users/agent/Projects/TeraoHunter/work/free_df_par_3_4_reps.pkl', 'rb'))
cov = []
for S in reps:
    if len(S) < 14:
        continue
    mask = np.zeros(len(pts), bool)
    mask[list(S)] = True
    if all(mask[L].any() for L in lines):
        cov.append(S)
sample = rng.sample(cov, 25)
t0 = time.time()
ok = 0
for S in sample:
    fr, g = check(q, l, P[list(S)])
    good = fr and g == predicted(q, l, len(S))
    ok += good
    if not good:
        print('MISMATCH PG(3,3)', len(S), fr, g, predicted(q, l, len(S)), flush=True)
print('PG(3,3): %d/%d sampled covering classes free with predicted exponents (%.0fs)' % (ok, len(sample), time.time() - t0), flush=True)

# PG(4,2): random minimal-ish blocking sets (complements of random caps) and random supersets
q, l = 2, 5
pts = proj_points(q, l)
lines = lines_of(q, l)
t0 = time.time()
ok = tot = 0
for t in range(15):
    S = set(pts)
    order = pts[:]
    rng.shuffle(order)
    for p in order[:rng.randint(5, len(pts))]:
        S.discard(p)
        if not all(L & S for L in lines):
            S.add(p)
    if len(S) <= 15:
        continue
    tot += 1
    fr, g = check(q, l, sorted(S))
    good = fr and g == predicted(q, l, len(S))
    ok += good
    if not good:
        print('MISMATCH PG(4,2)', len(S), fr, g, predicted(q, l, len(S)), flush=True)
print('PG(4,2): %d/%d random covering sets free with predicted exponents (%.0fs)' % (ok, tot, time.time() - t0), flush=True)
