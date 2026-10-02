# Stress test of the covering lemma: random hyperplane sets A in F_q^l whose normals meet every line
# of the dual space must be free with exponents (1, q, ..., q^(l-2), |A| - (1 + q + ... + q^(l-2))).
import sys, random, itertools, time
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.ffree import ExactModule
from covering_check import proj_points, lines_of

rng = random.Random(7)
for (q, l, trials) in [(5, 3, 40), (7, 3, 25), (2, 4, 60), (3, 4, 40), (2, 5, 25), (4, 3, 0)]:
    if q == 4:
        continue
    pts = proj_points(q, l)
    lines = lines_of(q, l)
    base = sum(q ** i for i in range(l - 1))
    ok = bad = 0
    t0 = time.time()
    for t in range(trials):
        # random minimal-ish blocking set: shuffle points, keep removing while still blocking
        S = set(pts)
        order = pts[:]
        rng.shuffle(order)
        target = rng.randint(0, len(pts))
        for p in order[:target]:
            S.discard(p)
            if not all(L & S for L in lines):
                S.add(p)
        if len(S) <= base:
            continue                     # non-essential (a hyperplane of the dual) -- skip
        N = np.array(sorted(S), dtype=np.int64)
        st, g = ExactModule(N, [1] * len(S), q).freeness()
        pred = sorted([q ** i for i in range(l - 1)] + [len(S) - base])
        if st == 'free' and sorted(g) == pred:
            ok += 1
        else:
            bad += 1
            print('MISMATCH', q, l, len(S), st, g, pred, flush=True)
    print('q=%d l=%d: %d agree, %d mismatch (%.0fs)' % (q, l, ok, bad, time.time() - t0), flush=True)
