# Random greedy complete caps in PG(n-1, 2) (points = nonzero vectors of F_2^n as ints).
# A cap: no three distinct points with a ^ b ^ c = 0.  Complete: every outside point is a ^ b, a, b in K.
import random, sys, json
from collections import Counter

n = int(sys.argv[1]); trials = int(sys.argv[2])
N = (1 << n) - 1
rng = random.Random(1)
sizes = Counter()
best = {}
for t in range(trials):
    pts = list(range(1, N + 1))
    rng.shuffle(pts)
    K, sec = [], set()
    for p in pts:
        if p in sec or p in K:
            continue
        for a in K:
            sec.add(a ^ p)
        K.append(p)
    # completeness check
    assert all(p in sec or p in K for p in range(1, N + 1))
    k = len(K)
    sizes[k] += 1
    if k not in best:
        best[k] = sorted(K)
print('n=%d complete cap sizes (random greedy):' % n, sorted(sizes.items()))
json.dump({str(k): v for k, v in best.items()}, open('/Users/agent/Projects/TeraoHunter/work/caps_pg%d2.json' % (n - 1), 'w'))
