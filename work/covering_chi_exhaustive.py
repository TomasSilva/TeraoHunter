# Exhaustive consistency check of the covering lemma on PG(3,3): every orbit class whose normals
# meet every line must have chi = (t-1)(t-3)(t-9)(t-(|A|-13)) (Terao factorisation).
import sys, pickle, itertools, time
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice, integer_roots

q, l = 3, 4
pts = [v for v in itertools.product(range(q), repeat=l) if any(v) and next(c for c in v if c) == 1]
idx = {v: i for i, v in enumerate(pts)}
P = np.array(pts)


def norm(v):
    v = [x % q for x in v]
    c = next(x for x in v if x)
    inv = pow(c, q - 2, q)
    return tuple((x * inv) % q for x in v)


lines = set()
for a, b in itertools.combinations(pts, 2):
    lines.add(frozenset(idx[norm([x + t * y for x, y in zip(a, b)])] for t in range(q)) | {idx[b]})
lines = [np.array(sorted(L)) for L in lines]
print('lines of PG(3,3):', len(lines))
reps = pickle.load(open('/Users/agent/Projects/TeraoHunter/work/free_df_par_3_4_reps.pkl', 'rb'))
t0 = time.time()
ncov = bad = 0
sizes = {}
for S in reps:
    if len(S) < 13:
        continue
    mask = np.zeros(len(pts), bool)
    mask[list(S)] = True
    if not all(mask[L].any() for L in lines):
        continue
    ncov += 1
    sizes[len(S)] = sizes.get(len(S), 0) + 1
    ex = integer_roots(Lattice(P[list(S)] % q, q).char_poly())
    if ex is None or sorted(ex) != sorted([1, 3, 9, len(S) - 13]):
        bad += 1
        print('MISMATCH', S, ex, flush=True)
print('covering classes: %d, chi mismatches: %d (%.0fs)' % (ncov, bad, time.time() - t0))
print('by size:', sorted(sizes.items()))
