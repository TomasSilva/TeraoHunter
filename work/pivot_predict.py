# Predictor sweep for the pivot construction  B = B' + L0(s),  L0 through X1 in B', m0 >> 0.
# Criterion (conjectured, checked numerically): free <=> b2_off == e1*e2 (generic Euler exps on L0).
import sys, itertools, random, time
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice, multi_freeness
from hunter.multi import wakamiko, generic_rank2_exponents, local_exponents
p = 32003
rng = random.Random(5)
L13 = []
for v in itertools.product([-1,0,1], repeat=3):
    if any(v) and next(c for c in v if c) == 1:
        L13.append(list(v))
assert len(L13) == 13

def canon(lat_points, n):
    # cheap isomorphism-invariant of a line arrangement: sorted point sizes + per-line profiles
    prof = [[] for _ in range(n)]
    for S in lat_points:
        for i in S: prof[i].append(len(S))
    return (n, tuple(sorted(len(S) for S in lat_points)), tuple(sorted(tuple(sorted(x)) for x in prof)))

seen = set()
lib = []
for size in range(5, 9):
    combos = list(itertools.combinations(range(13), size))
    rng.shuffle(combos)
    for c in combos[:4000]:
        N = [L13[i] for i in c]
        try:
            L = Lattice(np.array(N) % p, p)
        except ValueError:
            continue
        pts = [sorted(S) for S in L.flats[2]]
        if max(len(S) for S in pts) == size:   # pencil
            continue
        key = canon(pts, size)
        if key in seen: continue
        seen.add(key)
        lib.append((N, pts))
print('library size', len(lib), flush=True)

K = 3
cands = 0
t0 = time.time()
for N, pts in lib:
    n = len(N)
    for X1 in pts:
        X1s = set(X1)
        others = [i for i in range(n) if i not in X1s]
        if len(others) < 3:      # need >= 4 points on L0 for moduli
            continue
        offpts = [S for S in pts if S != X1]
        for m in itertools.product(range(1, K + 1), repeat=n):
            M1 = sum(m[i] for i in X1)
            Ls = [M1] + [m[i] for i in others]
            M = sum(Ls)
            if 2 * max(Ls) >= M:
                continue
            e = (M // 2, M - M // 2)
            b2off = 0
            for S in offpts:
                if len(S) <= 3:
                    w = wakamiko([m[i] for i in S])
                else:
                    w = local_exponents(np.array(N) % p, m, S, p)
                b2off += w[0] * w[1]
                if b2off > e[0] * e[1]:
                    break
            if b2off == e[0] * e[1]:
                cands += 1
                print('CAND', N, 'X1', X1, 'm', m, 'L0mults', Ls, 'e', e, flush=True)
print('done candidates', cands, time.time() - t0)
