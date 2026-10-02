# Parametrized connecting-line families: points built by construction sequences
#   'll': new point = intersection of two connecting lines (determined)
#   'l1': new point on a connecting line, P = A + t*B (free parameter t)
# A_t = all lines through >= 2 of the points.  For free, non-DF A_t (generic t): scan all t in F_p
# for members with the same labelled lattice that are NOT free (rank-3 Terao counterexample;
# rank 3 is verified only up to 14 lines, so n >= 15 matters).
import sys, random, time, json, itertools
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice, _cross3
from hunter.ffree import free_rank3_ziegler
from hunter.divisional import rank3_df
from hunter.moduli import moduli_dimension

p = int(sys.argv[1]); T = float(sys.argv[2]); seed = int(sys.argv[3])
KLO, KHI = int(sys.argv[4]), int(sys.argv[5])
rng = random.Random(seed)
out = open(f'/Users/agent/Projects/TeraoHunter/work/pointc_{p}_{seed}.jsonl', 'w')


def norm(v):
    v = [int(x) % p for x in v]
    for c in v:
        if c:
            inv = pow(c, p - 2, p)
            return tuple(x * inv % p for x in v)
    return None


def connecting(points):
    lines = {}
    for a, b in itertools.combinations(range(len(points)), 2):
        L = norm(_cross3(np.array(points[a]), np.array(points[b]), p))
        if L is None:
            return None
        lines.setdefault(L, set()).update([a, b])
    return lines


def replay(recipe, t):
    pts = [norm(v) for v in ([1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1])]
    for st in recipe:
        if st[0] == 'll':
            (a, b), (c, d) = st[1], st[2]
            L1 = _cross3(np.array(pts[a]), np.array(pts[b]), p)
            L2 = _cross3(np.array(pts[c]), np.array(pts[d]), p)
            P = norm(_cross3(L1, L2, p))
        else:
            (a, b) = st[1]
            P = norm([(pts[a][i] + t * pts[b][i]) % p for i in range(3)])
        if P is None or P in pts:
            return None
        pts.append(P)
    lines = connecting(pts)
    if lines is None:
        return None
    keys = sorted(lines)
    return np.array(keys, dtype=np.int64), [sorted(lines[k]) for k in keys]


def random_recipe(k):
    for _ in range(200):
        recipe, npts = [], 4
        pos = rng.randrange(0, k - 4)
        ok = True
        for s in range(k - 4):
            pairs = list(itertools.combinations(range(npts), 2))
            if s == pos:
                recipe.append(('l1', rng.choice(pairs)))
            else:
                a, b = rng.sample(pairs, 2)
                if len(set(a) | set(b)) < 3:
                    ok = False
                    break
                recipe.append(('ll', a, b))
            npts += 1
        if not ok:
            continue
        r1, r2 = replay(recipe, rng.randrange(2, p)), replay(recipe, rng.randrange(2, p))
        if r1 is None or r2 is None:
            continue
        # labelled structure: the point-sets of lines must agree (same combinatorics)
        if sorted(map(tuple, r1[1])) != sorted(map(tuple, r2[1])):
            continue
        return recipe
    return None


def arr_for(recipe, t, template):
    r = replay(recipe, t)
    if r is None:
        return None
    N, sets = r
    if sorted(map(tuple, sets)) != template:
        return None
    # order lines consistently by their point sets
    order = sorted(range(len(sets)), key=lambda i: tuple(sets[i]))
    return N[order]


t0 = time.time()
st = {'families': 0, 'chi': 0, 'free': 0, 'nondf_free': 0, 'scanned': 0, 'mixed': 0}
while time.time() - t0 < T:
    k = rng.randint(KLO, KHI)
    recipe = random_recipe(k)
    if recipe is None:
        continue
    st['families'] += 1
    t1 = rng.randrange(2, p)
    r = replay(recipe, t1)
    if r is None:
        continue
    template = sorted(map(tuple, r[1]))
    N = arr_for(recipe, t1, template)
    if N is None or N.shape[0] < 15:
        continue
    try:
        L = Lattice(N, p)
    except ValueError:
        continue
    if L.exponents() is None:
        continue
    st['chi'] += 1
    if free_rank3_ziegler(N, p, L) != 'free':
        continue
    st['free'] += 1
    if rank3_df(N, p, L):
        continue
    st['nondf_free'] += 1
    sig = L.signature()
    md, _ = moduli_dimension(N, p, L)
    bad, same = [], 0
    for t in range(2, p):
        Nt = arr_for(recipe, t, template)
        if Nt is None:
            continue
        try:
            Lt = Lattice(Nt, p)
        except ValueError:
            continue
        if Lt.signature() != sig:
            continue
        same += 1
        if free_rank3_ziegler(Nt, p, Lt) != 'free':
            bad.append(t)
    st['scanned'] += 1
    rec = {'p': p, 'k': k, 'recipe': recipe, 'n_lines': int(N.shape[0]), 'exps': L.exponents(),
           't': sorted(len(F) for F in L.flats[2] if len(F) > 2), 'moduli_tangent': int(md),
           'same_lattice_t': same, 'nonfree_t': bad}
    out.write(json.dumps(rec) + '\n'); out.flush()
    print('NON-DF FREE FAMILY', json.dumps(rec), flush=True)
    if bad:
        st['mixed'] += 1
        print('!!! RANK-3 TERAO COUNTEREXAMPLE CANDIDATE', flush=True)
print(st, time.time() - t0, flush=True)
