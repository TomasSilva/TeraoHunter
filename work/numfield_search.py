# Search for Terao counterexamples over quadratic fields K = Q(alpha): connecting-line arrangements
# A(P) of point sets P with coordinates in a small set E_K, computed modulo split primes P1, P2
# (lattice accepted only if identical mod both).  Keeps free, non-DF lattices (one record per lattice).
import sys, itertools, random, json, time
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np, pynauty
from hunter.core import Lattice, freeness
from hunter.divisional import rank3_df
from hunter.moduli import moduli_dimension

FIELDS = {'i': (1, 0, 1), 'sqrt2': (1, 0, -2), 'omega': (1, 1, 1), 'phi': (1, -1, -1),
          'sqrt-2': (1, 0, 2), 'sqrt3': (1, 0, -3), 'sqrt-7': (1, 0, 7), 'sqrt5': (1, 0, -5)}
name = sys.argv[1]; secs = int(sys.argv[2]); seed = int(sys.argv[3])
a2, a1, a0 = FIELDS[name]


def split_primes(lo, cnt):
    out, p = [], lo
    while len(out) < cnt:
        p += 1
        if all(p % d for d in range(2, int(p ** .5) + 1)):
            r = [x for x in range(p) if (a2 * x * x + a1 * x + a0) % p == 0]
            if r: out.append((p, r[0]))
    return out
(P1, r1), (P2, r2) = split_primes(30000, 1)[0], split_primes(50000, 1)[0]
# K-elements as pairs (u, v) = u + v*alpha
E = [(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (1, -1), (-1, 1), (-1, -1), (2, 0)]
def red(x, p, r): return (x[0] + x[1] * r) % p
def mulK(x, y):  # (u+v a)(s+t a), a^2 = -(a1 a + a0)/a2 with a2 = 1
    u, v = x; s, t = y
    w = v * t
    return (u * s - a0 * w, u * t + v * s - a1 * w)
pts = set()
for x in itertools.product(E, repeat=3):
    if all(c == (0, 0) for c in x): continue
    pts.add(x)
pts = list(pts)
def proj_key(x, p, r):
    v = [red(c, p, r) for c in x]
    if not any(v): return None
    k = next(c for c in v if c); inv = pow(k, p - 2, p)
    return tuple(c * inv % p for c in v)
uniq = {}
for x in pts:
    k = proj_key(x, P1, r1)
    if k and k not in uniq: uniq[k] = x
S = list(uniq.values())
alpha_map = lambda x: (mulK(x[0], (0, 1)), x[1], x[2])          # (x:y:z) -> (alpha x : y : z)
cyc = lambda x: (x[1], x[2], x[0])


def arrangement(P, p, r):
    V = [np.array([red(c, p, r) for c in x], dtype=np.int64) for x in P]
    lines = {}
    for i, j in itertools.combinations(range(len(V)), 2):
        u, w = V[i], V[j]
        l = np.array([u[1]*w[2]-u[2]*w[1], u[2]*w[0]-u[0]*w[2], u[0]*w[1]-u[1]*w[0]]) % p
        if not l.any(): return None
        k = next(c for c in l if c); l = tuple(int(c) * pow(int(k), p - 2, p) % p for c in l)
        lines[l] = 1
    return np.array(list(lines), dtype=np.int64)


def sig(L):
    return tuple(sorted(len(F) for F in L.flats[2]))


def cert(L, n):
    pts_ = [F for F in L.flats[2]]
    adj = {n + k: list(F) for k, F in enumerate(pts_)}
    g = pynauty.Graph(n + len(pts_), adjacency_dict=adj, vertex_coloring=[set(range(n)), set(range(n, n + len(pts_)))])
    return pynauty.certificate(g)


rng = random.Random(seed)
Esmall = [(0, 0), (1, 0), (-1, 0), (0, 1)]
R = list({proj_key(x, P1, r1): x for x in itertools.product(Esmall, repeat=3) if proj_key(x, P1, r1)}.values())
# special lines: coefficient vectors in E (same reduction machinery), deduped
LK = list({proj_key(x, P1, r1): x for x in itertools.product(E[:9], repeat=3) if proj_key(x, P1, r1)}.values())
LK1 = np.array([[red(c, P1, r1) for c in x] for x in LK], dtype=np.int64)
seen, t0, tried, stats = set(), time.time(), 0, {'chi': 0, 'free': 0, 'nondf': 0}
out = open('/Users/agent/Projects/TeraoHunter/work/numfield/hits_%s_%d.jsonl' % (name, seed), 'w')


def consider(idx_or_P, kind):
    global stats
    if kind == 'pts':
        N = arrangement(idx_or_P, P1, r1)
        if N is None: return
    else:
        N = LK1[list(idx_or_P)]
    if len(N) < 13 or np.linalg.matrix_rank(N.astype(float)) < 3: return
    try:
        L = Lattice(N, P1)
    except ValueError:
        return
    ex = L.exponents()
    if ex is None: return
    stats['chi'] += 1
    if rank3_df(N, P1, L): return
    if freeness(N, P1, L)['status'] != 'free': return
    stats['free'] += 1
    if kind == 'pts':
        N2 = arrangement(idx_or_P, P2, r2)
    else:
        N2 = np.array([[red(c, P2, r2) for c in LK[i]] for i in idx_or_P], dtype=np.int64)
    if N2 is None or len(N2) != len(N): return
    try:
        if sig(Lattice(N2, P2)) != sig(L): return
    except ValueError:
        return
    c = cert(L, len(N))
    if c in seen: return
    seen.add(c); stats['nondf'] += 1
    md, _ = moduli_dimension(N, P1, L)
    rec = {'field': name, 'kind': kind, 'data': idx_or_P if kind == 'pts' else [LK[i] for i in idx_or_P],
           'n': len(N), 'exps': list(ex), 't': sorted(len(F) for F in L.flats[2] if len(F) > 2), 'moduli': int(md)}
    out.write(json.dumps(rec) + '\n'); out.flush()


def through_count(Lset, j):
    # number of existing multiple points (pairwise intersections) the line j passes through
    cnt = set()
    lj = LK1[j]
    for a, b in itertools.combinations(Lset, 2):
        u, w = LK1[a], LK1[b]
        x = np.array([u[1]*w[2]-u[2]*w[1], u[2]*w[0]-u[0]*w[2], u[0]*w[1]-u[1]*w[0]]) % P1
        if int(x @ lj % P1) == 0 and x.any():
            k = next(c for c in x if c); cnt.add(tuple(int(c) * pow(int(k), P1 - 2, P1) % P1 for c in x))
    return len(cnt)


while time.time() - t0 < secs:
    tried += 1
    if rng.random() < 0.5:
        consider(rng.sample(R, min(len(R), rng.randint(6, 10))), 'pts')
    else:
        cur = rng.sample(range(len(LK)), 5)
        while len(cur) < 26:
            cand = rng.sample([j for j in range(len(LK)) if j not in cur], 25)
            sc = [(through_count(cur, j) + rng.random(), j) for j in cand]
            cur.append(max(sc)[1])
            if len(cur) >= 13: consider(list(cur), 'lines')
print(name, seed, 'tried', tried, stats, 'distinct free non-DF lattices', len(seen), flush=True)
