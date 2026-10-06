# Structured search over K = Q(alpha): deletions/additions of special (reflection-type) arrangements.
# Exact arithmetic in Z[alpha]; lattice accepted only if identical modulo two large split primes.
import sys, itertools, random, json, time
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np, pynauty
from hunter.core import Lattice, freeness
from hunter.divisional import rank3_df
from hunter.moduli import moduli_dimension
FIELDS = {'i': (0, 1), 'omega': (1, 1), 'phi': (-1, -1), 'sqrt2': (0, -2)}   # alpha^2 + a1 alpha + a0 = 0
name, base_name, secs, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
a1, a0 = FIELDS[name]
def mul(x, y):
    u, v = x; s, t = y; w = v * t
    return (u * s - a0 * w, u * t + v * s - a1 * w)
def sub(x, y): return (x[0] - y[0], x[1] - y[1])
def cross(a, b): return tuple(sub(mul(a[i], b[j]), mul(a[j], b[i])) for i, j in ((1, 2), (2, 0), (0, 1)))
Z, O, A = (0, 0), (1, 0), (0, 1)
def neg(x): return (-x[0], -x[1])
def pw(x, k):
    r = O
    for _ in range(k): r = mul(r, x)
    return r
def monomial(r, zeta, with_coords):
    L = [(O, Z, Z), (Z, O, Z), (Z, Z, O)] if with_coords else []
    for k in range(r):
        c = neg(pw(zeta, k))
        L += [(O, c, Z), (Z, O, c), (c, Z, O)]
    return L
phi = A
BASES = {
    'B3': [(O, Z, Z), (Z, O, Z), (Z, Z, O), (O, O, Z), (O, neg(O), Z), (Z, O, O), (Z, O, neg(O)), (O, Z, O), (O, Z, neg(O))],
    'G443': monomial(4, A, False), 'G413': monomial(4, A, True),          # i
    'G333': monomial(3, A, False), 'G313': monomial(3, A, True),          # omega (cube root)
    'G613': monomial(6, neg(mul(A, A)), True),                            # zeta6 = -omega^2
    'H3': [(O, Z, Z), (Z, O, Z), (Z, Z, O)] + [v for s1 in (1, -1) for s2 in (1, -1)
            for v in [(O, (s1, 0) if False else (s1 * 0, s1), (s2 * -1, s2)) ]],
}
# H3 normals: cyclic shifts of (1, +-phi, +-phi^{-1}), phi^{-1} = phi - 1
h = []
for s1 in (1, -1):
    for s2 in (1, -1):
        v = (O, (0, s1), (-s2, s2))
        h += [v, (v[2], v[0], v[1]), (v[1], v[2], v[0])]
BASES['H3'] = [(O, Z, Z), (Z, O, Z), (Z, Z, O)] + h
base = BASES[base_name]
def primes(lo):
    p = lo
    while True:
        p += 1
        if all(p % d for d in range(2, int(p ** .5) + 1)):
            r = [x for x in range(p) if (x * x + a1 * x + a0) % p == 0]
            if r: return p, r[0]
(P1, r1), (P2, r2) = primes(30000), primes(50000)
def red(v, p, r): return [(c[0] + c[1] * r) % p for c in v]
def key(v, p, r):
    w = red(v, p, r)
    if not any(w): return None
    k = next(c for c in w if c); inv = pow(k, p - 2, p); return tuple(c * inv % p for c in w)
base = list({key(v, P1, r1): v for v in base}.values())
# candidate new lines: through two intersection points of the base
ipts = list({key(cross(a, b), P1, r1): cross(a, b) for a, b in itertools.combinations(base, 2) if key(cross(a, b), P1, r1)}.values())
bk = {key(v, P1, r1) for v in base}
cand = {}
for X, Y in itertools.combinations(ipts, 2):
    l = cross(X, Y); k = key(l, P1, r1)
    if k and k not in bk: cand[k] = l
cand = list(cand.values())
print(name, base_name, 'base', len(base), 'intersection points', len(ipts), 'candidate lines', len(cand), flush=True)
def lat(V, p, r):
    N = np.array([red(v, p, r) for v in V], dtype=np.int64)
    return N, Lattice(N, p)
def cert(L, n):
    F = L.flats[2]; adj = {n + k: list(S) for k, S in enumerate(F)}
    return pynauty.certificate(pynauty.Graph(n + len(F), adjacency_dict=adj, vertex_coloring=[set(range(n)), set(range(n, n + len(F)))]))
rng = random.Random(seed); seen = set(); t0 = time.time(); st = {'tried': 0, 'chi': 0, 'free': 0, 'nondf_free': 0}
out = open('/Users/agent/Projects/TeraoHunter/work/numfield/deep_%s_%s_%d.jsonl' % (name, base_name, seed), 'w')
while time.time() - t0 < secs:
    st['tried'] += 1
    D = rng.sample(range(len(base)), rng.randint(0, min(4, len(base) - 4)))
    Ad = rng.sample(range(len(cand)), rng.randint(0, min(3, len(cand))))
    V = [v for i, v in enumerate(base) if i not in D] + [cand[j] for j in Ad]
    if len(V) < 10: continue
    try:
        N, L = lat(V, P1, r1)
    except ValueError:
        continue
    ex = L.exponents()
    if ex is None: continue
    st['chi'] += 1
    if freeness(N, P1, L)['status'] != 'free': continue
    st['free'] += 1
    if rank3_df(N, P1, L): continue
    try:
        N2, L2 = lat(V, P2, r2)
    except ValueError:
        continue
    if sorted(len(S) for S in L2.flats[2]) != sorted(len(S) for S in L.flats[2]): continue
    c = cert(L, len(N))
    if c in seen: continue
    seen.add(c); st['nondf_free'] += 1
    md, _ = moduli_dimension(N, P1, L)
    rec = {'field': name, 'base': base_name, 'kind': 'lines', 'data': V, 'n': len(V), 'exps': list(ex),
           't': sorted(len(S) for S in L.flats[2] if len(S) > 2), 'moduli': int(md)}
    out.write(json.dumps(rec) + '\n'); out.flush()
print(name, base_name, st, 'distinct free non-DF lattices:', len(seen), flush=True)
