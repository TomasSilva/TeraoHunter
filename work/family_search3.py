# 2-parameter families restricted to random lines; ALL multiplicity vectors in {1..K}^n filtered by
# the necessary condition b2(B,m) = sigma2(d) for integer exponents d (sum |m|), then freeness tests.
import sys, random, time, json, itertools
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice, multi_freeness, _cross3
from hunter.multi import wakamiko, local_exponents
p = 1009
seed = int(sys.argv[1]); rng = random.Random(seed)
NPAR = int(sys.argv[3]) if len(sys.argv) > 3 else 2
out = open(f'/Users/agent/Projects/TeraoHunter/work/family3_{seed}.jsonl', 'w')
def norm(v):
    v = [int(x) % p for x in v]
    for c in v:
        if c:
            inv = pow(c, p - 2, p)
            return tuple(x * inv % p for x in v)
    return None
def replay2(recipe, line, u):
    par = [(line[2*w] + line[2*w+1] * u) % p for w in range(NPAR)]
    lines = [norm(v) for v in ([1,0,0],[0,1,0],[0,0,1],[1,1,1])]
    for step in recipe:
        if step[0] == 'pp':
            (i, j), (k, l) = step[1], step[2]
            P = _cross3(np.array(lines[i]), np.array(lines[j]), p)
            Q = _cross3(np.array(lines[k]), np.array(lines[l]), p)
            v = norm(_cross3(P, Q, p))
        else:
            (i, j), w = step[1], step[2]
            v = norm([(lines[i][k] + par[w] * lines[j][k]) % p for k in range(3)])
        if v is None or v in lines:
            return None
        lines.append(v)
    return np.array(lines, dtype=np.int64)
def random_recipe(n):
    while True:
        recipe, k = [], 4
        pos = sorted(rng.sample(range(n - 4), NPAR))
        for s in range(n - 4):
            pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
            if s in pos:
                recipe.append(('p1', rng.choice(pairs), pos.index(s)))
            else:
                a, b = rng.sample(pairs, 2)
                recipe.append(('pp', a, b))
            k += 1
        line = [rng.randrange(1, p) for _ in range(2 * NPAR)]
        B1, B2 = replay2(recipe, line, rng.randrange(p)), replay2(recipe, line, rng.randrange(p))
        if B1 is None or B2 is None: continue
        try:
            L1, L2 = Lattice(B1, p), Lattice(B2, p)
        except ValueError:
            continue
        if L1.signature() != L2.signature(): continue
        if sum(1 for S in L1.flats[2] if len(S) >= 3) < 3: continue
        return recipe, line, L1, B1
def exps_for(M, b2):
    sols = []
    for d1 in range(0, M // 3 + 1):
        for d2 in range(d1, (M - d1) // 2 + 1):
            d3 = M - d1 - d2
            if d1 * d2 + d1 * d3 + d2 * d3 == b2:
                sols.append((d1, d2, d3))
    return sols
K = 4
t0 = time.time(); nfam = 0; ntest = 0
while time.time() - t0 < float(sys.argv[2]):
    n = rng.randint(7, 9)
    recipe, line, L, B1 = random_recipe(n)
    nfam += 1
    sig = L.signature()
    pts = [sorted(S) for S in L.flats[2]]
    big = [S for S in pts if len(S) >= 4]
    if big:
        continue   # keep to <= 3-line points: local exponents combinatorial (Wakamiko)
    for m in itertools.product(range(1, K + 1), repeat=n):
        M = sum(m)
        b2 = 0
        for S in pts:
            if len(S) == 2:
                b2 += m[S[0]] * m[S[1]]
            else:
                e = wakamiko([m[i] for i in S]); b2 += e[0] * e[1]
        cand = exps_for(M, b2)
        if not cand: continue
        ntest += 1
        ua, ub = rng.randrange(p), rng.randrange(p)
        Ba, Bb = replay2(recipe, line, ua), replay2(recipe, line, ub)
        if Ba is None or Bb is None: continue
        if multi_freeness(Ba, m, p)[0] != 'free': continue
        if multi_freeness(Bb, m, p)[0] != 'free': continue
        bad = []
        for u in range(p):
            B = replay2(recipe, line, u)
            if B is None: continue
            try:
                Lt = Lattice(B, p)
            except ValueError:
                continue
            if Lt.signature() != sig: continue
            st, g = multi_freeness(B, m, p)
            if st != 'free': bad.append((u, list(g)))
        rec = {'recipe': recipe, 'line': line, 'npar': NPAR, 'n': n, 'm': list(m), 'cand_exps': cand, 'nonfree_u': bad,
               'mult_points': [S for S in pts if len(S) >= 3]}
        out.write(json.dumps(rec) + '\n'); out.flush()
        if bad: print('HIT', json.dumps(rec), flush=True)
print('families', nfam, 'tests', ntest, 'time', time.time() - t0)
