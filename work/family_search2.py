# 2-parameter families (two 'p1' steps with params s,t), restricted to random lines
#   s = a + b u, t = c + d u.  Find (B_u, m) free for generic u, non-free at special u (same lattice).
import sys, random, time, json
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice, multi_freeness, _cross3
p = 1009
seed = int(sys.argv[1]); rng = random.Random(seed)
out = open(f'/Users/agent/Projects/TeraoHunter/work/family2_{seed}.jsonl', 'w')
def norm(v):
    v = [int(x) % p for x in v]
    for c in v:
        if c:
            inv = pow(c, p - 2, p)
            return tuple(x * inv % p for x in v)
    return None
def replay2(recipe, line, u):
    a, b, c, d = line
    par = [(a + b * u) % p, (c + d * u) % p]
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
        pos = sorted(rng.sample(range(n - 4), 2))
        for s in range(n - 4):
            pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
            if s in pos:
                recipe.append(('p1', rng.choice(pairs), pos.index(s)))
            else:
                a, b = rng.sample(pairs, 2)
                recipe.append(('pp', a, b))
            k += 1
        line = [rng.randrange(1, p) for _ in range(4)]
        B1, B2 = replay2(recipe, line, rng.randrange(p)), replay2(recipe, line, rng.randrange(p))
        if B1 is None or B2 is None: continue
        try:
            L1, L2 = Lattice(B1, p), Lattice(B2, p)
        except ValueError:
            continue
        if L1.signature() != L2.signature(): continue
        if sum(1 for S in L1.flats[2] if len(S) >= 3) < 3: continue
        return recipe, line, L1
t0 = time.time()
while time.time() - t0 < float(sys.argv[2]):
    n = rng.randint(7, 10)
    recipe, line, L = random_recipe(n)
    sig = L.signature()
    for _ in range(25):
        m = [rng.choice([1, 1, 2, 2, 3, 4]) for _ in range(n)]
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
        rec = {'recipe': recipe, 'line': line, 'n': n, 'm': m, 'nonfree_u': bad,
               'mult_points': [sorted(S) for S in L.flats[2] if len(S) >= 3]}
        out.write(json.dumps(rec) + '\n'); out.flush()
        if bad: print('HIT', json.dumps(rec), flush=True)
