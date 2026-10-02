# Random 1-parameter families of line arrangements (construction sequences) with
# multiplicities; find (B_t, m) free for generic t but non-free at special t (same lattice).
import sys, random, time, json
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice, multi_freeness, _cross3
p = 1009
seed = int(sys.argv[1]); rng = random.Random(seed)
out = open(f'/Users/agent/Projects/TeraoHunter/work/family_search_{seed}.jsonl', 'w')

def norm(v):
    v = [int(x) % p for x in v]
    for c in v:
        if c:
            inv = pow(c, p - 2, p)
            return tuple(x * inv % p for x in v)
    return None

def replay(recipe, t):
    lines = [norm(v) for v in ([1,0,0],[0,1,0],[0,0,1],[1,1,1])]
    for step in recipe:
        if step[0] == 'pp':
            (a, b), (c, d) = step[1], step[2]
            P = _cross3(np.array(lines[a]), np.array(lines[b]), p)
            Q = _cross3(np.array(lines[c]), np.array(lines[d]), p)
            v = norm(_cross3(P, Q, p))
        else:
            (a, b) = step[1]
            v = norm([(lines[a][k] + t * lines[b][k]) % p for k in range(3)])
        if v is None or v in lines:
            return None
        lines.append(v)
    return np.array(lines, dtype=np.int64)

def random_recipe(n):
    while True:
        recipe = []
        k = 4
        p1_at = rng.randrange(0, n - 4)
        ok = True
        for s in range(n - 4):
            pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
            if s == p1_at:
                recipe.append(('p1', rng.choice(pairs)))
            else:
                a, b = rng.sample(pairs, 2)
                recipe.append(('pp', a, b))
            k += 1
        B1, B2 = replay(recipe, rng.randrange(2, p)), replay(recipe, rng.randrange(2, p))
        if B1 is None or B2 is None:
            continue
        try:
            L1, L2 = Lattice(B1, p), Lattice(B2, p)
        except ValueError:
            continue
        if L1.signature() != L2.signature():
            continue
        if sum(1 for S in L1.flats[2] if len(S) >= 3) < 2:
            continue
        return recipe, L1

t0 = time.time()
nfam = 0
while time.time() - t0 < float(sys.argv[2]):
    n = rng.randint(6, 10)
    recipe, L = random_recipe(n)
    nfam += 1
    sig = L.signature()
    for _ in range(25):
        m = [rng.choice([1, 1, 2, 2, 3, 4]) for _ in range(n)]
        ta, tb = rng.randrange(2, p), rng.randrange(2, p)
        Ba, Bb = replay(recipe, ta), replay(recipe, tb)
        if Ba is None or Bb is None: continue
        sa, ga = multi_freeness(Ba, m, p)
        if sa != 'free': continue
        sb, gb = multi_freeness(Bb, m, p)
        if sb != 'free': continue
        bad = []
        for t in range(1, p):
            B = replay(recipe, t)
            if B is None: continue
            try:
                Lt = Lattice(B, p)
            except ValueError:
                continue
            if Lt.signature() != sig: continue
            st, g = multi_freeness(B, m, p)
            if st != 'free':
                bad.append((t, list(g)))
        rec = {'recipe': recipe, 'n': n, 'm': m, 'exps': ga, 'nonfree_t': bad,
               'mult_points': [sorted(S) for S in L.flats[2] if len(S) >= 3]}
        out.write(json.dumps(rec) + '\n'); out.flush()
        if bad:
            print('HIT', json.dumps(rec), flush=True)
print('families', nfam, 'time', time.time() - t0)
