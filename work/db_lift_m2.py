# Multiplicities with <= 2 heavy lines per multiple point, heavy-pair graph a forest.
import sys, json, random, itertools, time
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.db import records, lines_mod_p
from hunter.core import Lattice, multi_freeness
p = 1000003
info = {r['hash']: r for r in json.load(open('/Users/agent/Projects/TeraoHunter/work/db_moduli.json'))}
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
t0 = time.time()
stats = {'tested': 0, 'free': 0}
def is_forest(n, edges):
    par = list(range(n))
    def f(a):
        while par[a] != a:
            par[a] = par[par[a]]; a = par[a]
        return a
    for a, b in edges:
        ra, rb = f(a), f(b)
        if ra == rb: return False
        par[ra] = rb
    return True
for rec in records():
    r = info.get(rec['lattice_hash'])
    if not r or r['moduli_dim'] < 1 or rec['n'] > 16:
        continue
    N = np.array(lines_mod_p(rec, p), dtype=np.int64)
    L = Lattice(N, p)
    n = N.shape[0]
    mp = [sorted(S) for S in L.flats[2] if len(S) >= 3]
    tried = set()
    for _ in range(40):
        k = rng.randint(2, min(6, n))
        J = sorted(rng.sample(range(n), k))
        ok = True
        edges = []
        for S in mp:
            h = [i for i in S if i in J]
            if len(h) > 2: ok = False; break
            if len(h) == 2: edges.append(tuple(h))
        if not ok or not is_forest(n, edges): continue
        cs = tuple(rng.randint(2, 3) for _ in J)
        key = tuple(zip(J, cs))
        if key in tried: continue
        tried.add(key)
        m = [1] * n
        for i, c in key: m[i] = c
        st, g = multi_freeness(N, m, p)
        stats['tested'] += 1
        if st == 'free':
            stats['free'] += 1
            print('FREE', rec['lattice_hash'], 'n', rec['n'], 'moduli', r['moduli_dim'], 'heavy', key, 'pairs', edges, 'exps', g, flush=True)
print(stats, 'time %.0f' % (time.time() - t0))
