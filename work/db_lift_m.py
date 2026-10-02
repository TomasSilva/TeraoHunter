# For database arrangements C with moduli >= 1: lift-compatible multiplicities
# (independent set J, no two J-lines in a common multiple point), test (C,m) freeness.
import sys, json, random, itertools, time
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.db import records, lines_mod_p
from hunter.core import Lattice, multi_freeness
p = 1000003
info = {r['hash']: r for r in json.load(open('/Users/agent/Projects/TeraoHunter/work/db_moduli.json'))}
rng = random.Random(4)
t0 = time.time()
stats = {'tested': 0, 'free': 0}
hits = []
for rec in records():
    r = info.get(rec['lattice_hash'])
    if not r or r['moduli_dim'] < 1 or rec['n'] > 16:
        continue
    N = np.array(lines_mod_p(rec, p), dtype=np.int64)
    L = Lattice(N, p)
    n = N.shape[0]
    mp = [S for S in L.flats[2] if len(S) >= 3]
    conflict = [[False]*n for _ in range(n)]
    for S in mp:
        for a in S:
            for b in S:
                if a != b: conflict[a][b] = True
    # random independent sets of size 2..4
    tried = set()
    for _ in range(25):
        order = list(range(n)); rng.shuffle(order)
        J = []
        size = rng.randint(2, 4)
        for i in order:
            if all(not conflict[i][j] for j in J):
                J.append(i)
            if len(J) == size: break
        if len(J) < 2: continue
        cs = tuple(rng.randint(2, 3) for _ in J)
        key = (tuple(sorted(zip(J, cs))))
        if key in tried: continue
        tried.add(key)
        m = [1] * n
        for i, c in zip(J, cs): m[i] = c
        st, g = multi_freeness(N, m, p)
        stats['tested'] += 1
        if st == 'free':
            stats['free'] += 1
            hits.append((rec['lattice_hash'], rec['n'], r['moduli_dim'], key, g))
            print('FREE', rec['lattice_hash'], 'n', rec['n'], 'moduli', r['moduli_dim'], 'J,c', key, 'exps', g, flush=True)
print(stats, 'time %.0f' % (time.time() - t0))
