import sys, random, itertools, time, json
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice, freeness, local_freeness_l4
from hunter.divisional import df_witnesses
from hunter.gain import bias_normals
p = 2147483629
seed = int(sys.argv[1]); rng = random.Random(seed)
l = 4
pairs = [(i, j) for i in range(l) for j in range(i + 1, l)]
stats = {'tested': 0, 'chi': 0, 'locfree': 0, 'free': 0, 'nondf': 0}
t0 = time.time(); seen = set()
G = int(sys.argv[3]) if len(sys.argv) > 3 else 2
while time.time() - t0 < float(sys.argv[2]):
    coords = [i for i in range(l) if rng.random() < 0.5]
    edges = []
    for (i, j) in pairs:
        k = rng.choice([0, 1, 1, 2, 2, 3])
        gains = rng.sample(range(-G, G + 1), min(k, 2 * G + 1))
        edges += [(i, j, g) for g in gains]
    if len(coords) + len(edges) < 6: continue
    q = rng.randrange(2, p - 1)
    N = bias_normals(l, coords, edges, q, p)
    try:
        L = Lattice(N, p)
    except ValueError:
        continue
    if np.linalg.matrix_rank(N.astype(float)) < 4: continue
    stats['tested'] += 1
    if L.exponents() is None: continue
    stats['chi'] += 1
    if local_freeness_l4(N, p, L): continue
    stats['locfree'] += 1
    fr = freeness(N, p, L)
    if fr['status'] != 'free': continue
    stats['free'] += 1
    w = df_witnesses(N, p, L)
    if not w:
        stats['nondf'] += 1
        key = (tuple(coords), tuple(sorted(edges)))
        if key not in seen:
            seen.add(key)
            print('NON-DF FREE', json.dumps({'coords': coords, 'edges': edges, 'exps': fr['exps'], 'n': len(N)}), flush=True)
print(stats, time.time() - t0)
