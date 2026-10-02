import sys, time, json
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.db import records, lines_mod_p, field_d
from hunter.core import Lattice
from hunter.moduli import moduli_dimension
p = 1000003
out = []
t0 = time.time()
cnt = 0
for rec in records():
    cnt += 1
    N = lines_mod_p(rec, p)
    if N is None:
        continue
    N = np.array(N, dtype=np.int64)
    try:
        L = Lattice(N, p)
    except ValueError:
        continue
    # sanity: same t-vector as recorded
    tv = {}
    for S in L.flats[2]:
        if len(S) >= 2: tv[len(S)] = tv.get(len(S), 0) + 1
    rec_tv = {int(k): v for k, v in rec['t_vector'].items()}
    tv = {k: v for k, v in tv.items() if k >= 3 or k in rec_tv}
    ok = all(tv.get(k, 0) == v for k, v in rec_tv.items() if k >= 3)
    dim, _ = moduli_dimension(N, p, L)
    out.append({'hash': rec['lattice_hash'], 'n': rec['n'], 'exps': rec['exponents'], 'field': field_d(rec),
                'moduli_dim': int(dim), 'tv_ok': ok, 'supersolvable': rec['supersolvable']})
print('processed', cnt, 'reduced', len(out), 'time %.0f' % (time.time() - t0))
json.dump(out, open('/Users/agent/Projects/TeraoHunter/work/db_moduli.json', 'w'))
from collections import Counter
print('moduli dim histogram', sorted(Counter(r['moduli_dim'] for r in out).items()))
print('tv mismatches', sum(1 for r in out if not r['tv_ok']))
pos = [r for r in out if r['moduli_dim'] > 0]
print('positive moduli by n', sorted(Counter(r['n'] for r in pos).items()))
print('positive moduli non-supersolvable by n', sorted(Counter(r['n'] for r in pos if not r['supersolvable']).items()))
