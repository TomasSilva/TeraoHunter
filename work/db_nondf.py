import sys, json
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.db import records, lines_mod_p
from hunter.core import Lattice
from hunter.divisional import rank3_df
p = 1000003
info = {r['hash']: r for r in json.load(open('/Users/agent/Projects/TeraoHunter/work/db_moduli.json'))}
rows = []
for rec in records():
    N = lines_mod_p(rec, p)
    if N is None: continue
    N = np.array(N, dtype=np.int64)
    L = Lattice(N, p)
    df = rank3_df(N, p, L)
    r = info.get(rec['lattice_hash'], {})
    rows.append((rec['n'], df, r.get('moduli_dim', -1), rec['lattice_hash'], tuple(rec['exponents']), rec['field']['name'] if isinstance(rec['field'], dict) else 'QQ'))
from collections import Counter
print('total', len(rows), 'non-DF', sum(1 for r in rows if not r[1]))
print('non-DF by n:', sorted(Counter(r[0] for r in rows if not r[1]).items()))
print('non-DF with moduli>0:', [(r[0], r[2], r[4], r[5], r[3]) for r in rows if not r[1] and r[2] > 0])
json.dump([r for r in rows if not r[1]], open('/Users/agent/Projects/TeraoHunter/work/db_nondf.json', 'w'))
