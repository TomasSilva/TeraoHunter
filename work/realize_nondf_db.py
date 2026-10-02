import sys, json, time
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.db import records, lines_mod_p
from hunter.core import Lattice
from hunter.realize import Realizer
from hunter.ffree import free_rank3_ziegler
nondf = json.load(open('/Users/agent/Projects/TeraoHunter/work/db_nondf.json'))
hashes = {r[3]: r for r in nondf}
recs = {rec['lattice_hash']: rec for rec in records() if rec['lattice_hash'] in hashes}
P0 = 1000003
primes = [p for p in range(7, 120) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
for h, rec in recs.items():
    Nref = np.array(lines_mod_p(rec, P0), dtype=np.int64)
    Lref = Lattice(Nref, P0)
    tref = sorted(len(S) for S in Lref.flats[2] if len(S) > 2)
    print('lattice', h[:8], 'n', rec['n'], 'exps', rec['exponents'], 'field', rec['field'] if rec['field'] == 'QQ' else rec['field']['name'], flush=True)
    for p in primes:
        for root in (0, 1):
            N = lines_mod_p(rec, p, root_choice=root)
            if N is None: continue
            N = np.array(N, dtype=np.int64)
            try:
                L = Lattice(N, p)
            except ValueError:
                continue
            if sorted(len(S) for S in L.flats[2] if len(S) > 2) != tref or L.exponents() != tuple(rec['exponents']):
                continue
            t = time.time()
            try:
                sols = Realizer(N, p, max_solutions=20000).solutions()
            except ValueError as e:
                print('   p', p, 'error', e); break
            fr = [free_rank3_ziegler(s, p) for s in sols]
            nf = len(fr) - fr.count('free')
            print('   p %3d: %4d realizations, non-free %d%s (%.1fs)' % (p, len(sols), nf, '  <<<<<< MIXED' if nf and nf < len(fr) else '', time.time() - t), flush=True)
            break
