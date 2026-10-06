# Stage 2: for every distinct free non-DF lattice found by numfield_search.py, enumerate ALL
# realisations over two medium split primes q (lattice preserved mod q) and test each exactly.
# A Terao counterexample over K would show mixed freeness for almost all split q.
import sys, json, glob, itertools
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice
from hunter.ffree import free_rank3_ziegler
from hunter.realize import Realizer
FIELDS = {'i': (1, 0, 1), 'sqrt2': (1, 0, -2), 'omega': (1, 1, 1), 'phi': (1, -1, -1),
          'sqrt-2': (1, 0, 2), 'sqrt3': (1, 0, -3), 'sqrt-7': (1, 0, 7), 'sqrt5': (1, 0, -5)}
def roots(f, p): return [x for x in range(p) if (f[0]*x*x + f[1]*x + f[2]) % p == 0]
def isprime(p): return p > 1 and all(p % d for d in range(2, int(p**.5)+1))
def red(c, p, r): return (c[0] + c[1]*r) % p
def normals(rec, p, r):
    if rec['kind'] == 'lines':
        return np.array([[red(c, p, r) for c in l] for l in rec['data']], dtype=np.int64)
    V = [np.array([red(c, p, r) for c in x]) for x in rec['data']]
    L = {}
    for u, w in itertools.combinations(V, 2):
        l = np.array([u[1]*w[2]-u[2]*w[1], u[2]*w[0]-u[0]*w[2], u[0]*w[1]-u[1]*w[0]]) % p
        if not l.any(): return None
        k = next(c for c in l if c); L[tuple(int(c)*pow(int(k), p-2, p) % p for c in l)] = 1
    return np.array(list(L), dtype=np.int64)
sigof = lambda L: tuple(sorted(len(F) for F in L.flats[2]))
recs = [json.loads(l) for f in sorted(glob.glob('/Users/agent/Projects/TeraoHunter/work/numfield/hits_*.jsonl')) for l in open(f)]
print(len(recs), 'lattices to screen', flush=True)
for rec in recs:
    f = FIELDS[rec['field']]
    ref = normals(rec, 30011 if roots(f, 30011) else 30013, (roots(f, 30011) or roots(f, 30013))[0])
    done = []
    for q in range(60, 400):
        if len(done) == 2: break
        if not isprime(q) or not roots(f, q): continue
        for r in roots(f, q)[:1]:
            N = normals(rec, q, r)
            if N is None or len(N) != rec['n']: continue
            try:
                L = Lattice(N, q)
            except ValueError:
                continue
            if sorted(len(F) for F in L.flats[2] if len(F) > 2) != rec['t']: continue
            try:
                sols = Realizer(N, q, max_solutions=5000).solutions()
            except ValueError:
                continue
            fr = [free_rank3_ziegler(s, q) for s in sols]
            done.append((q, len(sols), fr.count('free')))
    mixed = any(0 < fr_ < tot for _, tot, fr_ in done)
    print('%s n=%d exps=%s moduli=%d t=%s  realisations(free/total): %s %s' % (
        rec['field'], rec['n'], rec['exps'], rec['moduli'], {k: rec['t'].count(k) for k in set(rec['t'])},
        ['F_%d: %d/%d' % (q, fr_, tot) for q, tot, fr_ in done], '!!! MIXED' if mixed else ''), flush=True)
