# Connecting-line arrangements A(P) = all lines through >= 2 points of a point set P in PG(2,q),
# for ALL k-point sets P up to PGL(3,q) (orderly generation, pynauty).  The known non-DF free
# 15-line lattice with moduli is A(P) for 6 special points.  For free non-DF A(P): enumerate all
# GF(q)-realizations of its lattice and test freeness.
import sys, itertools, time, json
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
import pynauty
from hunter.core import Lattice
from hunter.ffree import free_rank3_ziegler
from hunter.divisional import rank3_df
from hunter.realize import Realizer
from hunter.moduli import moduli_dimension

q = int(sys.argv[1]); KMAX = int(sys.argv[2])
pts = [v for v in itertools.product(range(q), repeat=3) if any(v) and next(c for c in v if c) == 1]
NP = len(pts)
lines = pts
inc = [[j for j, P in enumerate(pts) if sum(a * b for a, b in zip(L, P)) % q == 0] for L in lines]
pt_lines = [[i for i in range(NP) if j in set(inc[i])] for j in range(NP)]
adj = {NP + i: inc[i] for i in range(NP)}


def cert(P):
    P = set(P)
    g = pynauty.Graph(2 * NP, directed=False, adjacency_dict=adj,
                      vertex_coloring=[set(P), set(j for j in range(NP) if j not in P), set(range(NP, 2 * NP))])
    return pynauty.certificate(g)


levels = {0: {cert([]): ()}}
t0 = time.time()
out = open(f'/Users/agent/Projects/TeraoHunter/work/connecting_{q}.jsonl', 'w')
stats = {}
for k in range(1, KMAX + 1):
    nxt = {}
    for P in levels[k - 1].values():
        for j in range(NP):
            if j in P:
                continue
            T = tuple(sorted(P + (j,)))
            c = cert(T)
            if c not in nxt:
                nxt[c] = T
    levels[k] = nxt
    st = {'orbits': len(nxt), 'chi': 0, 'free': 0, 'nondf_free': 0, 'with_moduli': 0, 'mixed': 0}
    for P in nxt.values():
        if k < 4:
            continue
        Pset = set(P)
        Ls = [i for i in range(NP) if sum(1 for j in inc[i] if j in Pset) >= 2]
        N = np.array([lines[i] for i in Ls], dtype=np.int64)
        if len(Ls) < 4 or np.linalg.matrix_rank(N.astype(float)) < 3:
            continue
        try:
            L = Lattice(N, q)
        except ValueError:
            continue
        if L.exponents() is None:
            continue
        st['chi'] += 1
        if free_rank3_ziegler(N, q, L) != 'free':
            continue
        st['free'] += 1
        if rank3_df(N, q, L):
            continue
        st['nondf_free'] += 1
        md, _ = moduli_dimension(N, q, L)
        try:
            sols = Realizer(N, q, max_solutions=20000).solutions()
        except ValueError:
            sols = []
        frs = [free_rank3_ziegler(s, q) for s in sols]
        nf = len(frs) - frs.count('free')
        if md > 0:
            st['with_moduli'] += 1
        rec = {'q': q, 'k': k, 'P': list(P), 'n_lines': len(Ls), 'exps': L.exponents(),
               't': sorted(len(F) for F in L.flats[2] if len(F) > 2), 'moduli_tangent': int(md),
               'realizations': len(sols), 'nonfree_realizations': nf}
        out.write(json.dumps(rec) + '\n'); out.flush()
        if nf and nf < len(frs):
            st['mixed'] += 1
            print('!!! MIXED FREENESS (TC COUNTEREXAMPLE OVER GF(%d)):' % q, json.dumps(rec), flush=True)
    stats[k] = st
    print('k=%d %s  (%.0fs)' % (k, st, time.time() - t0), flush=True)
