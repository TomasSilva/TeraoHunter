# Rank-4 analogue of connecting lines: A(P) = all planes of PG(3,q) through >= 3 non-collinear
# points of a point set P, for ALL k-point sets P up to PGL(4,q) (orderly generation, pynauty).
# Localization of A(P) at a point Q of P = connecting-line arrangement of the projection of
# P \ {Q} from Q.  Filters: chi splits -> locally free -> non-DF -> free -> all realizations.
import sys, itertools, time, json
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
import pynauty
from hunter.core import Lattice, localization_l4, rank_mod
from hunter.ffree import free_rank3_ziegler, free_simple
from hunter.divisional import df_witnesses
from hunter.realize import Realizer
from hunter.moduli import moduli_dimension

q = int(sys.argv[1]); KMIN = int(sys.argv[2]); KMAX = int(sys.argv[3])
pts = [v for v in itertools.product(range(q), repeat=4) if any(v) and next(c for c in v if c) == 1]
NP = len(pts)
planes = pts
P_arr = np.array(pts, dtype=np.int64)
inc_m = (P_arr @ P_arr.T) % q == 0                      # inc_m[plane, point]
inc = [list(np.nonzero(inc_m[i])[0]) for i in range(NP)]
adj = {NP + i: inc[i] for i in range(NP)}
print('PG(3,%d): %d points/planes' % (q, NP), flush=True)


def cert(P):
    P = set(P)
    g = pynauty.Graph(2 * NP, directed=False, adjacency_dict=adj,
                      vertex_coloring=[set(P), set(j for j in range(NP) if j not in P), set(range(NP, 2 * NP))])
    return pynauty.certificate(g)


def locally_free(N, L):
    for S in L.flats[3]:
        if len(S) <= 3:
            continue
        if free_rank3_ziegler(localization_l4(N, S, q), q) != 'free':
            return False
    return True


levels = {0: {cert([]): ()}}
out = open(f'/Users/agent/Projects/TeraoHunter/work/spanned_{q}.jsonl', 'w')
t0 = time.time()
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
    st = {'orbits': len(nxt), 'chi': 0, 'locfree': 0, 'nondf': 0, 'free': 0, 'mixed': 0}
    if k >= KMIN:
        for P in nxt.values():
            counts = inc_m[:, list(P)].sum(axis=1)
            # planes through >= 3 points of P that are not all collinear: a plane through 3
            # collinear points also contains the line; require the points in it to span a plane
            cand = np.nonzero(counts >= 3)[0]
            Ls = []
            for i in cand:
                pp = [P_arr[j] for j in P if inc_m[i, j]]
                if rank_mod(np.array(pp, dtype=np.int64) % q, q) >= 3:
                    Ls.append(i)
            if len(Ls) < 5:
                continue
            N = P_arr[Ls] % q
            if rank_mod(N, q) < 4:
                continue
            try:
                L = Lattice(N, q)
            except ValueError:
                continue
            if L.exponents() is None:
                continue
            st['chi'] += 1
            if not locally_free(N, L):
                continue
            st['locfree'] += 1
            if df_witnesses(N, q, L):
                continue
            st['nondf'] += 1
            fr, g = free_simple(N, q, L)
            rec = {'q': q, 'k': k, 'P': [int(x) for x in P], 'n_planes': len(Ls), 'exps': L.exponents(), 'status': fr}
            if fr == 'free':
                st['free'] += 1
                try:
                    sols = Realizer(N, q, max_solutions=2000).solutions()
                except ValueError:
                    sols = []
                frs = [free_simple(s, q)[0] for s in sols]
                nf = len(frs) - frs.count('free')
                rec.update({'realizations': len(sols), 'nonfree_realizations': nf})
                print('FREE NON-DF', json.dumps(rec), flush=True)
                if nf and nf < len(frs):
                    st['mixed'] += 1
                    print('!!! MIXED FREENESS (TC COUNTEREXAMPLE OVER GF(%d), P^3)' % q, flush=True)
            out.write(json.dumps(rec) + '\n'); out.flush()
    print('k=%d %s  (%.0fs)' % (k, st, time.time() - t0), flush=True)
