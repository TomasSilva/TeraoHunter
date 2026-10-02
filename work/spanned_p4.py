# P^4 (rank 5): A(P) = hyperplanes of PG(4,q) spanned by >= 4 points of P (spanning a solid),
# for ALL k-point sets P up to PGL(5,q).  Filters: chi splits -> not DF (rank-5 DF via rank-4
# restrictions) -> free (exact over F_q).  Free non-DF ones are reported.
import sys, itertools, time, json
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
import pynauty
from hunter.core import rank_mod, Lattice
from hunter.lattice_gen import LatticeGen
from hunter.restrict import ziegler_restriction
from hunter.divisional import df_witnesses, poly_div
from hunter.ffree import ExactModule

q = int(sys.argv[1]); KMIN = int(sys.argv[2]); KMAX = int(sys.argv[3])
pts = [v for v in itertools.product(range(q), repeat=5) if any(v) and next(c for c in v if c) == 1]
NP = len(pts)
P_arr = np.array(pts, dtype=np.int64)
inc_m = (P_arr @ P_arr.T) % q == 0
adj = {NP + i: list(np.nonzero(inc_m[i])[0]) for i in range(NP)}
print('PG(4,%d): %d points/hyperplanes' % (q, NP), flush=True)


def cert(P):
    P = set(P)
    g = pynauty.Graph(2 * NP, directed=False, adjacency_dict=adj,
                      vertex_coloring=[set(P), set(j for j in range(NP) if j not in P), set(range(NP, 2 * NP))])
    return pynauty.certificate(g)


def df_rank5(N, chi):
    for i in range(N.shape[0]):
        R, mu, groups = ziegler_restriction(N, i, q)
        if R.shape[0] < 4 or rank_mod(R, q) < 4:
            continue
        try:
            LR = Lattice(R, q)
        except ValueError:
            continue
        if poly_div(chi, LR.char_poly()) is None:
            continue
        if df_witnesses(R, q, LR):
            return True
    return False


levels = {0: {cert([]): ()}}
t0 = time.time()
out = open(f'/Users/agent/Projects/TeraoHunter/work/spanned_p4_{q}.jsonl', 'w')
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
    st = {'orbits': len(nxt), 'chi': 0, 'nondf': 0, 'free': 0}
    if k >= KMIN:
        for P in nxt.values():
            counts = inc_m[:, list(P)].sum(axis=1)
            Ls = [i for i in np.nonzero(counts >= 4)[0]
                  if rank_mod(P_arr[[j for j in P if inc_m[i, j]]] % q, q) >= 4]
            if len(Ls) < 6:
                continue
            N = P_arr[Ls] % q
            if rank_mod(N, q) < 5:
                continue
            try:
                L = LatticeGen(N, q)
            except ValueError:
                continue
            ex = L.exponents()
            if ex is None:
                continue
            st['chi'] += 1
            if df_rank5(N, L.char_poly()):
                continue
            st['nondf'] += 1
            fr, g = ExactModule(N, [1] * N.shape[0], q).freeness()
            rec = {'q': q, 'k': k, 'P': [int(x) for x in P], 'n': len(Ls), 'exps': ex, 'status': fr, 'gens': g}
            out.write(json.dumps(rec) + '\n'); out.flush()
            if fr == 'free':
                st['free'] += 1
                print('FREE NON-DF in P^4', json.dumps(rec), flush=True)
    print('k=%d %s  (%.0fs)' % (k, st, time.time() - t0), flush=True)
