# Exhaustive "free => DF ?" over PG(l-1, q), parallel evaluation.
#   phase 1: orbit representatives of hyperplane subsets (pynauty, orderly generation)
#   phase 2: for each essential representative: chi, DF (Abe), exact freeness if non-DF & split chi
import sys, itertools, time, json, pickle, os
from multiprocessing import Pool
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np

q = int(sys.argv[1]); l = int(sys.argv[2]); W = int(sys.argv[3])
pts = [v for v in itertools.product(range(q), repeat=l) if any(v) and next(c for c in v if c) == 1]
NP = len(pts)
P_arr = np.array(pts, dtype=np.int64)
OUT = f'/Users/agent/Projects/TeraoHunter/work/free_df_par_{q}_{l}'


def phase1():
    import pynauty
    inc = (P_arr @ P_arr.T) % q == 0
    adj = {NP + i: list(np.nonzero(inc[i])[0]) for i in range(NP)}

    def cert(S):
        S = set(S)
        g = pynauty.Graph(2 * NP, directed=False, adjacency_dict=adj,
                          vertex_coloring=[set(range(NP)), set(NP + i for i in S), set(NP + i for i in range(NP) if i not in S)])
        return pynauty.certificate(g)
    levels = {0: {cert([]): ()}}
    reps = []
    t0 = time.time()
    for k in range(1, NP + 1):
        nxt = {}
        for S in levels[k - 1].values():
            for j in range(NP):
                if j in S:
                    continue
                T = tuple(sorted(S + (j,)))
                c = cert(T)
                if c not in nxt:
                    nxt[c] = T
        levels[k] = nxt
        levels[k - 1] = None
        reps.extend(nxt.values())
        print('k=%d orbits=%d total=%d (%.0fs)' % (k, len(nxt), len(reps), time.time() - t0), flush=True)
    pickle.dump(reps, open(OUT + '_reps.pkl', 'wb'))
    return reps


def evaluate(chunk):
    from hunter.core import Lattice, rank_mod, integer_roots
    from hunter.lattice_gen import LatticeGen
    from hunter.restrict import ziegler_restriction
    from hunter.divisional import poly_div
    from hunter.ffree import ExactModule
    memo = {}

    def chi_of(N):
        ll = N.shape[1]
        if ll == 2:
            return [1, -N.shape[0], N.shape[0] - 1]
        return (Lattice(N, q) if ll <= 4 else LatticeGen(N, q)).char_poly()

    def is_df(N):
        ll = N.shape[1]
        if ll <= 2:
            return True
        key = (ll, tuple(sorted(map(tuple, N.tolist()))))
        if key in memo:
            return memo[key]
        chi = chi_of(N)
        res = False
        if integer_roots(chi) is not None:
            for i in range(N.shape[0]):
                R, mu, g = ziegler_restriction(N, i, q)
                if rank_mod(R, q) < ll - 1:
                    continue
                if poly_div(chi, chi_of(R)) is None:
                    continue
                if is_df(R):
                    res = True
                    break
        memo[key] = res
        return res

    st = {'essential': 0, 'chi_split': 0, 'df': 0, 'nondf_split': 0, 'nondf_free': 0}
    nondf, found = [], []
    for S in chunk:
        N = P_arr[list(S)] % q
        if len(S) < l or rank_mod(N, q) < l:
            continue
        st['essential'] += 1
        chi = chi_of(N)
        if integer_roots(chi) is None:
            continue
        st['chi_split'] += 1
        if is_df(N):
            st['df'] += 1
            continue
        st['nondf_split'] += 1
        nondf.append(N.tolist())
        fr, g = ExactModule(N, [1] * len(S), q).freeness()
        if fr == 'free':
            st['nondf_free'] += 1
            found.append(N.tolist())
    return st, nondf, found


if __name__ == '__main__':
    if os.path.exists(OUT + '_reps.pkl'):
        reps = pickle.load(open(OUT + '_reps.pkl', 'rb'))
    else:
        reps = phase1()
    chunks = [reps[i::W * 8] for i in range(W * 8)]
    t0 = time.time()
    total = {'essential': 0, 'chi_split': 0, 'df': 0, 'nondf_split': 0, 'nondf_free': 0}
    allnondf, allfound = [], []
    with Pool(W) as pool:
        for i, (st, nondf, found) in enumerate(pool.imap_unordered(evaluate, chunks)):
            for k2 in total:
                total[k2] += st[k2]
            allnondf += nondf
            allfound += found
            print('chunk %d/%d %s (%.0fs)' % (i + 1, len(chunks), total, time.time() - t0), flush=True)
    json.dump({'q': q, 'l': l, 'summary': total, 'nondf_split': allnondf, 'free_nondf': allfound},
              open(OUT + '.json', 'w'))
    print('FINAL', total, flush=True)
