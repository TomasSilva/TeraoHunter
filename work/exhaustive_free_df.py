# Is every free arrangement over F_q divisionally free?  Exhaustive over all hyperplane subsets of
# PG(l-1, q) up to PGL (orderly generation with pynauty).  For every essential subset: DF (Abe's
# recursive definition, purely combinatorial); for non-DF ones with split chi: exact freeness.
import sys, itertools, time, json
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
import pynauty
from hunter.core import Lattice, rank_mod, integer_roots
from hunter.lattice_gen import LatticeGen
from hunter.restrict import ziegler_restriction
from hunter.divisional import poly_div
from hunter.ffree import ExactModule

q = int(sys.argv[1]); l = int(sys.argv[2]); NMAX_FREE = int(sys.argv[3]) if len(sys.argv) > 3 else 40
pts = [v for v in itertools.product(range(q), repeat=l) if any(v) and next(c for c in v if c) == 1]
NP = len(pts)
P_arr = np.array(pts, dtype=np.int64)
inc = (P_arr @ P_arr.T) % q == 0
adj = {NP + i: list(np.nonzero(inc[i])[0]) for i in range(NP)}


def cert(S):
    S = set(S)
    g = pynauty.Graph(2 * NP, directed=False, adjacency_dict=adj,
                      vertex_coloring=[set(range(NP)), set(NP + i for i in S), set(NP + i for i in range(NP) if i not in S)])
    return pynauty.certificate(g)


def chi_of(N):
    ll = N.shape[1]
    if ll == 2:
        k = N.shape[0]
        return [1, -k, k - 1]
    if ll <= 4:
        return Lattice(N, q).char_poly()
    return LatticeGen(N, q).char_poly()


_memo = {}


def is_df(N):
    ll = N.shape[1]
    if ll <= 2:
        return True
    key = (ll, tuple(sorted(map(tuple, N.tolist()))))
    if key in _memo:
        return _memo[key]
    chi = chi_of(N)
    res = False
    if integer_roots(chi) is not None:
        for i in range(N.shape[0]):
            R, mu, groups = ziegler_restriction(N, i, q)
            if rank_mod(R, q) < ll - 1:
                continue
            if poly_div(chi, chi_of(R)) is None:
                continue
            if is_df(R):
                res = True
                break
    _memo[key] = res
    return res


levels = {0: {cert([]): ()}}
t0 = time.time()
summary = {'orbits': 0, 'essential': 0, 'chi_split': 0, 'df': 0, 'nondf_split': 0, 'nondf_free': 0, 'skipped_big': 0}
found = []
nondf_list = []
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
    summary['orbits'] += len(nxt)
    for S in nxt.values():
        N = P_arr[list(S)] % q
        if k < l or rank_mod(N, q) < l:
            continue
        summary['essential'] += 1
        chi = chi_of(N)
        if integer_roots(chi) is None:
            continue
        summary['chi_split'] += 1
        if is_df(N):
            summary['df'] += 1
            continue
        summary['nondf_split'] += 1
        nondf_list.append(N.tolist())
        if k > NMAX_FREE:
            summary['skipped_big'] += 1
            continue
        st, g = ExactModule(N, [1] * k, q).freeness()
        if st == 'free':
            summary['nondf_free'] += 1
            found.append({'k': k, 'S': list(S), 'exps': integer_roots(chi), 'normals': N.tolist()})
            print('FREE BUT NOT DF over F_%d, rank %d:' % (q, l), json.dumps(found[-1]), flush=True)
    print('k=%d orbits=%d  %s  (%.0fs)' % (k, len(nxt), summary, time.time() - t0), flush=True)
json.dump({'q': q, 'l': l, 'summary': summary, 'found': found, 'nondf_split': nondf_list},
          open(f'/Users/agent/Projects/TeraoHunter/work/free_df_{q}_{l}.json', 'w'))
print('FINAL', summary, flush=True)
