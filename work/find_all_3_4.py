# Full exhaustive PG(3,3) free-vs-DF check, same chunks as free_df_parallel.py 3 4 8 (reps[i::64]).
# Exact local freeness (ExactModule on every rank-3 localisation) is tested before the exact
# freeness test of the whole arrangement; local freeness is necessary for freeness over any field.
import sys, time, json, pickle
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice, rank_mod, integer_roots, localization_l4
from hunter.restrict import ziegler_restriction
from hunter.divisional import poly_div
from hunter.ffree import ExactModule, free_rank3_ziegler
import itertools

q, l = 3, 4
pts = [v for v in itertools.product(range(q), repeat=l) if any(v) and next(c for c in v if c) == 1]
P_arr = np.array(pts, dtype=np.int64)
reps = pickle.load(open('/Users/agent/Projects/TeraoHunter/work/free_df_par_3_4_reps.pkl', 'rb'))
order = [int(x) for x in sys.argv[1].split(',')]
tag = sys.argv[2]
memo = {}


def chi_of(N):
    if N.shape[1] == 2:
        return [1, -N.shape[0], N.shape[0] - 1]
    return Lattice(N, q).char_poly()


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
            if rank_mod(R, q) < ll - 1 or poly_div(chi, chi_of(R)) is None:
                continue
            if is_df(R):
                res = True
                break
    memo[key] = res
    return res


def locally_free_exact(N):
    """Exact local freeness: every rank-3 localisation free (ExactModule), cross-checked
    against the Yoshinaga-Ziegler shortcut; disagreements are logged."""
    L = Lattice(N, q)
    ok = True
    for S in L.flats[3]:
        if len(S) <= 3:
            continue
        loc = localization_l4(N, S, q)
        a = ExactModule(loc, [1] * len(S), q).freeness()[0] == 'free'
        b = free_rank3_ziegler(loc, q) == 'free'
        if a != b:
            print('DISAGREE', json.dumps(loc.tolist()), a, b, file=out, flush=True)
        if not a:
            ok = False
            break
    return ok


out = open(f'/Users/agent/Projects/TeraoHunter/work/find_all_3_4_{tag}.out', 'a')
for ci in order:
    t0 = time.time()
    st = {'essential': 0, 'chi_split': 0, 'df': 0, 'nondf_split': 0, 'loc_free': 0, 'free': 0}
    tdf = tloc = tfree = 0.0
    for S in reps[ci::64]:
        N = P_arr[list(S)] % q
        if len(S) < l or rank_mod(N, q) < l:
            continue
        st['essential'] += 1
        t1 = time.time()
        chi = chi_of(N)
        if integer_roots(chi) is None:
            tdf += time.time() - t1
            continue
        st['chi_split'] += 1
        d = is_df(N)
        tdf += time.time() - t1
        if d:
            st['df'] += 1
            continue
        st['nondf_split'] += 1
        t1 = time.time()
        lf = locally_free_exact(N)
        tloc += time.time() - t1
        if not lf:
            continue
        st['loc_free'] += 1
        t1 = time.time()
        fr, g = ExactModule(N, [1] * len(S), q).freeness()
        tfree += time.time() - t1
        if fr == 'free':
            st['free'] += 1
            rec = {'chunk': ci, 'S': list(S), 'k': len(S), 'chi': chi, 'exps': integer_roots(chi),
                   'gens': g, 'normals': N.tolist()}
            print('HIT', json.dumps(rec), file=out, flush=True)
    print('chunk', ci, st, 'df %.0fs loc %.0fs free %.0fs total %.0fs' % (tdf, tloc, tfree, time.time() - t0),
          file=out, flush=True)
print('DONE', file=out, flush=True)
