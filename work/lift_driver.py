import sys, json, random, time
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice, freeness, local_freeness_l4, multi_freeness
from hunter.lift import LiftProblem

def try_lift(B, m, p, tries=200, seed=0, verbose=False):
    rng = random.Random(seed)
    P = LiftProblem(B, m, p)
    found = []
    for k in range(tries):
        pat = P.greedy(rng)
        if pat is None:
            continue
        M = P.system(pat)
        from hunter.core import rank_mod
        rk = rank_mod(M, p) if M.shape[0] else 0
        c = P.solve(pat, rng)
        A = P.arrangement(c)
        try:
            LA = Lattice(A, p)
        except ValueError:
            continue
        loc = local_freeness_l4(A, p, LA)
        fr = freeness(A, p, LA)
        found.append((pat, rk, len(P.planes) - rk, loc, fr['status'], fr['exps']))
        if verbose:
            print('  lift found: conditions', M.shape[0], 'rank', rk, 'sol dim', len(P.planes) - rk,
                  '| nonfree localizations', len(loc), '| A:', fr['status'], fr['exps'])
    return P, found

if __name__ == '__main__':
    p = 1009
    # sanity: X3 family with (3,3,3,1,1,1)
    def X3(t):
        return np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,t % p]], dtype=np.int64) % p
    # -omega for p=1009: t^2 - t + 1 = 0 -> 375, 635
    for t in [123, 375, 1]:
        print('X3 t =', t, 'multi:', multi_freeness(X3(t), [3,3,3,1,1,1], p))
        P, found = try_lift(X3(t), [3,3,3,1,1,1], p, tries=60, seed=1, verbose=False)
        print('   R_X:', P.R, '| lifts found:', len(found), set((f[4], len(f[3])) for f in found))
