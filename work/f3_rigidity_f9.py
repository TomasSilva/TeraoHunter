# Rigidity test: realizations over GF(9) (and GF(27)) of the rank-3 restrictions A^H of the
# 22-plane F_3 example.  If every restriction lattice has only the F_3-realization (up to
# PGL and Frobenius), the 22-plane lattice cannot give a Terao counterexample in characteristic 3.
import sys, json
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice
from hunter.restrict import ziegler_restriction
from hunter.fq import GF, realizations

rec = json.load(open('/Users/agent/Projects/TeraoHunter/work/f3_rank4_free_nondf.json'))
seen = set()
for h in rec:
    N = np.array(h['normals']) % 3
    for i in range(len(N)):
        R, mu, g = ziegler_restriction(N, i, 3)
        L = Lattice(R, 3)
        key = (R.shape[0], tuple(sorted(len(S) for S in L.flats[2])))
        if key in seen:
            continue
        seen.add(key)
        flats2 = [sorted(S) for S in L.flats[2]]
        out = []
        for qq in [4, 5, 7, 8, 13, 16, 25]:
            sols = realizations(GF(qq), R.shape[0], flats2, max_solutions=200)
            out.append((qq, len(sols)))
        print('restriction with %d lines, point sizes %s: frame-normalised realizations %s' % (
            R.shape[0], sorted(set(key[1])), out), flush=True)
