import sys, random, itertools
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice, freeness, local_freeness_l4
from hunter.divisional import df_witnesses
from hunter.gain import bias_normals
p = 2147483629
l = 4
pairs = [(i, j) for i in range(l) for j in range(i + 1, l)]
q = 123456789
for coords in [[], [0], [0, 1], [0, 1, 2], [0, 1, 2, 3]]:
    for S in [[0, 1], [-1, 0, 1], [0, 1, 2], [-1, 1], [1, 2], [0, 2], [-1, 0, 1, 2], [-2, -1, 0, 1, 2]]:
        edges = [(i, j, g) for (i, j) in pairs for g in S]
        N = bias_normals(l, coords, edges, q, p)
        try:
            L = Lattice(N, p)
        except ValueError as e:
            print(coords, S, 'error', e); continue
        loc = local_freeness_l4(N, p, L)
        fr = freeness(N, p, L)
        df = df_witnesses(N, p, L) if fr['status'] == 'free' else None
        print(f'coords={coords} gains={S} n={len(N)} chi={L.char_poly()} exps={L.exponents()} locally_nonfree={len(loc)} -> {fr["status"]}  DF-witnesses={len(df) if df is not None else "-"}', flush=True)
