import sys, itertools, time
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter'); sys.path.insert(0, '/Users/agent/Projects/TeraoHunter/work')
import numpy as np
from hunter.core import Lattice
from hunter.ffree import free_simple
from family26 import config, spanned
target = (26, {3: 24, 4: 23})
for p in [int(x) for x in sys.argv[1:]]:
    t0 = time.time(); hits = []
    for a in range(2, p):
        for b in range(2, p):
            for c in range(2, p):
                try:
                    N, sets = spanned(config(a, b, c, p), p)
                    if len(N) != 26: continue
                    L = Lattice(N, p)
                except ValueError:
                    continue
                r2 = sorted(len(S) for S in L.flats[2] if len(S) > 2)
                if r2.count(3) == 24 and r2.count(4) == 23:
                    hits.append((a, b, c, L.exponents()))
    print('p', p, 'hits', len(hits), hits[:6], '(%.0fs)' % (time.time() - t0), flush=True)
