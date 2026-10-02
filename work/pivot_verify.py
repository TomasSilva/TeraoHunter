# Verify predictor candidates: freeness at random s (big p), then scan all s over a small prime.
import sys, ast, re, time, random
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import multi_freeness, Lattice, module_generators
from hunter.multi import euler_restriction_dominant
rng = random.Random(11)
pat = re.compile(r"^CAND (\[\[.*\]\]) X1 (\[[^\]]*\]) m (\([^)]*\)) L0mults (\[[^\]]*\])")
P_BIG, P_SMALL = 32003, int(sys.argv[2]) if len(sys.argv) > 2 else 1009
done = set()
for line in open(sys.argv[1]):
    mm = pat.match(line)
    if not mm: continue
    N = ast.literal_eval(mm.group(1)); X1 = ast.literal_eval(mm.group(2)); m = ast.literal_eval(mm.group(3))
    Ls = ast.literal_eval(mm.group(4))
    key = (str(N), str(X1), m)
    if key in done: continue
    done.add(key)
    u, w = N[X1[0]], N[X1[1]]
    M = sum(Ls)
    m0 = M + 1
    def B(s, p):
        return np.array(N + [[(u[i] + s * w[i]) % p for i in range(3)]], dtype=np.int64) % p
    mult = list(m) + [m0]
    res = [multi_freeness(B(rng.randrange(2, P_BIG), P_BIG), mult, P_BIG) for _ in range(2)]
    if not any(r[0] == 'free' for r in res):
        print('PREDICTION FAILED', N, X1, m, res, flush=True)
        continue
    p = P_SMALL
    sig = Lattice(B(7, p), p).signature()
    bad, nsame, jumps = [], 0, []
    t0 = time.time()
    for s in range(1, p):
        try:
            L = Lattice(B(s, p), p)
        except ValueError:
            continue
        if L.signature() != sig:
            continue
        nsame += 1
        st, g = multi_freeness(B(s, p), mult, p)
        if st != 'free':
            nor, ms = euler_restriction_dominant(B(s, p), mult, len(N), p)
            eg, _ = module_generators(nor, ms, p)
            bad.append((s, st, tuple(g), 'euler', tuple(eg[:2])))
    print('FREE-GENERIC', N, 'X1', X1, 'm', m, 'm0', m0, 'L0mults', Ls, 'gens', res[0][1],
          '| same-lattice s', nsame, '| NONFREE', bad[:8], len(bad), '(%.0fs)' % (time.time() - t0), flush=True)
