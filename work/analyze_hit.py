import sys, json
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice, multi_freeness
from hunter.multi import local_exponents
sys.argv = sys.argv
p = 1009
exec(open('/Users/agent/Projects/TeraoHunter/work/family_search.py').read().split("t0 = time.time()")[0].split("seed = int")[0])
from hunter.core import _cross3
def norm(v):
    v = [int(x) % p for x in v]
    for c in v:
        if c:
            inv = pow(c, p - 2, p)
            return tuple(x * inv % p for x in v)
def replay(recipe, t):
    lines = [norm(v) for v in ([1,0,0],[0,1,0],[0,0,1],[1,1,1])]
    for step in recipe:
        if step[0] == 'pp':
            (a, b), (c, d) = step[1], step[2]
            P = _cross3(np.array(lines[a]), np.array(lines[b]), p)
            Q = _cross3(np.array(lines[c]), np.array(lines[d]), p)
            v = norm(_cross3(P, Q, p))
        else:
            (a, b) = step[1]
            v = norm([(lines[a][k] + t * lines[b][k]) % p for k in range(3)])
        if v is None or v in lines: return None
        lines.append(v)
    return np.array(lines, dtype=np.int64)
rec = json.loads(sys.stdin.read())
recipe, m = rec['recipe'], rec['m']
tg = 500
for t in [tg, 77] + [x[0] for x in rec['nonfree_t']]:
    B = replay(recipe, t)
    L = Lattice(B, p)
    loc = {tuple(sorted(S)): local_exponents(B, m, S, p) for S in L.flats[2] if len(S) >= 3}
    print('t', t, multi_freeness(B, m, p), 'local exps:', loc)
