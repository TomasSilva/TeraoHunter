import sys, json, random
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter'); sys.path.insert(0, '/Users/agent/Projects/TeraoHunter/work')
import numpy as np
from hunter.core import Lattice, multi_freeness, _cross3
from work_lift import try_lift
p = 1009
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
tries = int(sys.argv[1]) if len(sys.argv) > 1 else 100
for t in [500, 77] + [x[0] for x in rec['nonfree_t']]:
    B = replay(recipe, t)
    P, found = try_lift(B, m, p, tries=tries, seed=t)
    print('t', t, 'multi', multi_freeness(B, m, p), 'R_X', P.R, '| lifts:', len(found),
          sorted(set((f[4], len(f[3]), f[2]) for f in found)), flush=True)
