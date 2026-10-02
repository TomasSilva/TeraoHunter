import sys, random, itertools
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.ffree import ExactModule, free_simple
from hunter.m2 import run_m2
def m2(normals, mult, p):
    l = len(normals[0]); names = [f"x{i}" for i in range(l)]
    forms = ["+".join(f"({c})*{x}" for c, x in zip(v, names) if c % p) for v in normals]
    s = f'''needsPackage "HyperplaneArrangements";
R = ZZ/{p}[{",".join(names)}];
A = arrangement {{{",".join(forms)}}};
D = der(A, {{{",".join(map(str, mult))}}}); M = prune image D;
print("FREE " | toString isFreeModule M);
print("GENS " | toString flatten degrees source mingens image D);'''
    out = run_m2(s)
    return "FREE true" in out, [l for l in out.splitlines() if l.startswith("GENS")][0][5:]
braid = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]] + [[1 if k==i else (-1 if k==j else 0) for k in range(4)] for i in range(4) for j in range(i+1,4)]
B4 = [[1 if k==i else 0 for k in range(4)] for i in range(4)] + [[1 if k==i else (s if k==j else 0) for k in range(4)] for i in range(4) for j in range(i+1,4) for s in (1,-1)]
rng = random.Random(2)
cases = []
for p in [5, 7]:
    cases.append((p, braid, [1]*len(braid)))
    cases.append((p, B4, [1]*len(B4)))
    sub = rng.sample(B4, 12); cases.append((p, sub, [1]*12))
    A3 = [[1,0,0],[0,1,0],[0,0,1],[1,-1,0],[1,0,-1],[0,1,-1]]
    cases.append((p, A3, [2,2,2,2,2,2]))
    cases.append((p, A3, [1,2,3,1,2,1]))
    cases.append((p, [[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,2]], [3,3,3,1,1,1]))
    cases.append((p, braid, [2]*len(braid)))
bad = 0
for p, N, m in cases:
    Nm = [[c % p for c in v] for v in N]
    st, g = ExactModule(np.array(Nm), m, p).freeness()
    f2, g2 = m2(Nm, m, p)
    ok = (st == 'free') == f2
    bad += not ok
    print(p, len(N[0]), len(N), m[:6], '-> ours', st, g, '| M2', f2, g2, 'OK' if ok else 'MISMATCH', flush=True)
print('mismatches', bad)
