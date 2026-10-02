import sys, random
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import multi_freeness, module_generators, DEFAULT_P, freeness
from hunter.m2 import run_m2
p = DEFAULT_P

def m2_multi(normals, mult):
    l = len(normals[0])
    names = [f"x{i}" for i in range(l)]
    forms = []
    for v in normals:
        forms.append("+".join(f"({c})*{x}" for c, x in zip(v, names) if c != 0))
    s = f'''needsPackage "HyperplaneArrangements";
R = QQ[{",".join(names)}];
A = arrangement {{{",".join(forms)}}};
D = der(A, {{{",".join(map(str,mult))}}});
M = prune image D;
print("FREE " | toString isFreeModule M);
print("GENS " | toString flatten degrees source mingens image D);
'''
    out = run_m2(s)
    free = "FREE true" in out
    gens = [l for l in out.splitlines() if l.startswith("GENS")][0][5:]
    return free, gens

tests = [
  ([[1,0,0],[0,1,0],[0,0,1]], [2,3,4]),
  ([[1,0,0],[0,1,0],[0,0,1],[1,1,1]], [1,1,1,1]),
  ([[1,0,0],[0,1,0],[0,0,1],[1,1,1]], [1,1,1,3]),
  ([[1,0,0],[0,1,0],[0,0,1],[1,1,1]], [2,2,2,2]),
  ([[1,0,0],[0,1,0],[0,0,1],[1,-1,0],[1,0,-1],[0,1,-1]], [2,2,2,2,2,2]),
  ([[1,0,0],[0,1,0],[0,0,1],[1,-1,0],[1,0,-1],[0,1,-1]], [1,2,3,1,2,1]),
  ([[1,0],[0,1],[1,1],[1,2]], [2,2,2,2]),
  ([[1,0],[0,1],[1,1],[1,-1]], [2,2,2,2]),
  ([[1,0],[0,1],[1,1],[1,3]], [3,1,2,2]),
]
for normals, mult in tests:
    st, gens = multi_freeness(np.array(normals) % p, mult, p)
    f2, g2 = m2_multi(normals, mult)
    print(normals, mult, '->', st, gens, '| M2:', f2, g2)
    assert (st == 'free') == f2
