import sys, random, itertools
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.fq import GF, LineArr
from hunter.core import Lattice
from hunter.ffree import free_rank3_ziegler
from hunter.m2 import run_m2
rng = random.Random(3)
# 1. prime fields vs prime engine
bad = 0
for q in [7, 11]:
    F = GF(q); pts = F.points(3)
    for _ in range(200):
        S = rng.sample(pts, rng.randint(8, 20))
        try:
            A = LineArr(F, S)
        except ValueError:
            continue
        if A.exponents() is None: continue
        a = A.free()
        b = free_rank3_ziegler(np.array(S), q) == 'free'
        if a != b: bad += 1
print('prime-field mismatches:', bad)
# 2. GF(8), GF(9) vs Macaulay2 with the same modulus
def m2_free_fq(F, lines):
    p, k = F.p, F.k
    fpoly = '+'.join(f'({c})*s^{i}' for i, c in enumerate(F.f) if c)
    def elt(a):
        digits = [(a // p ** i) % p for i in range(k)]
        terms = [f'({d})*s^{i}' for i, d in enumerate(digits) if d]
        return '(' + ('+'.join(terms) if terms else '0') + ')'
    forms = ['+'.join(f'{elt(c)}*{x}' for c, x in zip(L, 'xyz') if c) for L in lines]
    s = f'''needsPackage "HyperplaneArrangements";
kk = GF(ZZ/{p}[s]/ideal({fpoly}), Variable => s);
R = kk[x,y,z];
A = arrangement {{{",".join(forms)}}};
print("FREE " | toString isFreeModule prune image der A);'''
    return 'FREE true' in run_m2(s)
for q in [8, 9]:
    F = GF(q); pts = F.points(3); cnt = 0; mism = 0
    while cnt < 12:
        S = rng.sample(pts, rng.randint(8, 18))
        try:
            A = LineArr(F, S)
        except ValueError:
            continue
        if A.exponents() is None and rng.random() < 0.85: continue
        cnt += 1
        ours = A.free()
        theirs = m2_free_fq(F, S)
        mism += ours != theirs
        print(q, len(S), A.exponents(), 'ours', ours, 'M2', theirs, '' if ours == theirs else 'MISMATCH', flush=True)
    print('GF(%d) mismatches' % q, mism)
