import sys, random, itertools, time
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice
from hunter.ffree import ExactModule, free_simple, free_rank3_ziegler
from hunter.m2 import run_m2
from math import comb

def m2_check(normals, p, kmax):
    l = len(normals[0]); names = [f"x{i}" for i in range(l)]
    forms = ["+".join(f"({c})*{x}" for c, x in zip(v, names) if c % p) for v in normals]
    s = f'''needsPackage "HyperplaneArrangements";
R = ZZ/{p}[{",".join(names)}];
A = arrangement {{{",".join(forms)}}};
D = der A; M = prune image D;
print("FREE " | toString isFreeModule M);
print("HF " | toString apply(0..{kmax}, k -> hilbertFunction(k, M)));'''
    out = run_m2(s)
    free = "FREE true" in out
    hf = [int(x) for x in [l for l in out.splitlines() if l.startswith('HF')][0][3:].strip(' (){}').split(',')]
    return free, hf

def pg_lines(p, l):
    out = []
    for v in itertools.product(range(p), repeat=l):
        if any(v) and next(c for c in v if c) == 1:
            out.append(list(v))
    return out

rng = random.Random(5)
mismatch = 0; tested = 0
for p in [5, 7, 11]:
    for l in [3, 4]:
        allv = pg_lines(p, l)
        cnt = 0; tries = 0
        while cnt < 6 and tries < 4000:
            tries += 1
            n = rng.randint(l + 2, min(len(allv), 14 if l == 3 else 11))
            S = rng.sample(allv, n)
            N = np.array(S, dtype=np.int64)
            if np.linalg.matrix_rank(N.astype(float)) < l: continue
            L = Lattice(N % p, p)
            ex = L.exponents()
            if ex is None and rng.random() < 0.8: continue   # prefer split chi
            cnt += 1
            st, g = free_simple(N, p, L)
            kmax = 6
            M = ExactModule(N, [1]*n, p)
            hf = [M.dim(k) for k in range(kmax + 1)]
            free_m2, hf_m2 = m2_check(S, p, kmax)
            z = free_rank3_ziegler(N, p, L) if l == 3 else '-'
            ok = (st == 'free') == free_m2 and hf == hf_m2 and (l == 4 or (z == 'free') == free_m2)
            tested += 1
            if not ok: mismatch += 1
            print(f'p={p} l={l} n={n} exps={ex} ours={st} {g} zieg={z} | M2 free={free_m2} | HF ok={hf == hf_m2}', 'OK' if ok else 'MISMATCH', flush=True)
print('tested', tested, 'mismatches', mismatch)
