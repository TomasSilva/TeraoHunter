import sys, time, random
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice, freeness, hilbert_DH, local_freeness_l4, DEFAULT_P
from hunter.m2 import m2_free

p = DEFAULT_P
def test(name, normals, kmax=None):
    N = np.array(normals, dtype=np.int64)
    t = time.time()
    L = Lattice(N % p, p)
    r = freeness(N % p, p, L)
    t1 = time.time() - t
    loc = local_freeness_l4(N % p, p, L) if N.shape[1] == 4 else None
    kk = kmax or (max(r['exps']) if r['exps'] else 6)
    hf = hilbert_DH(N % p, p, kk)
    m = m2_free([list(map(int, v)) for v in normals], hilbert_upto=kk)
    # M2 HF is of D(A); ours is D_H; D = S E + D_H  => HF_D(k) = C(k+l-2, l-1) + HF_DH(k)
    from math import comb
    l = N.shape[1]
    hfD = [ (comb(k - 1 + l - 1, l - 1) if k >= 1 else 0) + h for k, h in enumerate(hf)]
    print(f"{name}: chi={L.char_poly()} exps={r['exps']} status={r['status']} ({t1:.2f}s) locbad={len(loc) if loc is not None else '-'}")
    print(f"   ours HF_D={hfD}")
    print(f"   M2   HF_D={m.get('hf')} free={m['free']} gens={m['gens']}")
    assert hfD == m.get('hf'), "HILBERT MISMATCH"
    assert (r['status'] == 'free') == m['free'], "FREENESS MISMATCH"

# braid A4 in K^4 (x0 = 0 convention)
braid = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
for i in range(4):
    for j in range(i+1,4):
        v = [0]*4; v[i]=1; v[j]=-1; braid.append(v)
test("A4 braid", braid)
# B4 type
B4 = [[1 if k==i else 0 for k in range(4)] for i in range(4)]
for i in range(4):
    for j in range(i+1,4):
        for s in (1,-1):
            v=[0]*4; v[i]=1; v[j]=s; B4.append(v)
test("B4", B4)
# generic 6 planes
random.seed(1)
gen = [[random.randint(-9,9) for _ in range(4)] for _ in range(6)]
test("generic6", gen)
# random arrangements with small coefficients
for trial in range(6):
    S = set()
    while len(S) < random.randint(7, 11):
        v = tuple(random.randint(-1, 1) for _ in range(4))
        if any(v):
            # normalise sign
            f = next(c for c in v if c)
            v = tuple(c*f for c in v)
            S.add(v)
    test(f"rand{trial}", [list(v) for v in S])
# rank 3 tests
test("A3 lines", [[1,0,0],[0,1,0],[0,0,1],[1,-1,0],[1,0,-1],[0,1,-1]])
test("generic 5 lines", [[1,0,0],[0,1,0],[0,0,1],[1,1,1],[1,2,3]])
