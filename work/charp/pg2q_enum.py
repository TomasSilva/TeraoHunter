#!/usr/bin/env python3
"""Sub-arrangements of PG(2,q) (q = p^k) with n lines, p not dividing n, whose chibar has integer roots d1 <= d2 with d1 > q.
Such an arrangement is F_q-rational, so mdr <= q < d1: non-free. We list the lattices (up to isomorphism invariants).
Usage: pg2q_enum.py q [nmax]   (exhaustive over all n-subsets; q = 4 is feasible, q = 5 only for n close to 31)."""
import sys, itertools
from collections import defaultdict, Counter
q = int(sys.argv[1])
# build GF(q)
if q in (2,3,5,7,11,13):
    els = list(range(q)); add = lambda a,b: (a+b)%q; mul = lambda a,b: (a*b)%q; neg = lambda a: (-a)%q
    inv = lambda a: pow(a, q-2, q)
elif q == 4:
    # GF(4) = {0,1,w,w+1} as ints 0..3 with bits; w^2 = w+1
    els = [0,1,2,3]; add = lambda a,b: a ^ b
    def mul(a,b):
        r = 0
        for i in range(2):
            if (b>>i)&1: r ^= a<<i
        # reduce by x^2 + x + 1
        if r & 4: r ^= 0b111
        return r
    neg = lambda a: a
    inv = lambda a: [None,1,3,2][a]
elif q == 8:
    els = list(range(8)); add = lambda a,b: a ^ b
    def mul(a,b):
        r = 0
        for i in range(3):
            if (b>>i)&1: r ^= a<<i
        for s in (5,4,3):
            if r & (1<<s): r ^= 0b1011 << (s-3)
        return r
    neg = lambda a: a
    def inv(a):
        for b in range(1,8):
            if mul(a,b) == 1: return b
elif q == 9:
    # GF(9) = F_3[i]/(i^2+1): a + b i -> 3a + b
    els = list(range(9))
    def add(a,b): return 3*((a//3 + b//3)%3) + ((a%3 + b%3)%3)
    def mul(a,b):
        a0,a1 = a//3, a%3; b0,b1 = b//3, b%3
        return 3*((a0*b0 - a1*b1)%3) + ((a0*b1 + a1*b0)%3)
    def neg(a): return 3*((-(a//3))%3) + ((-(a%3))%3)
    def inv(a):
        for b in range(1,9):
            if mul(a,b) == 1: return b
else:
    raise SystemExit('unsupported q')
# points and lines of PG(2,q) as normalized triples
def normalize(v):
    for i in range(3):
        if v[i] != 0:
            iv = inv(v[i]); return tuple(mul(x, iv) for x in v)
    return None
pts = sorted(set(normalize(v) for v in itertools.product(els, repeat=3) if any(v)))
lines = pts[:]  # same coordinates (duality)
def incident(P, L):
    s = 0
    for a,b in zip(P,L): s = add(s, mul(a,b))
    return s == 0
NP = len(pts)
inc = [[j for j,L in enumerate(lines) if incident(P,L)] for P in pts]      # lines through each point
linepts = [[i for i,P in enumerate(pts) if incident(P,L)] for L in lines]  # points on each line
print(f'PG(2,{q}): {NP} points/lines, {len(inc[0])} lines per point', file=sys.stderr)
p = [r for r in range(2,q+1) if q % r == 0 and all(r % s for s in range(2,r))][0]
nmax = int(sys.argv[2]) if len(sys.argv) > 2 else NP
results = defaultdict(list)
import math
for n in range(2*q+3, nmax+1):
    if n % p == 0: continue
    cnt = 0; found = 0
    for S in itertools.combinations(range(NP), n):
        cnt += 1
        Sset = set(S)
        mult = [sum(1 for j in inc[i] if j in Sset) for i in range(NP)]
        b2 = sum(m-1 for m in mult if m >= 2)
        disc = (n-1)**2 - 4*(b2 - n + 1)
        if disc < 0: continue
        s = math.isqrt(disc)
        if s*s != disc or (n-1-s) % 2: continue
        d1 = (n-1-s)//2
        if d1 <= q: continue
        prof = tuple(sorted(Counter(m for m in mult if m >= 2).items()))
        # finer invariant: for each line, the sorted multiplicities of its points
        linv = tuple(sorted(tuple(sorted(mult[i] for i in linepts[j])) for j in S))
        key = (n, d1, prof, linv)
        if len(results[key]) < 1:
            results[key].append(S)
        found += 1
    print(f'n={n}: {cnt} subsets, {found} with integer roots and d1 > q, {sum(1 for k in results if k[0]==n)} invariant classes', file=sys.stderr)
import json
json.dump([dict(n=k[0], d1=k[1], profile=k[2], S=v[0]) for k,v in results.items()], open(f'pg2_{q}_cands.json','w'))
for k in sorted(results):
    print(k[0], 'class', (k[1], k[0]-1-k[1]), 'profile', k[2])
