#!/usr/bin/env python3
"""Additive three-pencil arrangements over F_p: lines y = a z (a in A), x = c z (c in B), x - y = b z (b in C), subsets of F_p,
optionally the common line z = 0; triple point iff c - a = b. The F_p-point has mdr <= p. We keep lattices with integer d1 > p,
p not dividing n, compute the realization space over F_p-bar (a linear space modulo the 3-dim symmetry group) and, when it is
positive-dimensional, test freeness of a random non-rational point over GF(p^k) in M2.  Usage: additive.py p"""
import sys, itertools, math, json, random
p = int(sys.argv[1]); random.seed(1)
from collections import defaultdict
def rank_mod(rows, m, mod):
    rows = [r[:] for r in rows]; rk = 0; col = 0
    for col in range(m):
        piv = None
        for i in range(rk, len(rows)):
            if rows[i][col] % mod: piv = i; break
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        inv = pow(rows[rk][col], -1, mod)
        rows[rk] = [(v*inv) % mod for v in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][col] % mod:
                f = rows[i][col]; rows[i] = [(rows[i][j] - f*rows[rk][j]) % mod for j in range(m)]
        rk += 1
    return rk
cands = {}
slots = 3*p
examined = 0
for zline in (0,1):
    for r in range(0, slots+1):
        n = slots - r + zline
        if n < 2*p + 3: break
        if n % p == 0: continue
        for O in itertools.combinations(range(slots), r):
            Os = set(O)
            A = [i for i in range(p) if i not in Os]; B = [i for i in range(p) if p+i not in Os]; C = [i for i in range(p) if 2*p+i not in Os]
            if not A or not B or not C: continue
            examined += 1
            Cs = set(C); t = 0; T = []
            for i,a in enumerate(A):
                for j,c in enumerate(B):
                    b = (c - a) % p
                    if b in Cs: t += 1; T.append((i,j,C.index(b)))
            m1, m2, m3 = len(A)+zline, len(B)+zline, len(C)+zline
            b2 = n*(n-1)//2 - sum((m-1)*(m-2)//2 for m in (m1,m2,m3)) - t
            disc = (n-1)**2 - 4*(b2-n+1)
            if disc < 0: continue
            s = math.isqrt(disc)
            if s*s != disc or (n-1-s) % 2: continue
            d1 = (n-1-s)//2
            if d1 <= p: continue
            # realization space over F_p-bar: unknowns a_i, c_j, b_k; equations c_j - a_i - b_k = 0 for triples;
            # symmetries: a -> l a + u, c -> l c + v, b -> l b + v - u  (3-dim)
            nA,nB,nC = len(A),len(B),len(C); mvars = nA+nB+nC
            rows = []
            for (i,j,k) in T:
                row = [0]*mvars; row[nA+j] = 1; row[i] = -1; row[nA+nB+k] = -1; rows.append(row)
            rk = rank_mod(rows, mvars, p)
            dim = mvars - rk - 3
            key = (n, d1, zline, tuple(sorted((len(A),len(B),len(C)))), t, dim)
            if key not in cands: cands[key] = dict(A=A,B=B,C=C,z=zline,n=n,d1=d1,d2=n-1-d1,t=t,dim=dim,T=T)
print(f'p={p}: {examined} examined; {len(cands)} candidate lattices (d1 > p): dims {sorted(set(k[5] for k in cands))}', file=sys.stderr)
for k,c in sorted(cands.items()):
    print(f"  n={c['n']} class=({c['d1']},{c['d2']}) z={c['z']} sizes={k[3]} t={c['t']} dim(real)={c['dim']}", file=sys.stderr)
json.dump(list(cands.values()), open(f'additive_p{p}.json','w'))
