#!/usr/bin/env python3
"""Sub-case (I), vector (5,5,5,1), via additive structures.
Net: A-lines y = a z (a in A), B-lines x = c z (c in B), C-lines x - y = b z (b in C); triple point (a,c,b) iff c - a = b.
A, B subsets of [0,D] containing 0, |A|=|B|=|C|=5, t = #{(a,c): c-a in C} in {16,17,18}.  For each net (up to isotopy of the
abstract triple system) we compute the dimension of its realization space V (solutions of c_j - a_i - b_k = 0), then all
possible extra lines L4 through P4 = [lam:1:0]: L4 passes through a set of net points (AB-cells, AC- or BC-doubles), at most
one per line, with the window s_H in [4,8] for every line and b2 = 80 (t = 10 + s_L4).  For rigid nets lam, d are determined;
Lemma F2 decides freeness: free iff lam not in {-1, 2, 1/2} u {roots of t^2-t+1}.  Non-rigid nets are reported."""
import itertools, sys, math
from fractions import Fraction as Fr
from collections import Counter, defaultdict
D = int(sys.argv[1]) if len(sys.argv) > 1 else 12
nets = {}
def canon(T):
    best = None
    for pr in itertools.permutations(range(5)):
        for pc in itertools.permutations(range(5)):
            relabel = {}; key = []
            for i in range(5):
                row = []
                for j in range(5):
                    k = T.get((pr[i], pc[j]))
                    if k is None: row.append(-1)
                    else:
                        if k not in relabel: relabel[k] = len(relabel)
                        row.append(relabel[k])
                key.append(tuple(row))
            key = tuple(key)
            if best is None or key < best: best = key
    return best
def rank_q(rows, ncols):
    M = [r[:] for r in rows]; rk = 0
    for col in range(ncols):
        piv = next((i for i in range(rk, len(M)) if M[i][col] != 0), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]
        inv = Fr(1) / M[rk][col]; M[rk] = [x*inv for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][col] != 0:
                f = M[i][col]; M[i] = [x - f*y for x, y in zip(M[i], M[rk])]
        rk += 1
    return rk
cnt = 0
for Ar in itertools.combinations(range(1, D+1), 4):
    A = (0,) + Ar
    for Br in itertools.combinations(range(1, D+1), 4):
        B = (0,) + Br
        r = Counter(c - a for a in A for c in B)
        vals = sorted(r.items(), key=lambda x: -x[1])
        if sum(v for _, v in vals[:5]) < 16: continue
        # choose C: 5 values with sum in {16,17,18}
        cnts = [v for _, v in vals]; keys = [k for k, _ in vals]
        def rec(start, k, lo, hi, chosen):
            if k == 0:
                s = sum(cnts[i] for i in chosen)
                if lo <= s <= hi: yield list(chosen)
                return
            for idx in range(start, len(vals)):
                if sum(cnts[idx:idx+k]) + sum(cnts[i] for i in chosen) < lo: break
                chosen.append(idx); yield from rec(idx+1, k-1, lo, hi, chosen); chosen.pop()
        for ch in rec(0, 5, 16, 18, []):
            C = tuple(keys[i] for i in ch)
            T = {}
            for i, a in enumerate(A):
                for j, c in enumerate(B):
                    if (c - a) in C: T[(i,j)] = C.index(c - a)
            rows = Counter(i for i,j in T); cols = Counter(j for i,j in T); syms = Counter(T.values())
            if min(rows[i] for i in range(5)) < 3 or min(cols[j] for j in range(5)) < 3 or min(syms[k] for k in range(5)) < 3: continue
            key = canon(T); cnt += 1
            if key in nets: continue
            # realization space dimension of the net
            eqs = []
            for (i,j), k in T.items():
                v = [Fr(0)]*15; v[5+j] = Fr(1); v[i] = Fr(-1); v[10+k] = Fr(-1); eqs.append(v)
            dimV = 15 - rank_q(eqs, 15)
            nets[key] = dict(A=A, B=B, C=C, T=T, t=len(T), dimV=dimV)
print(f'D={D}: {cnt} (A,B,C) structures, {len(nets)} nets up to isotopy; t distribution {sorted(Counter(n["t"] for n in nets.values()).items())}; dimV distribution {sorted(Counter(n["dimV"] for n in nets.values()).items())}', file=sys.stderr)
import pickle; pickle.dump(nets, open(f'nets551_D{D}.pkl','wb'))
for key, n in nets.items():
    print(f"t={n['t']} dimV={n['dimV']} A={n['A']} B={n['B']} C={n['C']}", file=sys.stderr)
