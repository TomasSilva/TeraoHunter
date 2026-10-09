#!/usr/bin/env python3
"""General 4-pencil enumeration for sub-case (I): a 4-point line H = {z=0} with pencils of sizes (nA,nB,nC,nD) through
P1=[1:0:0] (y = a z), P2=[0:1:0] (x = c z), P3=[1:1:0] (x - y = b z), P4=[lam:1:0] (x - lam y = d z).
Net = pencils A,B,C with additive structure (triple iff c - a = b).  D-lines pass through net points:
AB cell (a,c) at (c,a): c - lam a = d ; AC double (a,b) at (a+b, a): a + b - lam a = d ; BC double (c,b) at (c, c-b): c - lam(c-b) = d.
For rigid nets (dim V = 3) the normalized solution is rational; D-lines are determined by two of their points.
Usage: net4gen.py nA nB nC nD Dmax"""
import itertools, sys
from fractions import Fraction as Fr
from collections import Counter
nA, nB, nC, nD, Dmax = map(int, sys.argv[1:6])
pairs = nA*nB + nA*nC + nB*nC
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
smin = 1 + -(-(nA+nB+nC) // 3)   # each D-line: P4 plus >= ceil(net lines / 3) points
lattices = {}; nets_seen = 0; nonrigid = 0
stats = Counter()
for Ar in itertools.combinations(range(1, Dmax+1), nA-1):
    A = (0,) + Ar
    for Br in itertools.combinations(range(1, Dmax+1), nB-1):
        B = (0,) + Br
        r = Counter(c - a for a in A for c in B)
        keys = sorted(r)
        for C in itertools.combinations(keys, nC):
            Cs = set(C)
            T = {(i,j): C.index(c-a) for i,a in enumerate(A) for j,c in enumerate(B) if (c-a) in Cs}
            t = len(T)
            # b2 = 80 : t = nA+nB+nC+nD + pairs - 80 + sum(s_i - 1), s_i in [smin, 8]
            base = nA+nB+nC+nD + pairs - 80
            if not (base + nD*(smin-1) <= t <= base + nD*7): continue
            # per-line degree necessary conditions (delta >= 0)
            ta = Counter(i for i,j in T); tc = Counter(j for i,j in T); tb = Counter(T.values())
            if any(ta[i] < nB + nC - 7 for i in range(nA)): continue
            if any(tc[j] < nA + nC - 7 for j in range(nB)): continue
            if any(tb[k] < nA + nB - 7 for k in range(nC)): continue
            nets_seen += 1
            # realization space of the net
            NV = nA + nB + nC
            eqs = []
            for (i,j), k in T.items():
                v = [Fr(0)]*NV; v[nA+j] = Fr(1); v[i] = Fr(-1); v[nA+nB+k] = Fr(-1); eqs.append(v)
            dimV = NV - rank_q(eqs, NV)
            if dimV > 3:
                nonrigid += 1; stats['nonrigid nets'] += 1
                print(f'NONRIGID net dimV={dimV} t={t} A={A} B={B} C={C}', file=sys.stderr); continue
            # rigid: the given (A,B,C) is the realization (unique up to affine maps).  Net points with coordinates:
            pts = {}   # name -> (x, y, set of pencil-line ids)
            for i,a in enumerate(A):
                for j,c in enumerate(B):
                    name = ('AB', i, j); lines_ = {('A',i), ('B',j)}
                    if (i,j) in T: lines_.add(('C', T[(i,j)]))
                    pts[name] = (Fr(c), Fr(a), lines_)
            for i,a in enumerate(A):
                for k,b in enumerate(C):
                    if (a + b) in B: continue   # triple, already an AB cell
                    pts[('AC', i, k)] = (Fr(a+b), Fr(a), {('A',i), ('C',k)})
            for j,c in enumerate(B):
                for k,b in enumerate(C):
                    if (c - b) in A: continue
                    pts[('BC', j, k)] = (Fr(c), Fr(c-b), {('B',j), ('C',k)})
            names = list(pts)
            # candidate D-lines: through P4=[lam:1:0] i.e. x - lam y = d ; determined by two net points with different y
            # (if same y, lam infinite -> P4 = P1, excluded)
            cand = {}
            for n1, n2 in itertools.combinations(names, 2):
                x1,y1,_ = pts[n1]; x2,y2,_ = pts[n2]
                if y1 == y2: continue
                lam = (x1 - x2) / (y1 - y2); d = x1 - lam*y1
                if lam in (0, 1): continue
                cand.setdefault(lam, set()).add(d)
            for lam, ds in cand.items():
                # all D-lines for this lam: for each d, the set of net points on x - lam y = d
                dl = {}
                for d in ds:
                    on = [n for n in names if pts[n][0] - lam*pts[n][1] == d]
                    # check: at most one line of each pencil through the points on this D-line (automatic since distinct points on a line...)
                    covered = set().union(*[pts[n][2] for n in on])
                    if len(covered) != sum(len(pts[n][2]) for n in on): continue  # two points sharing a pencil line: impossible geometrically unless same point
                    s = 1 + len(on) + (nA + nB + nC - len(covered))   # P4 + net points + new double points with uncovered lines
                    if s <= 8: dl[d] = (on, covered, s)
                if len(dl) < nD: continue
                for dsel in itertools.combinations(sorted(dl), nD):
                    S = sum(dl[d][2] - 1 for d in dsel)
                    if t != base + S: continue
                    # per-line counts
                    ok = True
                    for i in range(nA):
                        delta = sum(1 for d in dsel if ('A',i) not in dl[d][1])
                        if 1 + nB + (nC - ta[i]) + delta > 8: ok = False; break
                    if not ok: continue
                    for j in range(nB):
                        delta = sum(1 for d in dsel if ('B',j) not in dl[d][1])
                        if 1 + nA + (nC - tc[j]) + delta > 8: ok = False; break
                    if not ok: continue
                    for k in range(nC):
                        delta = sum(1 for d in dsel if ('C',k) not in dl[d][1])
                        if 1 + (nA - tb[k]) + (nB - tb[k]) + tb[k] + delta > 8: ok = False; break
                    if not ok: continue
                    key = (A, B, C, lam, dsel)
                    lattices[key] = dict(t=t, lam=lam, s=[dl[d][2] for d in dsel])
                    stats[('lam', lam)] += 1
print(f'sizes ({nA},{nB},{nC},{nD}), Dmax={Dmax}: nets passing counts {nets_seen}, nonrigid {nonrigid}, lattices (rigid, with D-lines) {len(lattices)}', file=sys.stderr)
print('cross-ratios lam found:', sorted(Counter(v['lam'] for v in lattices.values()).items(), key=lambda x: -x[1])[:20], file=sys.stderr)
import pickle; pickle.dump(lattices, open(f'lat4_{nA}{nB}{nC}{nD}_D{Dmax}.pkl','wb'))
