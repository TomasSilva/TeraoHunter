#!/usr/bin/env python3
"""Group the FOUND configurations of pencils.c by lattice isomorphism class (nauty), report classes
with several inequivalent realisations or realisations at rational jump values."""
import sys, glob, re
from fractions import Fraction as Fr
from collections import defaultdict
from pynauty import Graph, certificate

def parse(line):
    t = line.split()
    sizes = tuple(int(x) for x in t[2:6]); p, q = t[7].split('/'); lam = Fr(int(p), int(q))
    i = t.index('A'); A = [Fr(x) for x in t[i+1:i+1+sizes[0]]]
    i = t.index('B'); B = [Fr(x) for x in t[i+1:i+1+sizes[1]]]
    i = t.index('C'); C = [Fr(x) for x in t[i+1:i+1+sizes[2]]]
    i = t.index('D'); D = [Fr(x) for x in t[i+1:i+1+sizes[3]]]
    return sizes, lam, A, B, C, D

def lines_of(sizes, lam, A, B, C, D):
    # line: (q, -p, -v) for direction (p,q); H = (0,0,1)
    dirs = [(Fr(1),Fr(0)), (Fr(0),Fr(1)), (Fr(1),Fr(1)), (lam, Fr(1))]
    L = [(Fr(0),Fr(0),Fr(1))]
    for d,(p,q) in enumerate(dirs):
        for v in [A,B,C,D][d]:
            L.append((q, -p, -v))
    return L

def cross(u, v):
    return (u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0])
def norm(v):
    for c in v:
        if c != 0:
            return tuple(x/c for x in v)
    return None

def lattice(L):
    n = len(L); pts = {}
    for i in range(n):
        for j in range(i+1, n):
            P = norm(cross(L[i], L[j]))
            pts.setdefault(P, set()).update((i, j))
    return [frozenset(s) for s in pts.values()]

def canon(n, pts):
    # bipartite graph: lines 0..n-1, points n.., coloured by multiplicity
    m = len(pts)
    g = Graph(n + m, directed=False)
    adj = {}
    for k, S in enumerate(pts):
        adj[n+k] = list(S)
    for k in range(n+m):
        g.connect_vertex(k, adj.get(k, []))
    mults = sorted(set(len(S) for S in pts))
    parts = [set(range(n))] + [set(n+k for k,S in enumerate(pts) if len(S)==mu) for mu in mults]
    g.set_vertex_coloring(parts)
    return certificate(g)

if __name__ == '__main__':
    classes = defaultdict(list)
    files = sorted(glob.glob('pencil_runs/*.out'))
    nfound = 0
    for f in files:
        for line in open(f):
            if not line.startswith('FOUND'): continue
            nfound += 1
            sizes, lam, A, B, C, D = parse(line)
            L = lines_of(sizes, lam, A, B, C, D)
            pts = lattice(L)
            b2 = sum(len(S)-1 for S in pts)
            assert b2 == 80, (b2, line)
            key = (sizes, canon(len(L), pts))
            classes[key].append((lam, A, B, C, D, f))
    print(f"configurations: {nfound}, lattice classes: {len(classes)}")
    S_rat = {Fr(-1), Fr(2), Fr(1,2)}
    for key, reals in sorted(classes.items(), key=lambda kv: (kv[0][0], len(kv[1]))):
        sizes = key[0]
        lams = sorted(set(r[0] for r in reals))
        # distinct realisations up to the affine normalisation used by pencils.c: compare the offset tuples
        distinct = set((r[0], tuple(r[1]), tuple(r[2]), tuple(r[3]), tuple(r[4])) for r in reals)
        flag = []
        if len(lams) > 1: flag.append("SEVERAL lambda")
        if len(distinct) > 1: flag.append(f"{len(distinct)} distinct realisations")
        if any(l in S_rat for l in lams): flag.append("lambda in rational jump set")
        print(f"sizes {sizes}: {len(reals)} hits, lambdas {[str(l) for l in lams]}, example A={[str(x) for x in reals[0][1]]} B={[str(x) for x in reals[0][2]]} C={[str(x) for x in reals[0][3]]} D={[str(x) for x in reals[0][4]]}  {' | '.join(flag)}")
