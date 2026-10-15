#!/usr/bin/env python3
"""Group the mu-flexible four-pencil hits (FLEX lines of ag8_flex) by lattice isomorphism class
(nauty on the line-point incidence graph); print one representative per class with counts."""
import sys, glob
from collections import OrderedDict
from pynauty import Graph, certificate

def mul8(a, b):
    r = 0
    for i in range(3):
        if b >> i & 1: r ^= a << i
    for i in (4, 3):
        if r >> i & 1: r ^= 0xB << (i - 3)
    return r
INV = {a: next(b for b in range(1, 8) if mul8(a, b) == 1) for a in range(1, 8)}

def lines(S):
    L = [(0, 0, 1)]
    for j, m in enumerate(S):
        for v in range(8):
            if m >> v & 1:
                L.append([(0, 1, v), (1, 0, v), (1, 1, v), (1, 2, v)][j])   # t = 2
    return L

def cross(u, v):
    X = mul8(u[1], v[2]) ^ mul8(u[2], v[1]); Y = mul8(u[2], v[0]) ^ mul8(u[0], v[2]); Z = mul8(u[0], v[1]) ^ mul8(u[1], v[0])
    for c in (X, Y, Z):
        if c:
            i = INV[c]; return (mul8(X, i), mul8(Y, i), mul8(Z, i))

def lattice(L):
    pts = {}
    for i in range(len(L)):
        for j in range(i + 1, len(L)):
            pts.setdefault(cross(L[i], L[j]), set()).update((i, j))
    return [frozenset(s) for s in pts.values()]

def canon(n, pts):
    m = len(pts); g = Graph(n + m, directed=False)
    for k, S_ in enumerate(pts): g.connect_vertex(n + k, sorted(S_))
    mults = sorted(set(len(S_) for S_ in pts))
    g.set_vertex_coloring([set(range(n))] + [set(n + k for k, S_ in enumerate(pts) if len(S_) == mu) for mu in mults])
    return certificate(g)

classes = OrderedDict()
for f in sorted(glob.glob(sys.argv[1] if len(sys.argv) > 1 else 'flex_*.out')):
    for line in open(f):
        if not line.startswith('FLEX'): continue
        t = line.split(); S = [int(x, 16) for x in t[1:5]]
        L = lines(S); pts = lattice(L)
        key = canon(len(L), pts)
        if key not in classes: classes[key] = [S, 0, line.strip()]
        classes[key][1] += 1
for k, (S, cnt, line) in enumerate(classes.values()):
    print(k, cnt, line)
print('classes', len(classes), file=sys.stderr)
