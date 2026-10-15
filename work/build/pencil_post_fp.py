#!/usr/bin/env python3
"""Group F_p census hits by lattice isomorphism class; write one incidence file per class for the
exact realisation-space computation over Q (pencil_real2.m2)."""
import sys, glob, os
from collections import defaultdict
from pynauty import Graph, certificate

def parse(line):
    t = line.split(); sizes = tuple(int(x) for x in t[2:6]); p = int(t[7]); mu = int(t[9])
    i = t.index('A'); A = [int(x) for x in t[i+1:i+1+sizes[0]]]
    i = t.index('B'); B = [int(x) for x in t[i+1:i+1+sizes[1]]]
    i = t.index('C'); C = [int(x) for x in t[i+1:i+1+sizes[2]]]
    i = t.index('D'); D = [int(x) for x in t[i+1:i+1+sizes[3]]]
    return sizes, p, mu, A, B, C, D

def lattice_fp(sizes, p, mu, A, B, C, D):
    dirs = [(1,0),(0,1),(1,1),(mu,1)]
    L = [(0,0,1)]
    for d,(pp,qq) in enumerate(dirs):
        for v in [A,B,C,D][d]:
            L.append((qq % p, (-pp) % p, (-v) % p))
    n = len(L); pts = {}
    for i in range(n):
        for j in range(i+1, n):
            u, v = L[i], L[j]
            X = (u[1]*v[2]-u[2]*v[1]) % p; Y = (u[2]*v[0]-u[0]*v[2]) % p; Z = (u[0]*v[1]-u[1]*v[0]) % p
            for c in (X,Y,Z):
                if c:
                    s = pow(c, p-2, p); key = (X*s % p, Y*s % p, Z*s % p); break
            pts.setdefault(key, set()).update((i,j))
    return n, [frozenset(s) for s in pts.values()]

def canon(n, pts):
    m = len(pts); g = Graph(n+m, directed=False)
    for k, S in enumerate(pts): g.connect_vertex(n+k, sorted(S))
    mults = sorted(set(len(S) for S in pts))
    g.set_vertex_coloring([set(range(n))] + [set(n+k for k,S in enumerate(pts) if len(S)==mu) for mu in mults])
    return certificate(g)

classes = defaultdict(list); n_hits = 0
for f in sorted(glob.glob(sys.argv[1] if len(sys.argv) > 1 else 'fp_runs/*.out')):
    for line in open(f):
        if not line.startswith('FOUND'): continue
        n_hits += 1
        sizes, p, mu, A, B, C, D = parse(line)
        n, pts = lattice_fp(sizes, p, mu, A, B, C, D)
        assert sum(len(S)-1 for S in pts) == 80
        key = (sizes, canon(n, pts))
        classes[key].append((p, mu, A, B, C, D, pts))
print(f"hits {n_hits}, lattice classes {len(classes)}")
os.makedirs('fp_classes', exist_ok=True)
for k,(key, reals) in enumerate(sorted(classes.items(), key=lambda kv: kv[0][0])):
    sizes = key[0]; ps = sorted(set(r[0] for r in reals)); mus = sorted(set((r[0], r[1]) for r in reals))
    pts = reals[0][6]; n = 1 + sum(sizes)
    mult = sorted((len(S) for S in pts), reverse=True)
    name = f"fp_classes/class{k:04d}_{''.join(map(str,sizes))}.txt"
    with open(name, 'w') as out:
        out.write(f"sizes {' '.join(map(str,sizes))}\n")
        out.write(f"realisations {len(reals)} primes {ps} mus {mus}\n")
        for S in pts:
            if len(S) >= 3: out.write("point " + " ".join(map(str, sorted(S))) + "\n")
    print(f"{name}: sizes {sizes}, {len(reals)} realisations over p={ps}, (p,mu)={mus[:6]}{'...' if len(mus)>6 else ''}, mults {dict((m, mult.count(m)) for m in set(mult))}")
