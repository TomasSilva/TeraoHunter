#!/usr/bin/env python3
"""Group the FOUND5 hits of pencils5x.c by lattice isomorphism class; write census5_classes.jsonl
(blocks = multiple points with >= 3 lines, lines 0 = H, 1..16 the pencils) for exact_tc_lattice.py."""
import sys, glob, json
from collections import defaultdict
from pynauty import Graph, certificate

def parse(line):
    t = line.split()
    S = [int(x) for x in t[2:7]]; p = int(t[8]); l1 = int(t[10]); l2 = int(t[12])
    offs = []; i = 13
    for name, s in zip("ABCED", S):
        assert t[i] == name, (t[i], name, line); offs.append([int(x) for x in t[i+1:i+1+s]]); i += 1 + s
    return S, p, l1, l2, offs

def arrangement(S, p, l1, l2, offs):
    dirs = [(1, 0), (0, 1), (1, 1), (l1, 1), (l2, 1)]
    L = [(0, 0, 1)]
    for (dp, dq), o in zip(dirs, offs):
        for v in o: L.append((dq % p, (-dp) % p, (-v) % p))
    return L

def lattice(L, p):
    n = len(L); pts = {}
    for i in range(n):
        for j in range(i + 1, n):
            u, v = L[i], L[j]
            X = (u[1]*v[2]-u[2]*v[1]) % p; Y = (u[2]*v[0]-u[0]*v[2]) % p; Z = (u[0]*v[1]-u[1]*v[0]) % p
            for c in (X, Y, Z):
                if c:
                    s = pow(c, p - 2, p); key = (X*s % p, Y*s % p, Z*s % p); break
            pts.setdefault(key, set()).update((i, j))
    return [frozenset(s) for s in pts.values()]

def canon(n, pts):
    m = len(pts); g = Graph(n + m, directed=False)
    for k, S_ in enumerate(pts): g.connect_vertex(n + k, sorted(S_))
    mults = sorted(set(len(S_) for S_ in pts))
    g.set_vertex_coloring([set(range(n))] + [set(n + k for k, S_ in enumerate(pts) if len(S_) == mu) for mu in mults])
    return certificate(g)

if __name__ == '__main__':
    files = sorted(glob.glob(sys.argv[1] if len(sys.argv) > 1 else 'census5_runs/*.out'))
    classes = {}; info = defaultdict(lambda: {"hits": 0, "primes": set(), "lams": set()})
    nhits = 0
    for f in files:
        for line in open(f):
            if not line.startswith('FOUND5'): continue
            nhits += 1
            S, p, l1, l2, offs = parse(line)
            L = arrangement(S, p, l1, l2, offs); n = len(L)
            pts = lattice(L, p)
            b2 = sum(len(P) - 1 for P in pts)
            sH = [sum(1 for P in pts if a in P) for a in range(n)]
            assert b2 == 80 and min(sH) >= 5 and max(sH) <= 8 and sH[0] == 5, (b2, sH, line)
            key = canon(n, pts)
            if key not in classes: classes[key] = (S, pts)
            info[key]["hits"] += 1; info[key]["primes"].add(p); info[key]["lams"].add((p, l1, l2))
    print(f"hits {nhits}, lattice classes {len(classes)}")
    with open('census5_classes.jsonl', 'w') as out:
        for k, (key, (S, pts)) in enumerate(classes.items()):
            blocks = [sorted(P) for P in pts if len(P) >= 3]
            mult = sorted((len(P) for P in pts), reverse=True)
            prof = {m: mult.count(m) for m in sorted(set(mult))}
            sH = sorted(sum(1 for P in pts if a in P) for a in range(17))
            rec = {"blocks": blocks, "sizes": S, "primes": sorted(info[key]["primes"]), "hits": info[key]["hits"],
                   "profile": prof, "points_per_line": sH}
            out.write(json.dumps(rec) + "\n")
            print(f"class {k}: H-vector {sorted(S, reverse=True)}, primes {sorted(info[key]['primes'])}, hits {info[key]['hits']}, points {prof}, s_H {dict((x, sH.count(x)) for x in sorted(set(sH)))}")
