#!/usr/bin/env python3
"""Lattice isomorphism classes of the 17-line window deletions of A(G26) and A(6,1,3); write jsonl for exact_tc_lattice.py"""
import json
from collections import defaultdict
from pynauty import Graph, certificate
def load_lattice(fn):
    return [tuple(int(x) for x in l.split()[1:]) for l in open(fn) if l.startswith('point')]
def canon(n, blocks):
    m = len(blocks); g = Graph(n+m, directed=False)
    for k,S in enumerate(blocks): g.connect_vertex(n+k, sorted(S))
    mults = sorted(set(len(S) for S in blocks))
    g.set_vertex_coloring([set(range(n))] + [set(n+k for k,S in enumerate(blocks) if len(S)==mu) for mu in mults])
    return certificate(g)
classes = {}; count = defaultdict(int); src = {}
for name, nl in (('g26', 21), ('fermat', 21)):
    pts = load_lattice(f'{name}_lattice.txt')
    for line in open(f'{name}_dels.txt'):
        dele = set(int(x) for x in line.split())
        keep = sorted(set(range(nl)) - dele); idx = {a:i for i,a in enumerate(keep)}
        blocks = [tuple(sorted(idx[a] for a in P if a in idx)) for P in pts]
        blocks = [B for B in blocks if len(B) >= 3]
        key = canon(17, blocks)
        if key not in classes:
            classes[key] = blocks; src[key] = (name, sorted(dele))
        count[key] += 1
print("deletions:", sum(count.values()), "lattice classes:", len(classes))
with open('del17.jsonl', 'w') as out:
    for key, blocks in classes.items():
        mult = sorted((len(B) for B in blocks), reverse=True)
        rec = {"blocks": [list(B) for B in blocks], "source": src[key][0], "deleted": src[key][1], "count": count[key],
               "profile": {str(m): mult.count(m) for m in sorted(set(mult))}}
        out.write(json.dumps(rec) + "\n")
for key in classes:
    b = classes[key]; mult = sorted((len(B) for B in b), reverse=True)
    print(src[key], count[key], {m: mult.count(m) for m in sorted(set(mult))})
