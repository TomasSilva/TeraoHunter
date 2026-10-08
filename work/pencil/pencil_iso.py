#!/usr/bin/env python3
"""Exact isomorphism classes of the triple systems found by pencil_enum.py (per pattern), including pencil permutations
allowed by K (P1: all of S3; P2a: swap A<->B)."""
import json, sys, itertools
from collections import defaultdict
entries = json.load(open(sys.argv[1]))
def canon(T, nA, nB, nC):
    best = None
    for sa in itertools.permutations(range(nA)):
        for sb in itertools.permutations(range(nB)):
            M = {}
            for i,j,k in T: M[(sa[i], sb[j])] = k
            relabel = {}; rows = []
            for i in range(nA):
                row = []
                for j in range(nB):
                    k = M.get((i,j))
                    if k is None: row.append(-1)
                    else:
                        if k not in relabel: relabel[k] = len(relabel)
                        row.append(relabel[k])
                rows.append(tuple(row))
            key = tuple(rows)
            if best is None or key < best: best = key
    return best
def canon_full(T, nA, nB, nC, perms):
    best = None
    for perm in perms:
        sizes = [nA,nB,nC]
        T2 = [tuple(t[perm[q]] for q in range(3)) for t in T]
        n2 = [sizes[perm[q]] for q in range(3)]
        key = (tuple(n2), canon(T2, *n2))
        if best is None or key < best: best = key
    return best
groups = defaultdict(list)
for e in entries:
    groups[e['pattern']].append(e)
for pname, es in groups.items():
    perms = list(itertools.permutations(range(3))) if pname == 'P1' else [(0,1,2),(1,0,2)]
    classes = defaultdict(list)
    seen = set()
    for e in es:
        key = (e['sigid'], tuple(map(tuple, e['T'])))
        if key in seen: continue
        seen.add(key)
        T = [tuple(t) for t in e['T']]
        nA, nB, nC = len(e['A']), len(e['B']), len(e['C'])
        c = canon_full(T, nA, nB, nC, perms)
        classes[c].append((e['sigid'], e['G'], e['dimreal'], e['torsion']))
    print(f'{pname}: {len(classes)} isomorphism classes of lattices among {len(seen)} (sigid,T) pairs')
    for c, members in classes.items():
        print('   class:', sorted(set(m[0] for m in members)), 'dim', members[0][2], 'torsion', members[0][3], 'found over', sorted(set(m[1] for m in members)))
