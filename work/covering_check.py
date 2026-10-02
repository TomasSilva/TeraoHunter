# Covering lemma check: does the arrangement A (normals = points of the dual PG(l-1,q)) meet every
# line of the dual space?  (Equivalently: every codim-2 subspace lies in a hyperplane of A.)
import sys, json, itertools
import numpy as np


def proj_points(q, l):
    return [v for v in itertools.product(range(q), repeat=l) if any(v) and next(c for c in v if c) == 1]


def normalize(v, q):
    v = [x % q for x in v]
    c = next(x for x in v if x)
    inv = pow(c, q - 2, q)
    return tuple((x * inv) % q for x in v)


def lines_of(q, l):
    pts = proj_points(q, l)
    seen, out = set(), []
    for a, b in itertools.combinations(pts, 2):
        L = frozenset([a, b] + [normalize([x + t * y for x, y in zip(a, b)], q) for t in range(1, q)])
        if L not in seen:
            seen.add(L)
            out.append(L)
    return out


def is_covering(normals, q):
    l = len(normals[0])
    S = set(normalize(v, q) for v in normals)
    return all(L & S for L in lines_of(q, l))


def private_lines(normals, q):
    """For each h in A*: number of dual lines through h meeting A* only in h (tangents)."""
    l = len(normals[0])
    S = set(normalize(v, q) for v in normals)
    tang = {h: 0 for h in S}
    for L in lines_of(q, l):
        I = L & S
        if len(I) == 1:
            tang[next(iter(I))] += 1
    return tang


if __name__ == '__main__':
    q = 3
    hits = []
    for f in ['work/find_hit_3_4_A.out', 'work/find_hit_3_4_B.out', 'work/find_all_3_4_C.out', 'work/find_all_3_4_D.out']:
        try:
            for line in open('/Users/agent/Projects/TeraoHunter/' + f):
                if line.startswith('HIT'):
                    hits.append(json.loads(line[4:]))
        except FileNotFoundError:
            pass
    for h in hits:
        t = private_lines(h['normals'], q)
        print(h['chunk'], h['k'], h['exps'], 'covering:', is_covering(h['normals'], q),
              'min tangents per point:', min(t.values()))
