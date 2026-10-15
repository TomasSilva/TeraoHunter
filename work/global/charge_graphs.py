"""
Count jump points of a 4-point multiarrangement on P^1 via critically fixed rational maps.

Claim (TC_GLOBAL.md, Theorem B): for Ziegler multiplicities k = (k1,..,k4), K = sum k even and
balanced, the jump points (cross-ratios where e1 drops from K/2 to K/2 - 1) correspond bijectively
to pairs (Phi, ordering), where Phi is a critically fixed rational map of degree D = K/2 - 1 whose
critical points are exactly the points with k_i >= 2 (local degree k_i), and the points with
k_i = 1 are distinct non-critical fixed points of Phi, chosen in order.

Critically fixed maps up to Moebius conjugacy <-> connected planar loopless multigraphs (plane maps)
with D - 1 edges, a vertex of valence v being a critical point of local degree v + 1
(Pilgrim-Tan blow-up; Cordwell-Gilbertson-Nuechterlein-Pilgrim-Pinella 2015; Hlushchanka 2019).
With labelled vertices, plane maps are counted up to isomorphism fixing the labels.

Prediction: #jump points = #plane maps * (D+1-c)(D-c)...(D-c-(4-c)+1), c = #{k_i >= 2}.
This script computes the prediction and compares with the interior degree of the jump polynomial
in jump4_all.out (boundary factors mu^a (mu-1)^b removed).
"""
import itertools, re, sys
from math import prod

def multigraphs(val):
    """All edge-multiplicity dicts on vertices 0..c-1 with given valences, connected, loopless."""
    c = len(val)
    pairs = [(i, j) for i in range(c) for j in range(i + 1, c)]
    res = []
    def rec(idx, rem, cur):
        if idx == len(pairs):
            if all(r == 0 for r in rem):
                res.append(dict(cur))
            return
        i, j = pairs[idx]
        # remaining pairs involving i after this one
        for m in range(min(rem[i], rem[j]) + 1):
            rem[i] -= m; rem[j] -= m
            cur[(i, j)] = m
            # prune: vertex i has no later pairs with smaller index i? check feasibility cheaply
            ok = True
            if ok:
                rec(idx + 1, rem, cur)
            rem[i] += m; rem[j] += m
        del cur[(i, j)]
    rec(0, list(val), {})
    out = []
    for g in res:
        # connectivity
        adj = {v: set() for v in range(c)}
        for (i, j), m in g.items():
            if m:
                adj[i].add(j); adj[j].add(i)
        seen = {0}; st = [0]
        while st:
            v = st.pop()
            for w in adj[v]:
                if w not in seen:
                    seen.add(w); st.append(w)
        if len(seen) == c:
            out.append(g)
    return out

def plane_maps(val, g):
    """Number of genus-0 rotation systems of multigraph g (labelled vertices) up to isomorphism."""
    c = len(val)
    darts_at = {v: [] for v in range(c)}
    alpha = {}
    eid = 0
    for (i, j), m in sorted(g.items()):
        for _ in range(m):
            a, b = (eid, 0), (eid, 1)
            darts_at[i].append(a); darts_at[j].append(b)
            alpha[a] = b; alpha[b] = a
            eid += 1
    E = eid
    if E == 0:
        return 1 if c == 1 else 0
    vert = {}
    for v in range(c):
        for d in darts_at[v]:
            vert[d] = v
    codes = set()
    orders = []
    for v in range(c):
        ds = darts_at[v]
        if len(ds) <= 1:
            orders.append([tuple(ds)])
        else:
            first = ds[0]
            orders.append([(first,) + p for p in itertools.permutations(ds[1:])])
    for combo in itertools.product(*orders):
        sigma = {}
        for cyc in combo:
            for k, d in enumerate(cyc):
                sigma[d] = cyc[(k + 1) % len(cyc)]
        # faces of phi = sigma o alpha
        seen = set(); F = 0
        for d in alpha:
            if d not in seen:
                F += 1
                x = d
                while x not in seen:
                    seen.add(x)
                    x = sigma[alpha[x]]
        if c - E + F != 2:
            continue
        # canonical code: BFS from each dart at vertex 0
        best = None
        for d0 in darts_at[0]:
            lab = {d0: 0}; q = [d0]; code = []
            k = 0
            while k < len(q):
                d = q[k]; k += 1
                for nb in (sigma[d], alpha[d]):
                    if nb not in lab:
                        lab[nb] = len(q); q.append(nb)
                code.append((vert[d], lab[sigma[d]], lab[alpha[d]]))
            code = tuple(code)
            if best is None or code < best:
                best = code
        codes.add(best)
    return len(codes)

def predicted(k):
    K = sum(k)
    if K % 2:
        return 0
    D = K // 2 - 1
    val = [x - 1 for x in k if x >= 2]
    c = len(val)
    if c == 0:
        return 0
    if c == 1:
        return 0  # a single critical point is impossible (Riemann-Hurwitz needs >= 2)
    nmaps = sum(plane_maps(val, g) for g in multigraphs(val))
    free = D + 1 - c          # non-critical fixed points of a critically fixed map
    need = 4 - c
    choices = prod(range(free - need + 1, free + 1)) if need > 0 else 1
    if c == 2:
        # the map is z^D; its Moebius automorphisms fixing both critical points are z -> eta z,
        # eta^(D-1) = 1, which permute the ordered choices freely
        return nmaps * choices // (D - 1)
    return nmaps * choices

def interior_degree(line):
    poly = line.split('|', 1)[1].strip()
    deg = 0
    for f in re.findall(r'\(([^()]*)\)(?:\^(\d+))?', poly):
        expr, e = f
        e = int(e) if e else 1
        if 'mu' not in expr:
            continue
        if expr.strip() in ('mu', 'mu-1'):
            continue
        m = re.findall(r'mu\^(\d+)', expr)
        d = max([int(x) for x in m] + [1])
        deg += d * e
    return deg

if __name__ == '__main__':
    Kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    bad = 0; tot = 0
    for line in open('/Users/agent/Projects/TeraoHunter/work/rational/jump4_all.out'):
        m = re.match(r'K (\d+) k \{([^}]*)\}', line)
        if not m:
            continue
        K = int(m.group(1)); k = [int(x) for x in m.group(2).split(',')]
        if K > Kmax or K % 2:
            continue
        p = predicted(k); d = interior_degree(line)
        tot += 1
        flag = '' if p == d else '   <-- MISMATCH'
        if p != d:
            bad += 1
        print(K, k, 'graphs-predicted', p, 'jump-poly', d, flag, flush=True)
    print('vectors', tot, 'mismatches', bad)
