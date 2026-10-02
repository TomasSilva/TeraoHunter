# Exhaustive Terao check over GF(q) (q prime) in rank 3:
#   all subsets of lines of PG(2,q) up to PGL(3,q)  (orderly generation with pynauty
#   canonical certificates of the coloured incidence graph),
#   lattice isomorphism class (pynauty certificate of the line / flat bipartite graph),
#   exact freeness over F_q (Ziegler-restriction criterion, any characteristic).
# Reports lattice classes realised by several PGL-orbits with different freeness.
import sys, itertools, time, json, pickle
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
import pynauty
from hunter.core import Lattice
from hunter.ffree import free_rank3_ziegler

q = int(sys.argv[1])
KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else None

pts = [v for v in itertools.product(range(q), repeat=3) if any(v) and next(c for c in v if c) == 1]
lines = list(pts)                         # self-dual coordinates
NP = len(pts)
inc = [[j for j, P in enumerate(pts) if sum(a * b for a, b in zip(L, P)) % q == 0] for L in lines]
print('PG(2,%d): %d points, %d lines, %d points per line' % (q, NP, NP, len(inc[0])), flush=True)

# coloured incidence graph: vertices 0..NP-1 points, NP..2NP-1 lines
adj = {NP + i: [j for j in inc[i]] for i in range(NP)}


def cert(S):
    S = set(S)
    g = pynauty.Graph(2 * NP, directed=False, adjacency_dict=adj,
                      vertex_coloring=[set(range(NP)),
                                       set(NP + i for i in S),
                                       set(NP + i for i in range(NP) if i not in S)])
    return pynauty.certificate(g)


def lattice_cert(S):
    S = sorted(S)
    flats = []
    for j in range(NP):
        through = [i for i in S if j in inc_set[i]]
        if len(through) >= 2:
            flats.append(through)
    idx = {L: t for t, L in enumerate(S)}
    nv = len(S) + len(flats)
    ad = {len(S) + f: [idx[L] for L in F] for f, F in enumerate(flats)}
    g = pynauty.Graph(nv, directed=False, adjacency_dict=ad,
                      vertex_coloring=[set(range(len(S))), set(range(len(S), nv))] if flats else [set(range(len(S)))])
    return pynauty.certificate(g), flats


inc_set = [set(x) for x in inc]
levels = {0: {cert([]): ()}}
t0 = time.time()
kmax = KMAX if KMAX is not None else NP
for k in range(1, kmax + 1):
    nxt = {}
    for S in levels[k - 1].values():
        for L in range(NP):
            if L in S:
                continue
            T = tuple(sorted(S + (L,)))
            c = cert(T)
            if c not in nxt:
                nxt[c] = T
    levels[k] = nxt
    print('k=%d orbits=%d  (%.0fs)' % (k, len(nxt), time.time() - t0), flush=True)
    if k > NP // 2 and KMAX is None:
        pass
pickle.dump({k: list(v.values()) for k, v in levels.items()},
            open(f'/Users/agent/Projects/TeraoHunter/work/pg2_{q}_orbits.pkl', 'wb'))

# lattices and freeness
classes = {}
nfree = 0
for k in range(3, kmax + 1):
    for S in levels[k].values():
        N = np.array([lines[i] for i in S], dtype=np.int64)
        if np.linalg.matrix_rank(N.astype(float)) < 3:
            continue
        lc, flats = lattice_cert(S)
        try:
            L = Lattice(N % q, q)
        except ValueError:
            continue
        st = free_rank3_ziegler(N, q, L)
        if st == 'free':
            nfree += 1
        classes.setdefault((k, lc), []).append((S, st))
mixed = {key: v for key, v in classes.items() if len(set(s for _, s in v)) > 1}
print('lattice classes:', len(classes), 'free orbits:', nfree,
      'classes with >1 orbit:', sum(1 for v in classes.values() if len(v) > 1), flush=True)
print('MIXED-FREENESS CLASSES (counterexamples over GF(%d)):' % q, len(mixed), flush=True)
for (k, lc), v in mixed.items():
    print('  n=%d' % k, [(S, s) for S, s in v], flush=True)
json.dump({'q': q, 'classes': len(classes), 'mixed': [[k, [(list(S), s) for S, s in v]] for (k, lc), v in mixed.items()]},
          open(f'/Users/agent/Projects/TeraoHunter/work/pg2_{q}_result.json', 'w'))
