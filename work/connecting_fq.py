# Connecting-line arrangements over ANY small field F_q (prime powers included), exhaustive over
# k-point sets up to PGammaL(3,q).  For free non-DF ones: all realizations of the lattice over F_q.
import sys, itertools, time, json
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import pynauty
from hunter.fq import GF, LineArr, realizations

q = int(sys.argv[1]); KMAX = int(sys.argv[2])
F = GF(q)
pts = F.points(3)
NP = len(pts)
inc = [[j for j, P in enumerate(pts) if F.dot(L, P) == 0] for L in pts]
adj = {NP + i: inc[i] for i in range(NP)}


def cert(P):
    P = set(P)
    g = pynauty.Graph(2 * NP, directed=False, adjacency_dict=adj,
                      vertex_coloring=[set(P), set(j for j in range(NP) if j not in P), set(range(NP, 2 * NP))])
    return pynauty.certificate(g)


levels = {0: {cert([]): ()}}
out = open(f'/Users/agent/Projects/TeraoHunter/work/connecting_fq_{q}.jsonl', 'w')
t0 = time.time()
for k in range(1, KMAX + 1):
    nxt = {}
    for P in levels[k - 1].values():
        for j in range(NP):
            if j in P:
                continue
            T = tuple(sorted(P + (j,)))
            c = cert(T)
            if c not in nxt:
                nxt[c] = T
    levels[k] = nxt
    st = {'orbits': len(nxt), 'chi': 0, 'free': 0, 'nondf_free': 0, 'multi_real': 0, 'mixed': 0}
    if k >= 4:
        for P in nxt.values():
            Pset = set(P)
            Ls = [pts[i] for i in range(NP) if sum(1 for j in inc[i] if j in Pset) >= 2]
            if len(Ls) < 4:
                continue
            try:
                A = LineArr(F, Ls)
            except ValueError:
                continue
            if len(A.flats) == 1:
                continue                     # pencil
            if A.exponents() is None:
                continue
            st['chi'] += 1
            if not A.free():
                continue
            st['free'] += 1
            if A.is_df():
                continue
            st['nondf_free'] += 1
            flats2 = [sorted(S) for S in A.flats]
            sols = realizations(F, A.n, flats2, max_solutions=5000)
            frs = [LineArr(F, s).free() for s in sols]
            nf = len(frs) - sum(frs)
            if len(sols) > 1:
                st['multi_real'] += 1
            rec = {'q': q, 'k': k, 'n_lines': A.n, 'exps': A.exponents(),
                   't': sorted(len(S) for S in A.flats if len(S) > 2), 'realizations': len(sols), 'nonfree': nf,
                   'P': list(P)}
            out.write(json.dumps(rec) + '\n'); out.flush()
            if nf and nf < len(frs):
                st['mixed'] += 1
                print('!!! MIXED FREENESS (TC COUNTEREXAMPLE OVER GF(%d)):' % q, json.dumps(rec), flush=True)
    print('k=%d %s  (%.0fs)' % (k, st, time.time() - t0), flush=True)
