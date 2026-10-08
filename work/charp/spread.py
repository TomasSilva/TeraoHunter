#!/usr/bin/env python3
"""One-parameter three-pencil families over Z (multiplicative: x = g^a y, z = g^c x, y = g^b z, K coord lines, triple iff a+c+b=0;
additive: y = a z, x = c z, x - y = b z, optional z, triple iff c - a = b) with integer-rooted chibar.
For each family: forbidden integers W (differences within a set, and a+c+b resp. c-a-b over non-triples); the reduction g -> element of
order N (mult.) resp. parameters mod p (add.) keeps the lattice iff N (resp. p) divides no w in W.  A prime power q = N+1 (resp. prime p)
with q < d1 and p not dividing n would make the F_q-point non-free (Frobenius); the family is then a counterexample iff its generic
member is free.  We list all families admitting such q, with the generic syzygy degrees computed in M2 (over QQ)."""
import itertools, math, sys, json
D = int(sys.argv[1]) if len(sys.argv) > 1 else 8
def isprimepower(q):
    if q < 2: return None
    for r in range(2, q+1):
        if q % r == 0:
            if all(r % s for s in range(2, r)):
                x = q
                while x % r == 0: x //= r
                return r if x == 1 else None
            return None
def roots(n, b2):
    disc = (n-1)**2 - 4*(b2-n+1)
    if disc < 0: return None
    s = math.isqrt(disc)
    if s*s != disc or (n-1-s) % 2: return None
    return (n-1-s)//2, (n-1+s)//2
out = []
def consider(kind, A, B, C, K, n, t, mults, Wset):
    b2 = n*(n-1)//2 - sum((m-1)*(m-2)//2 for m in mults) - t
    r = roots(n, b2)
    if r is None: return
    d1, d2 = r
    W = set(abs(w) for w in Wset if w != 0)
    # admissible q: multiplicative: N = q-1 divides no w; additive: p = q divides no w (q prime)
    adm = []
    for q in range(3, d1):
        pp = isprimepower(q)
        if pp is None or n % pp == 0: continue
        N = q - 1 if kind == 'mult' else q
        if kind == 'add' and pp != q: continue
        if all(w % N for w in W): adm.append(q)
    if adm:
        out.append(dict(kind=kind, A=A, B=B, C=C, K=K, n=n, t=t, d1=d1, d2=d2, q=adm, maxW=max(W) if W else 0))
rng = range(0, D+1)
for sA in range(3, 8):
    for Arest in itertools.combinations(range(1, D+1), sA-1):
        A = (0,)+Arest
        if max(A) != D and sA < 8: pass
        for sB in range(3, 8):
            for Brest in itertools.combinations(range(1, D+1), sB-1):
                B = (0,)+Brest
                # multiplicative: sums a+c; C from -(A+B): choose C = top-r sums (all tie-breaks), r in 3..7
                from collections import Counter
                cnt = Counter(a+c for a in A for c in B)
                vals = sorted(cnt.items(), key=lambda x: (-x[1], x[0]))
                for r in range(3, 8):
                    if r > len(vals): break
                    thr = vals[r-1][1]
                    above = [v for v,c in vals if c > thr]; ties = [v for v,c in vals if c == thr]
                    need = r - len(above)
                    for tie in itertools.combinations(ties, need):
                        Csum = sorted(above + list(tie))        # the sums s with a+c = s carrying triple points
                        # multiplicative: b = -s ; additive: with a' = a, c' = c, b = c - a -> use differences instead (handled below)
                        t = sum(cnt[s] for s in Csum)
                        Cm = [-s for s in Csum]
                        # forbidden: differences within A, B, Cm; a + c + b for non-triples
                        W = set()
                        for X in (A, B, Cm):
                            for u,v in itertools.combinations(X, 2): W.add(u-v)
                        for a in A:
                            for c in B:
                                for b in Cm:
                                    if a + c + b != 0: W.add(a+c+b)
                        for K in ('', 'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz'):
                            n = sA + sB + r + len(K)
                            mz = sA + ('x' in K) + ('y' in K); my = sB + ('x' in K) + ('z' in K); mx = r + ('y' in K) + ('z' in K)
                            consider('mult', list(A), list(B), Cm, K, n, t, (mz,my,mx), W)
                # additive: triple iff c - a = b : counts of differences
                cntd = Counter(c - a for a in A for c in B)
                valsd = sorted(cntd.items(), key=lambda x: (-x[1], x[0]))
                for r in range(3, 8):
                    if r > len(valsd): break
                    thr = valsd[r-1][1]
                    above = [v for v,c in valsd if c > thr]; ties = [v for v,c in valsd if c == thr]
                    need = r - len(above)
                    for tie in itertools.combinations(ties, need):
                        Cd = sorted(above + list(tie)); t = sum(cntd[s] for s in Cd)
                        W = set()
                        for X in (A, B, Cd):
                            for u,v in itertools.combinations(X, 2): W.add(u-v)
                        for a in A:
                            for c in B:
                                for b in Cd:
                                    if c - a - b != 0: W.add(c-a-b)
                        for zl in (0, 1):
                            n = sA + sB + r + zl
                            consider('add', list(A), list(B), Cd, 'z' if zl else '', n, t, (sA+zl, sB+zl, r+zl), W)
# dedupe
seen = set(); res = []
for o in out:
    key = (o['kind'], tuple(o['A']), tuple(o['B']), tuple(o['C']), o['K'])
    if key in seen: continue
    seen.add(key); res.append(o)
print(f'{len(res)} families with an admissible q < d1', file=sys.stderr)
res.sort(key=lambda o: (o['n'], o['d1']))
json.dump(res, open(f'spread_D{D}.json','w'))
for o in res[:60]:
    print(o['kind'], 'n', o['n'], 'class', (o['d1'], o['d2']), 'A', o['A'], 'B', o['B'], 'C', o['C'], 'K', o['K'], 't', o['t'], 'q', o['q'], 'maxW', o['maxW'])
