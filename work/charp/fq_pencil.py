#!/usr/bin/env python3
"""Three-pencil arrangements defined over F_q (A, B, C subsets of F_q^* = Z/(q-1), coordinate lines K) whose
lattice class (d1, d2) has d1 > q: the F_q-rational member has mdr <= q (Frobenius syzygy), hence is non-free.
For each such lattice we compute all realization components over the algebraic closure of F_p (Smith normal form
of the relation lattice; torsion prime to p) and test freeness of a generic point of every component in M2 over GF(p^k).
Any free component is a counterexample to Terao's conjecture over F_p-bar.
Usage: fq_pencil.py q [maxomit]"""
import sys, itertools, math, json, random
from collections import defaultdict
sys.path.insert(0,'.')
from pencil_enum import diagonalize, signature, is_prime
random.seed(2026)
q = int(sys.argv[1]); N = q - 1
p = [r for r in range(2, q+1) if is_prime(r) and q % r == 0][0]
k0 = 1
while p**k0 != q: k0 += 1
maxomit = int(sys.argv[2]) if len(sys.argv) > 2 else q - 3

def lattice_data(A, B, C, K):
    nA, nB, nC = len(A), len(B), len(C)
    n = nA + nB + nC + len(K)
    t = 0
    Cset = set(C)
    for a in A:
        for c in B:
            if (-(a+c)) % N in Cset: t += 1
    mz = nA + ('x' in K) + ('y' in K); my = nB + ('x' in K) + ('z' in K); mx = nC + ('y' in K) + ('z' in K)
    b2 = n*(n-1)//2 - sum((m-1)*(m-2)//2 for m in (mz,my,mx)) - t
    disc = (n-1)**2 - 4*(b2 - n + 1)
    if disc < 0: return None
    s = int(round(disc**0.5))
    if s*s != disc or (n-1-s) % 2: return None
    d1 = (n-1-s)//2; d2 = (n-1+s)//2
    return dict(n=n, t=t, b2=b2, d1=d1, d2=d2, m=(mz,my,mx))

slots = [('A',i) for i in range(N)] + [('B',i) for i in range(N)] + [('C',i) for i in range(N)]
cands = {}
count = 0
for K in ['', 'x','y','z','xy','xz','yz','xyz']:
    for r in range(0, maxomit+1):
        if 3*N - r + len(K) < 2*q + 3: continue
        for O in itertools.combinations(range(3*N), r):
            Oset = set(O)
            A = [i for i in range(N) if i not in Oset]
            B = [i for i in range(N) if (N+i) not in Oset]
            C = [i for i in range(N) if (2*N+i) not in Oset]
            if not A or not B or not C: continue
            d = lattice_data(A,B,C,K)
            count += 1
            if d is None or d['d1'] <= q or d['n'] % p == 0: continue
            # triple system
            T = set()
            Cidx = {c:kk for kk,c in enumerate(C)}
            for i,a in enumerate(A):
                for j,c in enumerate(B):
                    b = (-(a+c)) % N
                    if b in Cidx: T.add((i,j,Cidx[b]))
            sig = (K, signature(T, len(A), len(B), len(C)), d['n'], d['b2'], d['d1'])
            if sig not in cands:
                cands[sig] = dict(A=A,B=B,C=C,K=K,T=sorted(T), **d)
print(f'q={q} p={p}: {count} (A,B,C,K) examined, {len(cands)} candidate lattices with d1 > q', file=sys.stderr)

def components_charp(T, nA, nB, nC, p, cap=3000):
    m = nA+nB+nC
    rows = []
    for i,j,k in sorted(T):
        r = [0]*m; r[i]=1; r[nA+j]=1; r[nA+nB+k]=1; rows.append(r)
    r=[0]*m; r[0]=1; rows.append(r)
    r=[0]*m; r[nA]=1; rows.append(r)
    d, U = diagonalize(rows, m)
    # prime-to-p parts of the torsion
    dp = []
    for x in d:
        if x >= 2:
            while x % p == 0: x //= p
        dp.append(x)
    free = [i for i in range(m) if dp[i]==0]
    tors = [i for i in range(m) if dp[i]>=2]
    forb = []
    for P0,nn in ((0,nA),(nA,nB),(nA+nB,nC)):
        for i in range(nn):
            for j in range(i+1,nn):
                v=[0]*m; v[P0+i]=1; v[P0+j]=-1; forb.append(v)
    for i in range(nA):
        for j in range(nB):
            for k in range(nC):
                if (i,j,k) not in T:
                    v=[0]*m; v[i]=1; v[nA+j]=1; v[nA+nB+k]=1; forb.append(v)
    forbU = [[sum(v[l]*U[l][i] for l in range(m)) for i in range(m)] for v in forb]
    D = 1
    for i in tors: D = D*dp[i]//math.gcd(D,dp[i])
    ntor = 1
    for i in tors: ntor *= dp[i]
    comps = []
    it = itertools.product(*[range(dp[i]) for i in tors]) if ntor <= cap else (tuple(random.randrange(dp[i]) for i in tors) for _ in range(cap))
    seen = set()
    # Galois classes under Frobenius (multiplication by p)
    for s in it:
        orbit = []
        u = 1
        for _ in range(60):
            orbit.append(tuple((u*s[r]) % dp[tors[r]] for r in range(len(tors)))); u = (u*p) % D if D > 1 else 1
        key = min(orbit)
        if key in seen: continue
        seen.add(key); s = key
        degenerate = False
        for w in forbU:
            if all(w[i]==0 for i in free):
                num = sum(s[r]*w[i]*(D//dp[i]) for r,i in enumerate(tors))
                if num % D == 0: degenerate = True; break
        if not degenerate: comps.append(s)
    return dict(d=dp, U=U, free=free, tors=tors, D=D, ntor=ntor, comps=comps, m=m)

def field_k(D, p, k0):
    # smallest k multiple of k0 with D | p^k - 1 and p^k >= 10^4
    k = k0
    while (p**k - 1) % D != 0 or p**k < 10**4: k += k0
    return k

out = []
m2 = []
for idx, (sig, c) in enumerate(cands.items()):
    nA,nB,nC = len(c['A']),len(c['B']),len(c['C'])
    T = set(tuple(x) for x in c['T'])
    comp = components_charp(T, nA,nB,nC, p)
    k = field_k(comp['D'], p, k0)
    c['components'] = comp['comps']; c['torsion'] = [comp['d'][i] for i in comp['tors']]; c['dim'] = len(comp['free']); c['k'] = k
    print(f"cand {idx}: n={c['n']} class=({c['d1']},{c['d2']}) m={c['m']} t={c['t']} |A|,|B|,|C|={nA},{nB},{nC} K={c['K']!r} dim={c['dim']} torsion={c['torsion']} components={len(comp['comps'])} field GF({p}^{k})", file=sys.stderr)
    # M2 block per component
    for ci, s in enumerate(comp['comps']):
        lines = []
        def lau(l):
            num, den, const = [], [], []
            for i in range(comp['m']):
                ex = comp['U'][l][i]
                if ex == 0: continue
                if i in comp['free']:
                    (num if ex>0 else den).append(f'g_{comp["free"].index(i)}^{abs(ex)}')
                elif comp['d'][i] >= 2:
                    r = comp['tors'].index(i)
                    e = ((p**k - 1)//comp['d'][i]) * s[r] * ex
                    const.append(f'w^({e % (p**k-1)})')
            cst = '*'.join(const) if const else '1'
            return cst, ('*'.join(num) if num else '1'), ('*'.join(den) if den else '1')
        for l in range(nA):
            cst,nu,de = lau(l); lines.append(f'({de})*x - ({cst}*{nu})*y')
        for l in range(nA, nA+nB):
            cst,nu,de = lau(l); lines.append(f'({de})*z - ({cst}*{nu})*x')
        for l in range(nA+nB, nA+nB+nC):
            cst,nu,de = lau(l); lines.append(f'({de})*y - ({cst}*{nu})*z')
        for v in c['K']: lines.append(v)
        gdefs = ' '.join(f'g_{i} = w^({random.randrange(1, p**k-1)});' for i in range(len(comp['free'])))
        m2.append(f'''
<< "cand {idx} comp {ci} s={list(s)}: n={c['n']} class=({c['d1']},{c['d2']}) dim={c['dim']} " << endl;
kk = GF({p},{k}); w = kk_0; R = kk[x,y,z]; {gdefs}
lns = {{{', '.join(lines)}}};
if #unique lns != #lns then << "   COINCIDENT LINES" << endl;
Q = product lns; degs = mdrOf Q;
<< "   syz degrees " << degs << " mdr=" << min degs << (if #degs == 2 then "  FREE  <== COUNTEREXAMPLE?" else "  nonfree") << endl;
''')
    out.append(c)
json.dump(out, open(f'fq_cands_q{q}.json','w'))
with open(f'fq_mdr_q{q}.m2','w') as f:
    f.write('''mdrOf = Q -> (R := ring Q; jac := matrix{{diff(R_0,Q),diff(R_1,Q),diff(R_2,Q)}}; S := mingens image syz jac;
   sort apply(numcols S, i -> first (degrees source S)#i - (first degree Q - 1)));
''')
    f.write('\n'.join(m2))
print(f'{len(out)} candidates, {len(m2)} component tests written to fq_mdr_q{q}.m2', file=sys.stderr)
