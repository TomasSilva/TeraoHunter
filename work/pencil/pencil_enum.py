#!/usr/bin/env python3
"""Three-pencil 18-line lattices with chibar = (t-8)(t-9) and maximal multiplicity 7.
Lines: A_a: x = a y (through P_z=[0:0:1]), B_c: z = c x (through P_y=[0:1:0]), C_b: y = b z (through P_x=[1:0:0]),
plus coordinate lns K subset {x,y,z}.  Triple point A_a ∩ B_c ∩ C_b  iff  a*b*c = 1.
The realization space of such a lattice (abstract triple system T) is the group
  G = {(a,b,c) in (C*)^m : a_i c_j b_k = 1 for (i,j,k) in T} modulo the 2-torus (x,y,z) -> (λx, μy, νz);
components of G <-> torsion of Z^m / (row span of the relation matrix), computed by Smith normal form.
For every non-degenerate component we output a generic point over F_p for the M2 mdr check.
"""
import itertools, json, random, sys, math
from collections import defaultdict

random.seed(12345)
PATTERNS = [
  # name, |A|, |B|, |C|, K, t  (n = |A|+|B|+|C|+|K| = 18 ; sum_centres C(m-1,2) + t = 64)
  ('P1', 5,5,5,'xyz',19),
  ('P2a',6,6,4,'yz',24), ('P2b',6,5,5,'yz',24),
  ('P3a',6,6,5,'x',28), ('P3b',6,5,6,'x',29),
  ('P4', 6,6,6,'',34),
]
# ---------- abelian groups Z^f x prod Z/N: elements are tuples ----------
class Grp:
    def __init__(self, nfree, mods, D):
        self.nfree, self.mods, self.D = nfree, list(mods), D
        self.n = nfree + len(mods)
    def add(self, u, v):
        return tuple(u[i]+v[i] if i < self.nfree else (u[i]+v[i]) % self.mods[i-self.nfree] for i in range(self.n))
    def neg(self, u):
        return tuple(-u[i] if i < self.nfree else (-u[i]) % self.mods[i-self.nfree] for i in range(self.n))
    def box(self):
        ranges = [range(0, self.D+1)]*self.nfree + [range(N) for N in self.mods]
        return [tuple(v) for v in itertools.product(*ranges)]
    def zero(self): return tuple([0]*self.n)
    def name(self):
        return 'Z^%d' % self.nfree + ''.join('xZ/%d' % N for N in self.mods) + (' D=%d' % self.D if self.nfree else '')

def choose_C(rlist, nC, t):
    """rlist: list of (count, elem) sorted by count desc. yield lists of nC elems with sum of counts == t."""
    out = []
    cnts = [c for c,_ in rlist]
    # prefix maxima for pruning
    def rec(start, k, remaining, chosen):
        if k == 0:
            if remaining == 0: out.append(list(chosen))
            return
        for idx in range(start, len(rlist)):
            c = cnts[idx]
            if c > remaining: continue
            # best possible with k items from idx on: c + next (k-1) largest (list sorted desc)
            if sum(cnts[idx:idx+k]) < remaining: break
            chosen.append(rlist[idx][1]); rec(idx+1, k-1, remaining-c, chosen); chosen.pop()
    rec(0, nC, t, [])
    return out

def canonical_cyclic(A, N):
    S = sorted(a[0] for a in A)
    units = [u for u in range(1,N) if math.gcd(u,N)==1]
    best = tuple(S)
    for u in units:
        for v in range(N):
            img = tuple(sorted((u*a+v) % N for a in S))
            if img < best: return False
    return True

def triple_system(G, A, B, C):
    Cidx = {c:k for k,c in enumerate(C)}
    T = set()
    for i,a in enumerate(A):
        for j,b in enumerate(B):
            c = G.neg(G.add(a,b))
            if c in Cidx: T.add((i,j,Cidx[c]))
    return T

def signature(T, nA, nB, nC):
    dA = [0]*nA; dB=[0]*nB; dC=[0]*nC
    for i,j,k in T: dA[i]+=1; dB[j]+=1; dC[k]+=1
    prof = defaultdict(list)
    for i,j,k in T:
        ty = (dA[i], dB[j], dC[k])
        prof[('A',i)].append(ty); prof[('B',j)].append(ty); prof[('C',k)].append(ty)
    sig = []
    for P, n in (('A',nA),('B',nB),('C',nC)):
        sig.append(tuple(sorted(tuple(sorted(prof[(P,i)])) for i in range(n))))
    return tuple(sig)

# ---------- Smith-type diagonalization with column transform U ----------
def diagonalize(R, m):
    """R: list of integer rows (length m). Returns (d, U): d list of m diagonal entries (0 for free),
    U m x m unimodular (list of rows) with rowspan(R U) = rowspan(diag(d))."""
    M = [row[:] for row in R]; t = len(M)
    U = [[1 if i==j else 0 for j in range(m)] for i in range(m)]
    def colop(j, i, q):  # col_j -= q col_i  on M and U
        for row in M: row[j] -= q*row[i]
        for row in U: row[j] -= q*row[i]
    def colswap(i,j):
        for row in M: row[i],row[j] = row[j],row[i]
        for row in U: row[i],row[j] = row[j],row[i]
    k = 0
    while k < min(t,m):
        # find pivot of minimal |value| in submatrix
        best = None
        for i in range(k,t):
            for j in range(k,m):
                v = abs(M[i][j])
                if v and (best is None or v < best[0]): best = (v,i,j)
        if best is None: break
        _,pi,pj = best
        M[k],M[pi] = M[pi],M[k]; colswap(k,pj)
        while True:
            changed = False
            for i in range(k+1,t):
                if M[i][k]:
                    q = M[i][k] // M[k][k]
                    M[i] = [M[i][j] - q*M[k][j] for j in range(m)]
                    if M[i][k]:
                        M[k],M[i] = M[i],M[k]; changed = True
            for j in range(k+1,m):
                if M[k][j]:
                    q = M[k][j] // M[k][k]
                    colop(j,k,q)
                    if M[k][j]:
                        colswap(k,j); changed = True
            if not changed: break
        k += 1
    d = [abs(M[i][i]) if i < t else 0 for i in range(m)]
    return d, U

def components(T, nA, nB, nC, cap=4000):
    m = nA+nB+nC
    rows = []
    for i,j,k in sorted(T):
        r = [0]*m; r[i]=1; r[nA+j]=1; r[nA+nB+k]=1; rows.append(r)
    d, U = diagonalize(rows, m)
    free = [i for i in range(m) if d[i]==0]
    tors = [i for i in range(m) if d[i]>=2]
    # forbidden characters
    forb = []
    for P0,n in ((0,nA),(nA,nB),(nA+nB,nC)):
        for i in range(n):
            for j in range(i+1,n):
                v=[0]*m; v[P0+i]=1; v[P0+j]=-1; forb.append(v)
    for i in range(nA):
        for j in range(nB):
            for k in range(nC):
                if (i,j,k) not in T:
                    v=[0]*m; v[i]=1; v[nA+j]=1; v[nA+nB+k]=1; forb.append(v)
    forbU = []
    for v in forb:
        w = [sum(v[l]*U[l][i] for l in range(m)) for i in range(m)]
        forbU.append(w)
    # characters trivial on G^0 are those with w_i = 0 on free coords
    cand = [w for w in forbU if all(w[i]==0 for i in free)]
    D = 1
    for i in tors: D = D*d[i]//math.gcd(D,d[i])
    ntor = 1
    for i in tors: ntor *= d[i]
    comps = []
    seen = set()
    if ntor <= cap:
        iterator = itertools.product(*[range(d[i]) for i in tors])
    else:
        iterator = (tuple(random.randrange(d[i]) for i in tors) for _ in range(cap))
    units = [u for u in range(1,D+1) if math.gcd(u,D)==1]
    for s in iterator:
        # Galois canonical representative
        reps = []
        for u in units:
            reps.append(tuple((u*s[q]) % d[tors[q]] for q in range(len(tors))))
        key = min(reps)
        if key in seen: continue
        seen.add(key)
        s = key
        degenerate = False
        for w in cand:
            # sum_i s_i w_i / d_i integer ?
            num = 0
            for q,i in enumerate(tors):
                num += s[q]*w[i]*(D//d[i])
            if num % D == 0: degenerate = True; break
        if not degenerate:
            comps.append(s)
    return dict(d=d, U=U, free=free, tors=tors, D=D, ntor=ntor, comps=comps, dimreal=len(free)-2)

def is_prime(n):
    if n<2: return False
    for q in range(2,int(n**0.5)+1):
        if n%q==0: return False
    return True

def find_prime(D, lo=30011):
    p = lo - (lo-1)%D  # p ≡ 1 mod D
    while not is_prime(p): p += D
    return p

def prim_root_of_order(D, p):
    if D == 1: return 1
    fac = [q for q in range(2,D+1) if D%q==0 and is_prime(q)]
    while True:
        h = random.randrange(2,p-1)
        z = pow(h,(p-1)//D,p)
        if all(pow(z, D//q, p) != 1 for q in fac): return z

def realize(T, nA, nB, nC, K, comp, s, p, z):
    """generic point of component s over F_p. returns list of lns [cx,cy,cz]."""
    d, U, free, tors = comp['d'], comp['U'], comp['free'], comp['tors']
    m = nA+nB+nC
    eps = [1]*m
    for q,i in enumerate(tors): eps[i] = pow(z, (comp['D']//d[i])*s[q], p)
    for i in free: eps[i] = random.randrange(2,p-1)
    vals = []
    for l in range(m):
        v = 1
        for i in range(m):
            e = U[l][i]
            if e: v = v*pow(eps[i], e % (p-1), p) % p
        vals.append(v)
    a = vals[:nA]; c = vals[nA:nA+nB]; b = vals[nA+nB:]
    lns = []
    for ai in a: lns.append([1, (-ai)%p, 0])        # x - a y
    for cj in c: lns.append([(-cj)%p, 0, 1])        # z - c x
    for bk in b: lns.append([0, 1, (-bk)%p])        # y - b z
    if 'x' in K: lns.append([1,0,0])
    if 'y' in K: lns.append([0,1,0])
    if 'z' in K: lns.append([0,0,1])
    return lns

def profile(lns, p):
    """multiplicity profile of the arrangement over F_p (points = normalized cross products)."""
    pts = defaultdict(set)
    n = len(lns)
    def norm(v):
        v = [x%p for x in v]
        for i in range(3):
            if v[i]:
                inv = pow(v[i], p-2, p); return tuple(x*inv%p for x in v)
        return None
    for i in range(n):
        for j in range(i+1,n):
            u,v = lns[i],lns[j]
            w = [u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]]
            P = norm(w)
            if P is None: return None  # coincident lns
            pts[P].add(i); pts[P].add(j)
    prof = defaultdict(int)
    for P,S in pts.items(): prof[len(S)] += 1
    b2 = sum((len(S)-1) for S in pts.values())
    return dict(prof), b2

def main():
    groups = []
    wide = '--wide' in sys.argv
    if wide: sys.argv.remove('--wide')
    if wide:
        for N in range(25, 31): groups.append(Grp(0,[N],0))
        groups.append(Grp(1,[7],1)); groups.append(Grp(1,[8],1)); groups.append(Grp(1,[9],1))
        groups.append(Grp(2,[],3)); groups.append(Grp(2,[2],1)); groups.append(Grp(1,[],14))
    else:
      for N in range(3, 25): groups.append(Grp(0,[N],0))
      groups.append(Grp(1,[],10))
      groups.append(Grp(1,[2],5)); groups.append(Grp(1,[3],3)); groups.append(Grp(1,[4],2)); groups.append(Grp(1,[5],1)); groups.append(Grp(1,[6],1))
      groups.append(Grp(2,[],2))
    found = {}   # (pattern, sig) -> list of (T, Gname, A,B,C)
    sel = sys.argv[1:] if len(sys.argv) > 1 else [q[0] for q in PATTERNS]
    suffix = (('_' + '_'.join(sel)) if len(sys.argv) > 1 else '') + ('_wide' if wide else '')
    for (pname,nA,nB,nC,K,t) in PATTERNS:
        if pname not in sel: continue
        for G in groups:
            box = G.box(); zero = G.zero()
            others = [e for e in box if e != zero]
            nfound = 0
            for Arest in itertools.combinations(others, nA-1):
                A = [zero]+list(Arest)
                if G.nfree == 0 and not canonical_cyclic(A, G.mods[0]): continue
                for Brest in itertools.combinations(others, nB-1):
                    B = [zero]+list(Brest)
                    r = defaultdict(int)
                    for a in A:
                        for b in B: r[G.neg(G.add(a,b))] += 1
                    rlist = sorted(((cnt,e) for e,cnt in r.items()), key=lambda x:-x[0])
                    if sum(c for c,_ in rlist[:nC]) < t: continue
                    for C in choose_C(rlist, nC, t):
                        T = triple_system(G,A,B,C)
                        if len(T) != t: continue
                        sig = signature(T,nA,nB,nC)
                        key = (pname, sig)
                        if key not in found: found[key] = []
                        if len(found[key]) < 3:
                            found[key].append((sorted(T), G.name(), A, B, C)); nfound += 1
            print(f'{pname} {G.name():>14}: cumulative signature classes = {sum(1 for k in found if k[0]==pname)}', file=sys.stderr)
    # components and realizations
    entries = []
    for (pname, sig), reps in found.items():
        T, Gname, A, B, C = reps[0]
        pat = [q for q in PATTERNS if q[0]==pname][0]
        _, nA,nB,nC,K,t = pat
        Tset = set(tuple(x) for x in T)
        comp = components(Tset, nA,nB,nC)
        tors = [comp['d'][i] for i in comp['tors']]
        print(f'{pname} sig#{len(entries)} found over {Gname}: dim(real)={comp["dimreal"]} torsion={tors} (|Tor|={comp["ntor"]}) non-degenerate Galois classes: {len(comp["comps"])}', file=sys.stderr)
        p = find_prime(comp['D']); z = prim_root_of_order(comp['D'], p)
        for s in comp['comps']:
            lns = realize(Tset,nA,nB,nC,K,comp,s,p,z)
            pr = profile(lns,p)
            ok = pr is not None and pr[1]==89 and pr[0].get(3,0)==t and max(pr[0])<=7
            entries.append(dict(pattern=pname, sigid=len(entries), T=T, G=Gname, A=[list(a) for a in A], B=[list(b) for b in B], C=[list(c) for c in C],
                                dimreal=comp['dimreal'], torsion=tors, comp=list(s), p=p, lns=lns, profile=(pr[0] if pr else None), b2=(pr[1] if pr else None), ok=ok))
            if not ok: print('   WARNING: profile mismatch', pr, file=sys.stderr)
    with open('pencil_entries%s.json' % suffix,'w') as f: json.dump(entries, f)
    # M2 batch
    with open('pencil_mdr%s.m2' % suffix,'w') as f:
        f.write('entries = {\n')
        for e in entries:
            f.write('  {%d, %s},\n' % (e['p'], str(e['lns']).replace('[','{').replace(']','}')))
        f.write('};\n')
        f.write('''
scan(#ents, idx -> (
  e := ents#idx; p := e#0; kk := ZZ/p; R := kk[x,y,z];
  lns := apply(e#1, v -> v#0*x + v#1*y + v#2*z); Q := product lns; n := #lns;
  jac := matrix{{diff(x,Q), diff(y,Q), diff(z,Q)}};
  S := mingens image syz jac;
  degs := sort apply(numcols S, i -> first (degrees source S)#i - (n-1));
  << "entry " << idx << " p=" << p << " syzdegs=" << degs << " mdr=" << min degs << (if #degs == 2 then " FREE" else " nonfree") << endl;
));
''')
    print(f'{len(entries)} entries written', file=sys.stderr)

if __name__ == '__main__':
    main()
