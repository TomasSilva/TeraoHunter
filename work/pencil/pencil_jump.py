#!/usr/bin/env python3
"""For each three-pencil lattice found by pencil_enum.py: normalized realization group G' = G ∩ {a_0 = c_0 = 1} = G/T2,
all non-degenerate Galois classes of components; for each component an M2 block that
 (1) computes mdr at a random point of the component, and
 (2) if the component has positive dimension, computes the non-free locus inside it by Yoshinaga's criterion:
     A free  <=>  the Ziegler multirestriction on a line H has exponents (8,9)  <=>  the 17 x 16 matrix of D(m)_7 has full rank;
     the locus V(minors) is saturated by the degeneration locus (forbidden characters).
Usage: pencil_jump.py entries.json  ->  writes jump_<sigid>.m2 for every lattice (one per signature id group)."""
import json, sys, random, math, itertools
sys.path.insert(0, '.')
from pencil_enum import diagonalize, find_prime, prim_root_of_order, is_prime, PATTERNS

random.seed(777)

def normalized_components(T, nA, nB, nC, cap=4000):
    m = nA+nB+nC
    rows = []
    for i,j,k in sorted(T):
        r=[0]*m; r[i]=1; r[nA+j]=1; r[nA+nB+k]=1; rows.append(r)
    r=[0]*m; r[0]=1; rows.append(r)          # a_0 = 1
    r=[0]*m; r[nA]=1; rows.append(r)         # c_0 = 1
    d, U = diagonalize(rows, m)
    free = [i for i in range(m) if d[i]==0]
    tors = [i for i in range(m) if d[i]>=2]
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
    forbU = [[sum(v[l]*U[l][i] for l in range(m)) for i in range(m)] for v in forb]
    D = 1
    for i in tors: D = D*d[i]//math.gcd(D,d[i])
    ntor = 1
    for i in tors: ntor *= d[i]
    units = [u for u in range(1,D+1) if math.gcd(u,D)==1]
    comps, seen = [], set()
    it = itertools.product(*[range(d[i]) for i in tors]) if ntor <= cap else (tuple(random.randrange(d[i]) for i in tors) for _ in range(cap))
    for s in it:
        key = min(tuple((u*s[q]) % d[tors[q]] for q in range(len(tors))) for u in units)
        if key in seen: continue
        seen.add(key); s = key
        degenerate = False
        for w in forbU:
            if all(w[i]==0 for i in free):
                num = sum(s[q]*w[i]*(D//d[i]) for q,i in enumerate(tors))
                if num % D == 0: degenerate = True; break
        if not degenerate: comps.append(s)
    return dict(d=d,U=U,free=free,tors=tors,D=D,ntor=ntor,comps=comps,forbU=forbU,m=m)

def m2_script(T, nA, nB, nC, K, comp, tag):
    d,U,free,tors,D,m = comp['d'],comp['U'],comp['free'],comp['tors'],comp['D'],comp['m']
    p = find_prime(D, lo=30011); z = prim_root_of_order(D, p)
    e = len(free)
    out = []
    out.append(f'-- lattice {tag}: |A|,|B|,|C| = {nA},{nB},{nC}, K = {K!r}, t = {len(T)}; dim(realization component) = {e}; torsion = {[d[i] for i in tors]}')
    out.append(f'p = {p}; kk = ZZ/p; zD = {z}_kk;')
    gvars = ', '.join(f'g_{i}' for i in range(e)) if e else ''
    out.append('R0 = kk[x,y,z];')
    out.append(f'P = kk[{gvars}];' if e else 'P = kk;')
    out.append('R = P[x,y,z];')
    out.append(r"""
mdrOf = Q -> (jac := matrix{{diff((ring Q)_0,Q),diff((ring Q)_1,Q),diff((ring Q)_2,Q)}}; S := mingens image syz jac;
   sort apply(numcols S, i -> first (degrees source S)#i - (first degree Q - 1)));
-- Yoshinaga: A is free iff the Ziegler multirestriction on H has exponents (8,9) iff the 17x16 matrix of D(m)_7 has full rank.
-- lns = linear forms over R = P[x,y,z]; h = index of H (a coordinate line); spec: P -> kk numerical specialization for grouping.
zieglerMatrix = (lns, h, spec) -> (
   H := lns#h; others := delete(H, lns);
   vars3 := {x,y,z};
   k := position(vars3, v -> coefficient(v, H) != 0);
   vk := vars3#k; rest := delete(vk, vars3);
   invc := 1/(lift(coefficient(vk, H), kk));
   forms := apply(others, L -> sub(L, {vk => -(H - coefficient(vk,H)*vk) * invc}));
   U := rest#0; V := rest#1;
   pq := apply(forms, f -> (coefficient(U,f), coefficient(V,f)));
   pqn := apply(pq, w -> (spec(w#0), spec(w#1)));
   groups := {};
   scan(#pq, i -> (
      j := position(groups, G -> (w := pqn#(G#0); w#0*pqn#i#1 - w#1*pqn#i#0 == 0));
      if j === null then groups = append(groups, {i}) else groups = replace(j, append(groups#j, i), groups)));
   mults := apply(groups, G -> #G);
   S := P[symbol u, symbol v]; u := S_0; v := S_1;
   mons := flatten entries basis(7, S);
   rowsL := {};
   scan(groups, G -> (
      pp := pq#(G#0)#0; qq := pq#(G#0)#1; mu := #G;
      for i from 0 to mu-1 do (
         der := u^i * v^(mu-1-i);
         rowsL = append(rowsL, flatten {
            apply(mons, mo -> lift(sub(pp * diff(der, mo), {u => qq, v => -pp}), P)),
            apply(mons, mo -> lift(sub(qq * diff(der, mo), {u => qq, v => -pp}), P)) });
      )));
   (matrix rowsL, mults));
""")
    def laurent(l, s):
        num, den = [], []
        const = 1
        for i in range(m):
            ex = U[l][i]
            if ex == 0: continue
            if i in free:
                gi = f'g_{free.index(i)}'
                (num if ex>0 else den).append(f'{gi}^{abs(ex)}')
            elif d[i] >= 2:
                q = tors.index(i); const = (const * pow(z, ((D//d[i])*s[q]*ex) % D, p)) % p
        return const, ('*'.join(num) if num else '1'), ('*'.join(den) if den else '1')
    blocks = []
    for ci, s in enumerate(comp['comps']):
        L = []
        for l in range(nA):
            c,nu,de = laurent(l, s); L.append(f'({de})*x - ({c}*{nu})*y')
        for l in range(nA, nA+nB):
            c,nu,de = laurent(l, s); L.append(f'({de})*z - ({c}*{nu})*x')
        for l in range(nA+nB, m):
            c,nu,de = laurent(l, s); L.append(f'({de})*y - ({c}*{nu})*z')
        if 'x' in K: L.append('x')
        if 'y' in K: L.append('y')
        if 'z' in K: L.append('z')
        degs = []
        for w in comp['forbU']:
            wf = [w[i] for i in free]
            if all(v==0 for v in wf): continue
            c = 1
            for q,i in enumerate(tors): c = c*pow(z, ((D//d[i])*s[q]*w[i]) % D, p) % p
            pos = '*'.join(f'g_{free.index(i)}^{w[i]}' for i in free if w[i]>0) or '1'
            neg = '*'.join(f'g_{free.index(i)}^{-w[i]}' for i in free if w[i]<0) or '1'
            degs.append(f'({c}*{pos} - {neg})')
        hidx = m if K else 0
        vals = [random.randrange(2,p-1) for i in range(e)]
        subs = ', '.join(['x => R0_0', 'y => R0_1', 'z => R0_2'] + [f'g_{i} => {vals[i]}_kk' for i in range(e)])
        specvals = ', '.join(f'{vals[i]}_kk' for i in range(e))
        blk = f"""
<< "=== component {ci} s={list(s)} ===" << endl;
lns = {{{', '.join(L)}}};
Q = product lns;
Qs = sub(Q, {{{subs}}});
degsS = mdrOf Qs;
<< "random point: syz degrees " << degsS << " mdr=" << min degsS << (if #degsS == 2 then "  FREE" else "  nonfree") << endl;
"""
        if e > 0:
            blk += f"""
spec = map(kk, P, {{{specvals}}});
(N, mults) = zieglerMatrix(lns, {hidx}, spec);
<< "Ziegler multiplicities on H: " << mults << "  matrix " << numrows N << "x" << numcols N << endl;
I = minors(16, N);
degen = product {{{', '.join(degs) if degs else '1_P'}}};
J = saturate(I, degen);
<< "non-free locus inside the component (saturated by degenerations): " << (if J == ideal 1_P then "EMPTY" else toString mingens J) << endl;
if J != ideal 1_P then << "   dim " << dim J << " degree " << degree J << " decomposition " << apply(decompose J, q -> toString mingens q) << endl;
<< "   (for reference: gcd of the 16-minors before saturation has degree " << (if I == 0 then "INF" else first degree gcd flatten entries gens I) << ")" << endl;
"""
        blocks.append(blk)
    return '\n'.join(out) + '\n'.join(blocks) + '\n'

def main():
    entries = json.load(open(sys.argv[1]))
    lattices = {}
    for en in entries:
        key = (en['pattern'], tuple(tuple(t) for t in en['T']))
        lattices.setdefault(key, en)
    only = set(int(a) for a in sys.argv[2:]) if len(sys.argv) > 2 else None
    for (pname, T), en in lattices.items():
        if only is not None and en['sigid'] not in only: continue
        pat = [q for q in PATTERNS if q[0]==pname][0]; _,nA,nB,nC,K,t = pat
        comp = normalized_components(set(T), nA,nB,nC)
        tag = f"{pname}#{en['sigid']}"
        code = m2_script(set(T), nA,nB,nC,K, comp, tag)
        fn = f'jump_{pname}_{en["sigid"]}.m2'
        open(fn,'w').write(code)
        print(f'{tag}: dim={len(comp["free"])} torsion={[comp["d"][i] for i in comp["tors"]]} components(Galois classes, non-degenerate)={len(comp["comps"])} -> {fn}')

if __name__ == '__main__': main()
