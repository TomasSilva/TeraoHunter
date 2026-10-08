#!/usr/bin/env python3
"""Realization scheme over F_p-bar of the lattice of an F_q-rational sub-arrangement of PG(2,q) (from pg2q_enum.py),
and mdr at a point of every non-degenerate component (M2).  Usage: pg2q_real.py q cands.json"""
import sys, json, itertools, random
sys.argv_ = sys.argv
q = int(sys.argv[1]); cands = json.load(open(sys.argv[2]))
exec(open('pg2q_enum.py').read().split('# points and lines')[0].split('q = int(sys.argv[1])')[1])  # field ops
def normalize(v):
    for i in range(3):
        if v[i] != 0:
            iv = inv(v[i]); return tuple(mul(x, iv) for x in v)
pts = sorted(set(normalize(v) for v in itertools.product(els, repeat=3) if any(v)))
def incident(P, L):
    s = 0
    for a,b in zip(P,L): s = add(s, mul(a,b))
    return s == 0
p = [r for r in range(2,q+1) if q % r == 0 and all(r % s for s in range(2,r))][0]
random.seed(5)
for ci, c in enumerate(cands):
    S = c['S']; n = len(S); L = [pts[j] for j in S]
    # multiple points: sets of line indices (within S) through each point of PG(2,q)
    mpts = []
    for P in pts:
        thr = [i for i in range(n) if incident(P, L[i])]
        if len(thr) >= 2: mpts.append(thr)
    multi = [t for t in mpts if len(t) >= 3]
    # frame: 4 lines, no three concurrent
    conc = set()
    for t in multi:
        for tr in itertools.combinations(t, 3): conc.add(tr)
    frame = None
    for F in itertools.combinations(range(n), 4):
        if all(tr not in conc for tr in itertools.combinations(F, 3)): frame = F; break
    others = [i for i in range(n) if i not in frame]
    names = {}
    for k, i in enumerate(others): names[i] = k
    varlist = ', '.join(f'a_{k}, b_{k}, c_{k}' for k in range(len(others)))
    def form(i):
        if i == frame[0]: return 'x'
        if i == frame[1]: return 'y'
        if i == frame[2]: return 'z'
        if i == frame[3]: return 'x+y+z'
        k = names[i]; return f'a_{k}*x + b_{k}*y + c_{k}*z'
    # incidence equations: for each multiple point, 3x3 minors of the coefficient matrix of its lines
    eqs = []
    for t in multi:
        eqs.append(f'minors(3, matrix apply({{{", ".join(form(i) for i in t)}}}, l -> {{coefficient(x,l), coefficient(y,l), coefficient(z,l)}}))')
    # forbidden: any triple not in a multiple point must have nonzero determinant; any pair distinct
    forb = []
    for tr in itertools.combinations(range(n), 3):
        if tr not in conc: forb.append(f'det matrix apply({{{", ".join(form(i) for i in tr)}}}, l -> {{coefficient(x,l), coefficient(y,l), coefficient(z,l)}})')
    # coordinates of every line in the frame where the 4 frame lines are x, y, z, x+y+z (over F_q):
    # write l4 = sum lam_i l_i, then coordinates of v are mu with v = sum mu_i (lam_i l_i)
    def solve3(M, v):
        # solve M^T mu = v over GF(q) by brute force (q^3 small)
        for mu in itertools.product(els, repeat=3):
            w = [0,0,0]
            for i in range(3):
                for j in range(3):
                    w[j] = add(w[j], mul(mu[i], M[i][j]))
            if tuple(w) == tuple(v): return mu
    l1,l2,l3,l4 = [L[i] for i in frame]
    lam = solve3([l1,l2,l3], l4)
    basis = [tuple(mul(lam[i], l[j]) for j in range(3)) for i,l in enumerate([l1,l2,l3])]
    coords = {i: solve3(basis, L[i]) for i in others}
    charts = []
    for i in others:
        k = names[i]; mu = coords[i]
        pos = [j for j in range(3) if mu[j] != 0][0]
        charts.append(['a','b','c'][pos] + f'_{k} - 1')
    code = f'''
-- lattice {ci}: n={n}, class ({c['d1']},{n-1-c['d1']}), profile {c['profile']}; frame lines {list(frame)}
kk = ZZ/{p};
R = kk[{varlist}][x,y,z];
P = coefficientRing R;
I = ideal(0_P);
scan({{{", ".join(eqs)}}}, m -> I = I + sub(ideal m, P));
I = I + ideal({{{", ".join(charts)}}});
I = trim I;
<< "lattice {ci}: n={n} class=({c['d1']},{n-1-c['d1']}) frame={list(frame)}  dim(realization scheme in charts)=" << dim I << endl;
forb = apply({{{", ".join(forb)}}}, f -> sub(f, P));
comps = decompose I;
<< "   components over F_{p}: " << #comps << endl;
mdrOf = Q -> (RR := ring Q; jac := matrix{{{{diff(RR_0,Q),diff(RR_1,Q),diff(RR_2,Q)}}}}; SS := mingens image syz jac; sort apply(numcols SS, i -> first (degrees source SS)#i - (first degree Q - 1)));
scan(comps, Pc -> (
   degen := any(forb, f -> f % Pc == 0);
   << "   component: dim " << dim Pc << " degree " << degree Pc << (if degen then "  DEGENERATE (forced extra incidence)" else "  non-degenerate") << endl;
   if not degen and dim Pc == 0 then (
      -- residue field = GF({p}^k); build the arrangement over it and compute mdr
      Lf := toField(P/Pc);
      S := Lf[x,y,z];
      lns := apply({{{", ".join(form(i) for i in range(n))}}}, l -> sub(l, S));
      degs := mdrOf product lns;
      << "      residue field degree " << degree Pc << ": syz degrees " << degs << " mdr=" << min degs << (if #degs == 2 then "  FREE <== COUNTEREXAMPLE?" else "  nonfree") << endl;
   );
   if not degen and dim Pc > 0 then << "      POSITIVE-DIMENSIONAL non-degenerate component: " << toString mingens Pc << endl;
));
'''
    open(f'real_pg{q}_{ci}.m2','w').write(code)
    print(f'wrote real_pg{q}_{ci}.m2 (n={n}, {len(others)} free lines, {len(multi)} multiple points, {len(forb)} forbidden triples)')
