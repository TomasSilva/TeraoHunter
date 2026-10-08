-- lattice P1#3: |A|,|B|,|C| = 5,5,5, K = 'xyz', t = 19; dim(realization component) = 1; torsion = []
kk = QQ;
R0 = kk[x,y,z];
P = kk[g_0];
R = P[x,y,z];

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

<< "=== component 0 s=[] ===" << endl;
lns = {(1)*x - (1*1)*y, (1)*x - (1*g_0^1)*y, (1)*x - (1*g_0^2)*y, (1)*x - (1*g_0^3)*y, (1)*x - (1*g_0^4)*y, (1)*z - (1*1)*x, (1)*z - (1*g_0^1)*x, (1)*z - (1*g_0^2)*x, (1)*z - (1*g_0^3)*x, (1)*z - (1*g_0^4)*x, (g_0^4)*y - (1*1)*z, (g_0^3)*y - (1*1)*z, (g_0^5)*y - (1*1)*z, (g_0^2)*y - (1*1)*z, (g_0^6)*y - (1*1)*z, x, y, z};
Q = product lns;
Qs = sub(Q, {x => R0_0, y => R0_1, z => R0_2, g_0 => 7516_kk});
degsS = mdrOf Qs;
<< "random point: syz degrees " << degsS << " mdr=" << min degsS << (if #degsS == 2 then "  FREE" else "  nonfree") << endl;

spec = map(kk, P, {7516_kk});
(N, mults) = zieglerMatrix(lns, 15, spec);
<< "Ziegler multiplicities on H: " << mults << "  matrix " << numrows N << "x" << numcols N << endl;
I = minors(16, N);
degen = product {(1*1 - g_0^1), (1*1 - g_0^2), (1*1 - g_0^3), (1*1 - g_0^4), (1*1 - g_0^1), (1*1 - g_0^2), (1*1 - g_0^3), (1*1 - g_0^1), (1*1 - g_0^2), (1*1 - g_0^1), (1*1 - g_0^1), (1*1 - g_0^2), (1*1 - g_0^3), (1*1 - g_0^4), (1*1 - g_0^1), (1*1 - g_0^2), (1*1 - g_0^3), (1*1 - g_0^1), (1*1 - g_0^2), (1*1 - g_0^1), (1*1 - g_0^1), (1*g_0^1 - 1), (1*1 - g_0^2), (1*g_0^2 - 1), (1*g_0^2 - 1), (1*1 - g_0^1), (1*g_0^3 - 1), (1*1 - g_0^3), (1*g_0^1 - 1), (1*g_0^4 - 1), (1*1 - g_0^4), (1*1 - g_0^3), (1*1 - g_0^5), (1*1 - g_0^2), (1*1 - g_0^6), (1*1 - g_0^3), (1*1 - g_0^2), (1*1 - g_0^4), (1*1 - g_0^1), (1*1 - g_0^5), (1*1 - g_0^2), (1*1 - g_0^1), (1*1 - g_0^3), (1*1 - g_0^4), (1*1 - g_0^1), (1*1 - g_0^2), (1*g_0^1 - 1), (1*1 - g_0^3), (1*g_0^1 - 1), (1*1 - g_0^1), (1*g_0^2 - 1), (1*1 - g_0^2), (1*1 - g_0^3), (1*1 - g_0^2), (1*1 - g_0^4), (1*1 - g_0^1), (1*1 - g_0^5), (1*1 - g_0^2), (1*1 - g_0^1), (1*1 - g_0^3), (1*1 - g_0^4), (1*1 - g_0^1), (1*1 - g_0^2), (1*g_0^1 - 1), (1*1 - g_0^3), (1*g_0^1 - 1), (1*1 - g_0^1), (1*g_0^2 - 1), (1*1 - g_0^2), (1*g_0^1 - 1), (1*g_0^2 - 1), (1*g_0^3 - 1), (1*1 - g_0^1), (1*1 - g_0^2), (1*1 - g_0^1), (1*1 - g_0^3), (1*1 - g_0^4), (1*1 - g_0^1), (1*1 - g_0^2), (1*g_0^1 - 1), (1*1 - g_0^3), (1*g_0^1 - 1), (1*1 - g_0^1), (1*g_0^2 - 1), (1*1 - g_0^2), (1*g_0^1 - 1), (1*g_0^2 - 1), (1*g_0^3 - 1), (1*1 - g_0^1), (1*g_0^2 - 1), (1*g_0^3 - 1), (1*g_0^1 - 1), (1*g_0^4 - 1), (1*1 - g_0^1), (1*1 - g_0^2), (1*g_0^1 - 1), (1*1 - g_0^3), (1*g_0^1 - 1), (1*1 - g_0^1), (1*g_0^2 - 1), (1*1 - g_0^2), (1*g_0^1 - 1), (1*g_0^2 - 1), (1*g_0^3 - 1), (1*1 - g_0^1), (1*g_0^2 - 1), (1*g_0^3 - 1), (1*g_0^1 - 1), (1*g_0^4 - 1), (1*g_0^3 - 1), (1*g_0^4 - 1), (1*g_0^2 - 1), (1*g_0^5 - 1), (1*g_0^1 - 1), (1*g_0^1 - 1), (1*1 - g_0^1), (1*g_0^2 - 1), (1*1 - g_0^2), (1*g_0^1 - 1), (1*g_0^2 - 1), (1*g_0^3 - 1), (1*1 - g_0^1), (1*g_0^2 - 1), (1*g_0^3 - 1), (1*g_0^1 - 1), (1*g_0^4 - 1), (1*g_0^3 - 1), (1*g_0^4 - 1), (1*g_0^2 - 1), (1*g_0^5 - 1), (1*g_0^1 - 1), (1*g_0^4 - 1), (1*g_0^5 - 1), (1*g_0^3 - 1), (1*g_0^6 - 1), (1*g_0^2 - 1)};
J = saturate(I, degen);
<< "non-free locus inside the component (saturated by degenerations): " << (if J == ideal 1_P then "EMPTY" else toString mingens J) << endl;
if J != ideal 1_P then << "   dim " << dim J << " degree " << degree J << " decomposition " << apply(decompose J, q -> toString mingens q) << endl;
<< "   (for reference: gcd of the 16-minors before saturation has degree " << (if I == 0 then "INF" else first degree gcd flatten entries gens I) << ")" << endl;

