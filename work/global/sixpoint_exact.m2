-- Exact version: the (6,10) positions of each 6-point vector form a zero-dimensional scheme whose
-- rational univariate representation (msolve) has eliminant P(T) of degree = number of positions.
-- Rational positions <-> linear factors of P over QQ. The factor degrees are the Galois orbit sizes
-- of the corresponding critically fixed maps of degree 6.
needsPackage "Msolve";
vecs = apply(select(partitions 16, p -> #p == 6 and max toList p <= 5), p -> reverse sort toList p);
e = 6; nunk = 2*(e+1);
R = QQ[v_1..v_nunk, t_1..t_3, zz];
S = R[u, w];
for k in vecs do (
  pts := {"inf", 0_R, 1_R, t_1, t_2, t_3};
  f := sum(e+1, j -> v_(j+1) * u^j * w^(e-j)); g := sum(e+1, j -> v_(e+2+j) * u^j * w^(e-j));
  eqs := {};
  for i from 0 to 5 do (
    p := pts#i; m := k#i;
    for a from 0 to m-1 do (
      if p === "inf" then eqs = append(eqs, sub(diff(w^a, g), {u => 1_S, w => 0_S}))
      else eqs = append(eqs, sub(diff(u^a, f - sub(p,S)*g), {u => sub(p,S), w => 1_S}))));
  eqs = apply(eqs, q -> lift(q, R));
  setRandomSeed 7;
  nrm := sum(nunk, j -> (random(-5,5) + 0_QQ) * v_(j+1)) - 1;
  disc := t_1*t_2*t_3*(t_1-1)*(t_2-1)*(t_3-1)*(t_1-t_2)*(t_1-t_3)*(t_2-t_3);
  I := ideal(eqs | {nrm, zz*disc - 1});
  rur := msolveRUR I;
  P := rur#"findRootsUniPoly";
  fl := toList factor P;
  degs := sort flatten apply(select(fl, q -> first degree(q#0) > 0), q -> toList(q#1 : first degree(q#0)));
  nlin := #select(degs, d -> d == 1);
  print("k = " | toString k | ": positions " | toString(rur#"degree") | ", rational " | toString nlin |
        ", Galois orbit sizes " | toString degs));
