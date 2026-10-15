-- Real and rational (6,10) positions of 6-point lines (critically fixed maps of degree 6), via msolve.
-- A real solution is tested for rationality by rational reconstruction of the isolating-box midpoints
-- (denominators <= 10^4) followed by an exact check that the point lies on the variety.
needsPackage "Msolve";
vecs = apply(select(partitions 16, p -> #p == 6 and max toList p <= 5), p -> reverse sort toList p);
e = 6; nunk = 2*(e+1);
R = QQ[v_1..v_nunk, t_1..t_3, zz];
S = R[u, w];
ratApprox = x -> (   -- best rational approximation with denominator <= 10^4 (continued fractions)
  a := floor x; p0 := 1; q0 := 0; p1 := a; q1 := 1; y := x - a; cnt := 0;
  while y != 0 and cnt < 40 do (
    y = 1/y; a = floor y; y = y - a;
    p2 := a*p1 + p0; q2 := a*q1 + q0;
    if q2 > 10000 then break;
    (p0, q0, p1, q1) = (p1, q1, p2, q2); cnt = cnt + 1);
  p1/q1);
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
  sols := msolveRealSolutions(I, QQi);
  rats := {};
  for sl in sols do (
    cand := apply({nunk, nunk+1, nunk+2}, i -> ratApprox((sl#i#0 + sl#i#1)/2));
    Jc := I + ideal(t_1 - cand#0, t_2 - cand#1, t_3 - cand#2);
    if Jc != ideal(1_R) then rats = append(rats, cand));
  print("k = " | toString k | ": real positions " | toString(#sols) | ", rational positions " | toString(#rats) |
        (if #rats > 0 then ": (t1,t2,t3) = " | toString unique rats else "")));
