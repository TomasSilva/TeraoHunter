-- Count the (6,10) jump positions of a 6-point line (Ziegler multiplicities k, K = 16) by solving the
-- kernel system directly, and compare with the planar-graph prediction of TC_GLOBAL.md (critically fixed
-- maps of degree 6 with five-edge charge graphs). Over ZZ/32003; points at inf, 0, 1, t1, t2, t3.
kk = ZZ/32003; setRandomSeed 7;
countJumps = k -> (
  e := 6; nunk := 2*(e+1);
  R := kk[t_1..t_3, v_1..v_nunk];
  pts := {"inf", 0_R, 1_R, t_1, t_2, t_3};
  S := R[u, w];
  f := sum(e+1, j -> v_(j+1) * u^j * w^(e-j)); g := sum(e+1, j -> v_(e+2+j) * u^j * w^(e-j));
  eqs := {};
  for i from 0 to 5 do (
    p := pts#i; m := k#i;
    for a from 0 to m-1 do (
      if p === "inf" then eqs = append(eqs, sub(diff(w^a, g), {u => 1_S, w => 0_S}))
      else eqs = append(eqs, sub(diff(u^a, f - sub(p,S)*g), {u => sub(p,S), w => 1_S}))));
  eqs = apply(eqs, q -> lift(q, R));
  norm := sum(nunk, j -> random(kk) * v_(j+1)) - 1;
  I := ideal(eqs | {norm});
  disc := t_1*t_2*t_3*(t_1-1)*(t_2-1)*(t_3-1)*(t_1-t_2)*(t_1-t_3)*(t_2-t_3);
  J := saturate(I, disc);
  (dim J, degree J));
for k in {{3,3,3,3,2,2}, {2,2,3,3,3,3}} do (
  print("k = " | toString k | "  (dim, degree) = " | toString countJumps k));
