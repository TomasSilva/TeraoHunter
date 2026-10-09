-- Faster version over ZZ/32003: e = 7 by the determinant, e = 6 by two random 14x16 combinations (sufficient for emptiness;
-- if both vanish somewhere non-degenerate we report "possible" and confirm with the true minors).
kk = ZZ/32003;
setRandomSeed 11;
buildM = (mu, e, P, S, pts) -> (
  nunk := 2*(e+1); rows := {};
  u := S_0; v := S_1;
  for i from 0 to #mu-1 do (
    p := pts#i; m := mu#i;
    for k from 0 to m-1 do (
      row := for j from 0 to nunk-1 list (
        fj := if j <= e then u^j*v^(e-j) else 0_S;
        gj := if j <= e then 0_S else u^(j-e-1)*v^(2*e+1-j);
        if p === "inf" then lift(sub(diff(v^k, gj), {u => 1_S, v => 0_S}), P)
        else lift(sub(diff(u^k, fj - p*gj), {u => sub(p,S), v => 1_S}), P));
      rows = append(rows, row)));
  matrix rows);
jumpTest = (mu) -> (
  s := #mu; nv := s - 3;
  P := if nv > 0 then kk[t_1..t_nv] else kk;
  S := P[u,v];
  pts := {"inf", 0_P, 1_P} | (if nv > 0 then apply(nv, i -> P_i) else {});
  disc := 1_P;
  if nv > 0 then (scan(nv, i -> disc = disc * P_i * (P_i - 1)); for i from 0 to nv-1 do for j from i+1 to nv-1 do disc = disc * (P_i - P_j));
  M7 := buildM(mu, 7, P, S, pts);
  d7 := det M7;
  r7 := if d7 == 0 then "ZERO(generic jump)" else if nv == 0 then "EMPTY(constant)" else (
     J := saturate(ideal d7, disc); if J == ideal 1_P then "EMPTY" else ("nonempty(dim " | toString dim J | ")"));
  M6 := buildM(mu, 6, P, S, pts);
  r6 := if nv == 0 then (if rank M6 < 14 then "nonempty" else "EMPTY") else (
     L1 := random(P^14, P^16); L2 := random(P^14, P^16);
     I6 := ideal(det(L1*M6), det(L2*M6));
     J6 := saturate(I6, disc);
     if J6 == ideal 1_P then "EMPTY" else (
        -- confirm with the true minors on the candidate locus
        J6t := saturate(minors(14, M6) + J6, disc);
        if J6t == ideal 1_P then "EMPTY(after minors)" else ("nonempty(dim " | toString dim J6t | ")")));
  (r7, r6));
partsOf = (n, k, mx) -> ( if k == 0 then (if n == 0 then {{}} else {}) else flatten for a from 1 to min(mx, n) list apply(partsOf(n - a, k - 1, a), p -> prepend(a, p)) );
for s from 4 to 8 do for mu in partsOf(16, s, 5) do (
  r := jumpTest(mu);
  << "s=" << s << " mu=" << mu << "  e1<=7: " << r#0 << "  |  e1<=6: " << r#1 << endl << flush;
);
