-- 8-point lines (and 7-point) : is e1 <= 7 possible at all?  determinant of the 16x16 degree-7 matrix, saturated by the discriminant.
kk = ZZ/32003;
buildM = (mu, e, P, S, pts) -> (
  nunk := 2*(e+1); rows := {}; u := S_0; v := S_1;
  for i from 0 to #mu-1 do ( p := pts#i; m := mu#i;
    for k from 0 to m-1 do ( row := for j from 0 to nunk-1 list (
        fj := if j <= e then u^j*v^(e-j) else 0_S; gj := if j <= e then 0_S else u^(j-e-1)*v^(2*e+1-j);
        if p === "inf" then lift(sub(diff(v^k, gj), {u => 1_S, v => 0_S}), P) else lift(sub(diff(u^k, fj - p*gj), {u => sub(p,S), v => 1_S}), P));
      rows = append(rows, row)));
  matrix rows);
partsOf = (n, k, mx) -> ( if k == 0 then (if n == 0 then {{}} else {}) else flatten for a from 1 to min(mx, n) list apply(partsOf(n - a, k - 1, a), p -> prepend(a, p)) );
for s in {8, 7} do for mu in partsOf(16, s, 5) do (
  nv := s - 3; P := kk[t_1..t_nv]; S := P[u,v];
  pts := {"inf", 0_P, 1_P} | apply(nv, i -> P_i);
  M := buildM(mu, 7, P, S, pts); d := det M;
  -- non-degenerate factor exists iff d is not a unit times a product of t_i, t_i - 1, t_i - t_j
  disc := 1_P; scan(nv, i -> disc = disc * P_i * (P_i - 1)); for i from 0 to nv-1 do for j from i+1 to nv-1 do disc = disc * (P_i - P_j);
  r := if d == 0 then "ZERO(generic jump)" else (J := saturate(ideal d, disc); if J == ideal 1_P then "EMPTY" else "nonempty(hypersurface, degree " | toString(first degree d) | ")");
  << "s=" << s << " mu=" << mu << "  e1<=7: " << r << endl << flush;
);
