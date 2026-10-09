-- s = 7, 8: random specialization of all but one modulus; a non-degenerate root of the univariate determinant proves
-- that the (7,9) locus is nonempty for that vector (sufficient for non-emptiness; emptiness would need the full computation).
kk = ZZ/32003; setRandomSeed 5;
buildM = (mu, e, P, S, pts) -> (
  nunk := 2*(e+1); rows := {}; u := S_0; v := S_1;
  for i from 0 to #mu-1 do ( p := pts#i; m := mu#i;
    for k from 0 to m-1 do ( row := for j from 0 to nunk-1 list (
        fj := if j <= e then u^j*v^(e-j) else 0_S; gj := if j <= e then 0_S else u^(j-e-1)*v^(2*e+1-j);
        if p === "inf" then lift(sub(diff(v^k, gj), {u => 1_S, v => 0_S}), P) else lift(sub(diff(u^k, fj - p*gj), {u => sub(p,S), v => 1_S}), P));
      rows = append(rows, row)));
  matrix rows);
partsOf = (n, k, mx) -> ( if k == 0 then (if n == 0 then {{}} else {}) else flatten for a from 1 to min(mx, n) list apply(partsOf(n - a, k - 1, a), p -> prepend(a, p)) );
P = kk[t]; S = P[u,v];
for s in {7, 8} do for mu in partsOf(16, s, 5) do (
  found := false; tries := 0;
  while not found and tries < 3 do ( tries = tries + 1;
    vals := apply(s-4, i -> random kk);
    pts := {"inf", 0_P, 1_P} | apply(vals, c -> sub(c, P)) | {t};
    M := buildM(mu, 7, P, S, pts); d := det M;
    if d == 0 then (found = true; << "s=" << s << " mu=" << mu << ": determinant identically zero in the last modulus (generic jump along a line)" << endl; break);
    degen := product({t, t-1} | apply(vals, c -> t - c));
    g := d; for i from 1 to 20 do (gg := gcd(g, degen); if gg == 1 then break; g = g // gg);
    if first degree g > 0 then (found = true; << "s=" << s << " mu=" << mu << ": (7,9) locus nonempty (non-degenerate factor of degree " << first degree g << " in the last modulus)" << endl << flush);
  );
  if not found then << "s=" << s << " mu=" << mu << ": no non-degenerate root found in 3 random specializations (EMPTY?)" << endl << flush;
);
