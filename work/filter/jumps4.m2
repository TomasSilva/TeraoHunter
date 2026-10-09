-- explicit jumping cross-ratios (s = 4) and jump curves (s = 5) for e1 <= 7, over QQ
buildM = (mu, e, P, S, pts) -> (
  nunk := 2*(e+1); rows := {}; u := S_0; v := S_1;
  for i from 0 to #mu-1 do ( p := pts#i; m := mu#i;
    for k from 0 to m-1 do ( row := for j from 0 to nunk-1 list (
        fj := if j <= e then u^j*v^(e-j) else 0_S; gj := if j <= e then 0_S else u^(j-e-1)*v^(2*e+1-j);
        if p === "inf" then lift(sub(diff(v^k, gj), {u => 1_S, v => 0_S}), P) else lift(sub(diff(u^k, fj - p*gj), {u => sub(p,S), v => 1_S}), P));
      rows = append(rows, row)));
  matrix rows);
partsOf = (n, k, mx) -> ( if k == 0 then (if n == 0 then {{}} else {}) else flatten for a from 1 to min(mx, n) list apply(partsOf(n - a, k - 1, a), p -> prepend(a, p)) );
P = QQ[t]; S = P[u,v];
for mu in partsOf(16, 4, 5) do (
  M := buildM(mu, 7, P, S, {"inf", 0_P, 1_P, t}); d := det M;
  << "s=4 mu=" << mu << ": det = " << factor d << endl;
);
P = QQ[t1,t2]; S = P[u,v];
for mu in partsOf(16, 5, 5) do (
  M := buildM(mu, 7, P, S, {"inf", 0_P, 1_P, t1, t2}); d := det M;
  f := factor d;
  << "s=5 mu=" << mu << ": det factors (degrees): " << apply(toList f, g -> (first degree g#0, g#1)) << endl;
  nondeg := select(toList f, g -> (h := g#0; not (h == t1 or h == t2 or h == t1-1 or h == t2-1 or h == t1-t2 or h == t2 - t1 or first degree h == 0)));
  << "     non-degenerate factors: " << apply(nondeg, g -> g#0) << endl;
);
