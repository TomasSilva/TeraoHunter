-- Jump loci of 4-point multiarrangements on P^1 over QQ.
-- Points x, y, x - y, x - mu*y with Ziegler multiplicities k = (k1,k2,k3,k4), K = sum.
-- D_e(mu) = { f d_x + g d_y of degree e : x^k1 | f, y^k2 | g, (x-y)^k3 | f - g, (x-mu y)^k4 | f - mu g }.
-- e1(mu) = least e with D_e(mu) != 0.  {mu : e1(mu) <= e} = zeros of the gcd of maximal minors of the
-- degree-e condition matrix (when it has at least as many rows as columns).
R = QQ[mu]; S = R[x, y];
condMatrix = (k, e) -> (
   fb := apply(e + 1, i -> x^(e - i) * y^i);
   cols := 2 * (e + 1);
   rows := {};
   for i from 0 to e do if i > e - k#0 then rows = append(rows, apply(cols, c -> if c == i then 1_R else 0_R));
   for i from 0 to e do if i < k#1 then rows = append(rows, apply(cols, c -> if c == e + 1 + i then 1_R else 0_R));
   for j from 0 to k#2 - 1 do rows = append(rows, apply(cols, c -> (
        mon := if c <= e then fb#c else -fb#(c - e - 1);
        lift(sub(diff(x^j, sub(mon, {y => 1_S})), {x => 1_S}), R))));
   for j from 0 to k#3 - 1 do rows = append(rows, apply(cols, c -> (
        mon := if c <= e then fb#c else -(mu * fb#(c - e - 1));
        lift(sub(diff(x^j, sub(mon, {y => 1_S})), {x => mu_S}), R))));
   matrix rows);
genericE1 = k -> (e := 0; while true do (
      M := condMatrix(k, e); Mr := sub(M, {mu => 7/3});
      if rank Mr < numcols Mr then return e; e = e + 1));
jumpPoly = (k, e) -> (
   M := condMatrix(k, e); c := numcols M;
   if numrows M < c then return 0_R;
   gcd flatten entries gens minors(c, M));
ratRoots = g -> (
   fl := select(toList factor g, p -> first degree(p#0) == 1);
   unique apply(fl, p -> (h := p#0; - (coefficient(1_R, h)) / (coefficient(mu, h)))));
