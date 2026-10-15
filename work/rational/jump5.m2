-- Jump locus of a 5-point multiarrangement: points x, y, x - y, x - a y, x - b y, multiplicities k (K even):
-- F(a,b) = det of the degree-(K/2 - 1) condition matrix; its zero set (off a,b in {0,1}, a = b) is where e1 < K/2.
R = QQ[a, b]; S = R[x, y];
condMatrix5 = (k, e) -> (
   fb := apply(e + 1, i -> x^(e - i) * y^i);
   cols := 2 * (e + 1);
   rows := {};
   for i from 0 to e do if i > e - k#0 then rows = append(rows, apply(cols, c -> if c == i then 1_R else 0_R));
   for i from 0 to e do if i < k#1 then rows = append(rows, apply(cols, c -> if c == e + 1 + i then 1_R else 0_R));
   pts := {(1_R, 1_S), (a, a_S), (b, b_S)};       -- (x - t y): t = 1, a, b ; condition on f - t g at x = t
   for q from 0 to 2 do (
      t := (pts#q)#0; tS := (pts#q)#1;
      for j from 0 to k#(q + 2) - 1 do rows = append(rows, apply(cols, c -> (
         mon := if c <= e then fb#c else -(tS * fb#(c - e - 1));
         lift(sub(diff(x^j, sub(mon, {y => 1_S})), {x => tS}), R)))));
   matrix rows);
jumpCurve = k -> (K := sum k; e := K // 2 - 1; M := condMatrix5(k, e);
   if numrows M != numcols M then error "not square";
   det M);
