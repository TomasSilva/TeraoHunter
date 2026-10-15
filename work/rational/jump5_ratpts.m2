-- For the eleven 5-point H-vectors of the (8,8) window (partitions of 16 into five parts <= 5):
-- the jump curve F_k(a,b) = 0 and its rational points with a of height <= Hmax (b unrestricted).
load "jump5.m2";
Hmax = value(scriptCommandLine#1);
avals = unique flatten for q from 1 to Hmax list for p from -Hmax to Hmax list (if gcd(p, q) == 1 then p/q else continue);
avals = select(avals, v -> v != 0 and v != 1);
T = QQ[bb];
vecs = {{5,5,4,1,1},{5,5,3,2,1},{5,5,2,2,2},{5,4,4,2,1},{5,4,3,3,1},{5,4,3,2,2},{5,3,3,3,2},{4,4,4,3,1},{4,4,4,2,2},{4,4,3,3,2},{4,3,3,3,3}};
for k in vecs do (
   F := jumpCurve k;
   core := product select(toList factor F, p -> first degree(p#0) > 1 or (degree(a, p#0) == 1 and degree(b, p#0) == 1 and p#0 != a - b));
   G := product apply(select(toList factor F, p -> first degree(p#0) > 1), p -> p#0);
   pts := {};
   for v in avals do (
      g := sub(G, {a => v, b => bb});
      if g == 0 then (pts = append(pts, (v, "ALL")); continue);
      for p in toList factor g do if first degree(p#0) == 1 then (
         h := p#0; r := -(coefficient(1_T, h)) / coefficient(bb, h);
         if r != 0 and r != 1 and r != v then pts = append(pts, (v, r))));
   print("k " | toString k | ": non-degenerate part degree " | toString(first degree G) | ", rational points with height(a) <= " | toString Hmax | ": " | toString(#pts));
   print("   " | toString take(sort pts, 40));
);
