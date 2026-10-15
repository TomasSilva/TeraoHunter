-- real points on the 5-point jump curves: for a on a fine grid, count real roots b of F(a,b) (non-degenerate)
load "jump5.m2";
T = QQ[bb];
for k in {{5,3,3,3,2},{4,4,4,2,2},{4,4,3,3,2},{4,3,3,3,3}} do (
   F := jumpCurve k;
   G := product apply(select(toList factor F, p -> first degree(p#0) > 1), p -> p#0);
   cnt := 0; samples := {};
   for i from -60 to 60 do (
      v := i/10 + 1/37;
      g := sub(G, {a => v, b => bb});
      -- Sturm-free count: sign changes on a fine b-grid in [-20, 20] plus large-|b| check
      vals := apply(toList(-4000..4000), j -> sub(g, {bb => j/200 + 1/101}));
      ch := #select(toList(0..#vals-2), j -> (vals#j) * (vals#(j+1)) < 0);
      if ch > 0 then (cnt = cnt + 1; if #samples < 3 then samples = append(samples, (v, ch))));
   print("k " | toString k | ": a-values (of 121) with real roots b in [-20,20]: " | toString cnt | " e.g. " | toString samples));
