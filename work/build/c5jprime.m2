S = QQ[x,y,z];
pen = (dp, dq, offs) -> apply(offs, v -> dq*x - dp*y - v*z);
arrs = {
 ("R17a", {z} | pen(1,0,{0,1,2,3}) | pen(0,1,{0,2,4}) | pen(1,1,{2,3,5}) | pen(4/3,1,{4}) | pen(2,1,{2,4,5,6,8})),
 ("R17b", {z} | pen(1,0,{0,1,2,3,4}) | pen(0,1,{0,1,2,4}) | pen(1,1,{4}) | pen(-1/2,1,{1}) | pen(-1,1,{-2,-1,0,1,2}))};
for A in arrs do (
  f := product(A#1);
  Sy := syz matrix{{diff(x,f), diff(y,f), diff(z,f)}};
  cols := entries transpose Sy; dg := apply(cols, c -> first degree c#0);
  th := cols#(position(dg, d -> d == min dg));
  Xf := h -> th#0*diff(x,h) + th#1*diff(y,h) + th#2*diff(z,h);
  J := det matrix{{x,y,z}, th, apply(th, Xf)};
  Jp := J // f;
  print("== " | A#0 | ": deg J = " | toString(first degree J) | ", Q divides J: " | toString(J % f == 0) | ", deg J' = " | toString(first degree Jp));
  for fac in toList factor Jp do if first degree(fac#0) > 0 then print("   factor of J': " | toString(fac#0) | "  ^" | toString(fac#1) | "  (degree " | toString(first degree fac#0) | ")");
  IZ := saturate ideal th; print("   Bourbaki point: " | toString IZ);
  -- is the conic part smooth? singular points of J'
  sing := saturate(ideal(Jp, diff(x,Jp), diff(y,Jp), diff(z,Jp)));
  print("   Sing(J') degree " | toString degree sing | ", contains P0: " | toString(isSubset(sing, IZ)));
);
