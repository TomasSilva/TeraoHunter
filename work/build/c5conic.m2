S = QQ[x,y,z];
pen = (dp, dq, offs) -> apply(offs, v -> dq*x - dp*y - v*z);
data = {
 ("R17a", {z} | pen(1,0,{0,1,2,3}) | pen(0,1,{0,2,4}) | pen(1,1,{2,3,5}) | pen(4/3,1,{4}) | pen(2,1,{2,4,5,6,8}), 2*x^2-2*x*y+2*y^2-11*x*z+10*y*z+20*z^2),
 ("R17b", {z} | pen(1,0,{0,1,2,3,4}) | pen(0,1,{0,1,2,4}) | pen(1,1,{4}) | pen(-1/2,1,{1}) | pen(-1,1,{-2,-1,0,1,2}), x^2+x*y+y^2-2*x*z+2*y*z+2*z^2)};
vec = g -> {coefficient(x,g), coefficient(y,g), coefficient(z,g)};
cr = (p,q) -> {p#1*q#2-p#2*q#1, p#2*q#0-p#0*q#2, p#0*q#1-p#1*q#0};
for D in data do (
  A := D#1; Cq := D#2; n := #A; vs := apply(A, vec);
  pts := unique apply(subsets(n,2), p -> sort unique({p#0,p#1} | select(toList(0..n-1), k -> not member(k,p) and det matrix{vs#(p#0),vs#(p#1),vs#k} == 0)));
  onC := select(pts, P -> (c := cr(vs#(P#0), vs#(P#1)); sub(Cq, {x => c#0, y => c#1, z => c#2}) == 0));
  print("== " | D#0 | ": intersection points on the conic: " | toString(#onC) | ", multiplicities " | toString sort apply(onC, P -> #P));
  -- lines of A tangent to the conic
  tang := select(toList(0..n-1), i -> (h := A#i; I := ideal(Cq, h); degree saturate(I) == 1 or (degree radical saturate I) == 1));
  print("   lines of A tangent to the conic: " | toString(#tang));
);
