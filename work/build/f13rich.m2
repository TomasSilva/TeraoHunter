load "famtools.m2";
F13 = apply({
 {-l0, -l0, l1+l0}, {0_R, 2*l0, l1+l0}, {l0, -l0, l1+l0},
 {l0, -l0, l1-2*l0}, {-l0, -l0, l1-2*l0}, {0_R, 2*l0, l1-2*l0},
 {0_R, 2*l0, -2*l1+l0}, {l0, -l0, -2*l1+l0}, {-l0, -l0, -2*l1+l0},
 {(l0-l1)*l0, (l1+l0)*l0, -l1^2+l1*l0-l0^2},
 {l1*l0, (l1-2*l0)*l0, -l1^2+l1*l0-l0^2},
 {-l0^2, (l0-2*l1)*l0, -l1^2+l1*l0-l0^2},
 {0_R, 0_R, 1_R}}, prim);
(pts, degF) = latticeOf F13;
P = apply(pts, S -> crossv(F13#(S#0), F13#(S#1)));
w = apply(pts, S -> #S - 1);
tal = new MutableHashTable; seen = new MutableHashTable;
for pr in subsets(#pts, 2) do (
  i := pr#0; j := pr#1;
  if #(set(pts#i) * set(pts#j)) > 0 then continue;
  h := crossv(P#i, P#j);
  if seen#?h then continue; seen#h = true;
  W := sum(select(toList(0..#pts-1), k -> dot3(h, P#k) == 0), k -> w#k);
  tal#W = (if tal#?W then tal#W else 0) + 1);
print("F13 external lines through two unjoined multiple points, by generic weight: " | toString sort pairs new HashTable from tal);
