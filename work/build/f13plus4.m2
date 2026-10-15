load "famtools.m2";
-- ACKN 13-line family over QQ (x rescaled), lambda = l1/l0
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
print("F13: points " | toString tally w | ", degenerations " | toString keys degF);
-- candidate lines of generic weight 5 through two unjoined multiple points
cands = {};
seen = new MutableHashTable;
for pr in subsets(#pts, 2) do (
  i := pr#0; j := pr#1;
  if #(set(pts#i) * set(pts#j)) > 0 then continue;
  h := crossv(P#i, P#j);
  if seen#?h then continue; seen#h = true;
  on := select(toList(0..#pts-1), k -> dot3(h, P#k) == 0);
  W := sum(on, k -> w#k);
  if W == 5 then cands = append(cands, (h, on)));
print("candidate lines with |F13 meet l| = 8: " | toString(#cands));
for c in cands do print("   " | toString(c#0) | "  through points of weights " | toString apply(c#1, k -> w#k));
-- compatibility: intersection of two candidates lies on a line of F13 identically
nc = #cands;
compat = new MutableHashTable;
for pr in subsets(nc, 2) do (
  X := crossv((cands#(pr#0))#0, (cands#(pr#1))#0);
  if X == {0_R,0_R,0_R} then continue;
  compat#pr = any(F13, m -> dot3(m, X) == 0));
edges = select(subsets(nc,2), pr -> compat#?pr and compat#pr);
print("compatible pairs: " | toString(#edges) | " of " | toString(binomial(nc,2)));
cliques4 = select(subsets(nc, 4), S -> all(subsets(S,2), pr -> compat#?pr and compat#pr));
print("4-cliques: " | toString(#cliques4));
for S in cliques4 do (
  A := F13 | apply(S, i -> (cands#i)#0);
  print("== clique " | toString S);
  (ptsA, degA) = latticeOf A;
  ok := windowReport(A, ptsA);
  if not ok then continue;
  Th := thetaH(A, ptsA, 0);
  print("   Theta_H1 degree " | toString(first degree Th) | (if Th == 0 then "  ZERO: generically NON-FREE" else ""));
  if Th != 0 then (
    inter := classifyTheta(A, ptsA, degA, 0, Th);
    print("   interior factors: " | toString inter));
  hid := hiddenLines(A, ptsA, degA, 9);
  print("   hidden lines of weight >= 9 (interior): " | toString(#hid));
  for x in hid do print("      " | toString x));
