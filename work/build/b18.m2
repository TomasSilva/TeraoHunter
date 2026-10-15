load "famtools.m2";
-- The 15-line family F15(lambda) (free (7,7)): lines joining Q0..Q5.
Q = {{1_R,0_R,0_R}, {0_R,1_R,0_R}, {0_R,0_R,1_R}, {1_R,1_R,1_R},
     {l1^2, l1*l0, 2*l1*l0-l0^2}, {l1*(2*l1-l0), l1^2, (2*l1-l0)*l0}};
F = apply(subsets(6,2), p -> crossv(Q#(p#0), Q#(p#1)));
(pts, degF) = latticeOf F;
P = apply(pts, S -> crossv(F#(S#0), F#(S#1)));
w = apply(pts, S -> #S - 1);
print("F15: points " | toString tally w | ", degenerations " | toString keys degF);
-- candidate lines: through two unjoined multiple points, generic weight 7 or 8 (8 or 7 points on F)
cands = {}; seen = new MutableHashTable;
for pr in subsets(#pts, 2) do (
  i := pr#0; j := pr#1;
  if #(set(pts#i) * set(pts#j)) > 0 then continue;
  h := crossv(P#i, P#j);
  if seen#?h then continue; seen#h = true;
  on := select(toList(0..#pts-1), k -> dot3(h, P#k) == 0);
  W := sum(on, k -> w#k);
  if W == 7 or W == 8 then cands = append(cands, (h, on, 15 - W)));
nc = #cands;
print("candidate lines: " | toString nc | "  (points on F: " | toString tally apply(cands, c -> c#2) | ")");
-- pairwise: intersection on F?
onF = new MutableHashTable;
for pr in subsets(nc,2) do (
  X := crossv((cands#(pr#0))#0, (cands#(pr#1))#0);
  onF#pr = (X == {0_R,0_R,0_R}) or any(F, m -> dot3(m, X) == 0));
-- general Yoshinaga rank test on line a for exponents (d1, d2): free iff no derivation of degree d1-1 on the Ziegler restriction
zieglerRank = (H, ptsH, a, e) -> (
  n := #H;
  onH := select(ptsH, S -> member(a, S));
  k := position({0,1,2}, j -> (H#a)#j != 0);
  ij := delete(k, {0,1,2});
  rows := flatten apply(onH, S -> (
      b := first select(S, x -> x != a);
      Pt := crossv(H#a, H#b);
      Pi := Pt#(ij#0); Pj := Pt#(ij#1);
      m := #S - 1;
      apply(m, aa -> (
         Dd := u -> if u >= aa and (e-u) >= (m-1-aa) then ff(u,aa)*ff(e-u,m-1-aa)*Pi^(u-aa)*Pj^(e-u-(m-1-aa)) else 0_R;
         apply(e+1, u -> Pj*Dd(u)) | apply(e+1, u -> -Pi*Dd(u))))));
  M := matrix rows;
  (rank sub(M, frac R), 2*(e+1)));
-- triples
good = {};
for T in subsets(nc, 3) do (
  s := apply(T, i -> (cands#i)#2);
  newpts := 0;
  prs := subsets(T,2);
  X := apply(prs, pr -> crossv((cands#(pr#0))#0, (cands#(pr#1))#0));
  offF := select(toList(0..2), i -> not onF#(sort prs#i));
  -- distinct new points among those off F
  newpts = #unique apply(offF, i -> X#i);
  if sum s + newpts != 24 then continue;
  B := F | apply(T, i -> (cands#i)#0);
  (ptsB, degB) = latticeOf B;
  b2 := sum(ptsB, S -> #S - 1);
  if b2 != 87 then continue;
  mB := max apply(ptsB, S -> #S);
  sB := apply(18, a -> #select(ptsB, S -> member(a,S)));
  if mB > 7 or max sB > 9 then continue;
  (rk, need) := zieglerRank(B, ptsB, 0, 6);
  fr := (rk == need);
  if not fr then continue;
  good = append(good, (T, ptsB, sB));
  print("B free (7,10): lines " | toString T | " points on F " | toString s | ", mults " | toString tally apply(ptsB, S -> #S) | ", s_H " | toString tally sB));
print("free B found: " | toString(#good));
-- A2 = B minus a 7-point line
surv = 0;
for g in good do (
  (T, ptsB, sB) = g;
  B := F | apply(T, i -> (cands#i)#0);
  for a from 0 to 17 do (
    if sB#a != 7 then continue;
    A := drop(B, {a,a});
    (ptsA, degA) = latticeOf A;
    wA := apply(ptsA, S -> #S - 1);
    PA := apply(ptsA, S -> crossv(A#(S#0), A#(S#1)));
    mA := max apply(ptsA, S -> #S);
    sA := apply(17, c -> #select(ptsA, S -> member(c,S)));
    reasons := {};
    if mA > 6 then reasons = append(reasons, "mult>6");
    if min sA < 4 or max sA > 8 then reasons = append(reasons, "s_H");
    if not (member(4,sA) or member(5,sA) or member(8,sA)) then reasons = append(reasons, "no 4/5/8 line");
    -- external lines: weights of lines through unjoined pairs (generic)
    bad8 := 0; rich := 0; sig := 0;
    seenL := new MutableHashTable;
    for pr in subsets(#ptsA,2) do (
      i := pr#0; j := pr#1;
      if #(set(ptsA#i) * set(ptsA#j)) > 0 then continue;
      if wA#i + wA#j > 8 then sig = sig + 1;
      h := crossv(PA#i, PA#j);
      if seenL#?h then continue; seenL#h = true;
      W := sum(select(toList(0..#ptsA-1), k -> dot3(h, PA#k) == 0), k -> wA#k);
      if W == 8 then bad8 = bad8 + 1;
      if W >= 9 then rich = rich + 1);
    if sig > 0 then reasons = append(reasons, "unjoined pair sigma>8");
    if bad8 > 0 then reasons = append(reasons, "external line of weight 8 (" | toString bad8 | ")");
    if rich > 4 then reasons = append(reasons, "too many rich lines");
    print("   A2 = B(" | toString T | ") minus line " | toString a | ": mults " | toString tally apply(ptsA, S -> #S) | ", s_H " | toString tally sA
          | ", rich lines " | toString rich | "  -> " | (if #reasons == 0 then "SURVIVES" else toString reasons));
    if #reasons == 0 then surv = surv + 1));
print("surviving A2 candidates: " | toString surv);
