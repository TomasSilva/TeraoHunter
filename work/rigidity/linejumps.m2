-- Per-line jump functions of a one-parameter family of line arrangements (balanced case d1 = d2).
-- Family over P^1 with parameter (l0:l1).  For every line H the Ziegler multirestriction
-- (A^H, m^H) has |m^H| = n-1 = 2d; its exponents are (d,d) iff the 2d x 2d matrix of degree-(d-1)
-- derivations is nonsingular.  ThetaH = gcd over the three coordinate charts of its determinant.
-- Every irreducible factor of ThetaH is classified as
--   collision : a new concurrency through H (two points of H merge),
--   remote    : a new concurrency away from H (the lattice coarsens, H's points stay distinct),
--   interior  : neither -- a non-free member with the generic lattice (would contradict TC).
R = QQ[l0,l1];
prim = v -> (g := fold(gcd, v); if g == 0 then v else apply(v, c -> c // g));
crossv = (u,v) -> {u#1*v#2-u#2*v#1, u#2*v#0-u#0*v#2, u#0*v#1-u#1*v#0};
det3 = (a,b,c) -> det matrix{a,b,c};
ff = (x,k) -> product(k, r -> x - r);
monic = q -> q * (1/(leadCoefficient q));
irr = F -> select(apply(toList factor F, p -> {monic(p#0), p#1}), p -> first degree(p#0) > 0);
mult = (F,q) -> (c := 0; G := F; while G % q == 0 do (G = G // q; c = c + 1); c);

analyseFamily = (name, H) -> (
  n := #H;
  dd := (n-1)//2;  e := dd - 1;
  print("== " | name | ": " | toString n | " lines, balanced exponents (" | toString dd | "," | toString dd | ")");
  -- generic concurrencies
  trip := subsets(n,3);
  D3 := new MutableHashTable;
  for t in trip do D3#t = det3(H#(t#0), H#(t#1), H#(t#2));
  pts := unique apply(subsets(n,2), p -> sort unique({p#0,p#1} | select(toList(0..n-1), k -> not member(k,p) and D3#(sort{p#0,p#1,k}) == 0)));
  mlt := tally apply(pts, S -> #S);
  print("generic lattice: points by multiplicity " | toString mlt);
  b2 := sum(pts, S -> #S - 1);
  print("b2 = " | toString b2 | ",  chi0(0) = b2-n+1 = " | toString(b2-n+1) | "  (d^2 = " | toString(dd^2) | ")");
  -- degeneration factors: triples that are not generically concurrent
  degF := new MutableHashTable;   -- factor -> list of triples
  for t in trip do if D3#t != 0 then for p in irr(D3#t) do (
      q := p#0; degF#q = (if degF#?q then degF#q else {}) | {t});
  print("degeneration factors (new concurrencies): " | toString(#keys degF));
  for q in sort keys degF do print("   " | toString q | "   in " | toString(#(degF#q)) | " triples");
  results := {};
  for a from 0 to n-1 do (
    onH := select(pts, S -> member(a, S));
    dets := {};
    for k from 0 to 2 do if (H#a)#k != 0 then (
      ij := delete(k, {0,1,2});
      rows := flatten apply(onH, S -> (
          b := first select(S, x -> x != a);
          P := prim crossv(H#a, H#b);
          Pi := P#(ij#0); Pj := P#(ij#1);
          m := #S - 1;
          apply(m, aa -> (
             Dd := u -> if u >= aa and (e-u) >= (m-1-aa) then ff(u,aa)*ff(e-u,m-1-aa)*Pi^(u-aa)*Pj^(e-u-(m-1-aa)) else 0_R;
             apply(e+1, u -> Pj*Dd(u)) | apply(e+1, u -> -Pi*Dd(u))))));
      M := matrix rows;
      dets = dets | {det M});
    Th := fold(gcd, dets);
    sH := #onH;
    coll := select(keys degF, q -> any(degF#q, t -> member(a, t)));
    remo := select(keys degF, q -> not member(q, coll));
    fs := irr Th;
    typ := q -> if member(q, coll) then "collision" else if member(q, remo) then "remote" else "INTERIOR";
    interior := select(fs, p -> typ(p#0) == "INTERIOR");
    missed := select(remo, q -> Th % q != 0);
    print("-- line " | toString(a+1) | ": s_H = " | toString sH | ", multiplicities " | toString sort apply(onH, S -> #S-1)
          | ", deg ThetaH = " | toString(first degree Th));
    for p in fs do print("     " | toString(p#0) | " ^" | toString(p#1) | "   " | typ(p#0));
    print("     remote degenerations not in ThetaH (boundary-jump lemma predicts none): " | toString missed);
    print("     interior factors (TC violations on the family): " | toString interior);
    results = results | {(a+1, first degree Th, #interior, #missed)});
  results);
