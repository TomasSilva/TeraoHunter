-- Specialise a family at an algebraic parameter value and report: lattice (b2), exponents of the
-- Ziegler multirestriction on every line, and extra collinearities of intersection points
-- (lines outside A through >= 3 multiple points).
crossK = (u,v) -> {u#1*v#2-u#2*v#1, u#2*v#0-u#0*v#2, u#0*v#1-u#1*v#0};
normK = v -> (i := position(v, c -> c != 0); apply(v, c -> c / v#i));
analyseAt = (name, HK) -> (
  n := #HK; K := ring first first HK;
  trip := subsets(n,3);
  conc := t -> det matrix{HK#(t#0), HK#(t#1), HK#(t#2)} == 0;
  pts := unique apply(subsets(n,2), p -> sort unique({p#0,p#1} | select(toList(0..n-1), k -> not member(k,p) and conc(sort{p#0,p#1,k}))));
  b2 := sum(pts, S -> #S - 1);
  print("== " | name | ": points " | toString tally apply(pts, S -> #S) | ", b2 = " | toString b2 | ", chi0(0) = " | toString(b2-n+1));
  -- multirestriction exponents: smallest e with a nonzero derivation of degree e
  S := K[s,t];
  for a from 0 to n-1 do (
    onH := select(pts, S0 -> member(a, S0));
    k := position({0,1,2}, j -> (HK#a)#j != 0);
    ij := delete(k, {0,1,2});
    forms := apply(onH, S0 -> (b := first select(S0, x -> x != a); P := crossK(HK#a, HK#b);
                                ((P#(ij#1))*s - (P#(ij#0))*t, #S0 - 1)));
    tot := sum(forms, f -> f#1);
    e1 := first select(toList(0..tot), e -> (
        mons := flatten entries basis(e, S);
        if #mons == 0 then false else (
        unk := #mons;
        rows := flatten apply(forms, fm -> (al := fm#0; m := fm#1;
             cs := {coefficient(s, al), coefficient(t, al)};
             -- theta = f d_s + g d_t; theta(al) = cs0*f + cs1*g must be divisible by al^m
             k0 := min(m-1, e);
             apply(k0+1, aa -> apply(2*unk, u -> (
                 F := (if u < unk then cs#0 * mons#u else cs#1 * mons#(u-unk));
                 D := diff(s^aa * t^(k0-aa), F);
                 -- evaluate at the root of al: (s,t) = (-cs1, cs0)
                 sub(D, {s => -(cs#1), t => cs#0}))))));
        Mx := matrix apply(rows, r -> apply(r, c -> sub(c, K)));
        rank Mx < 2*unk)));
    print("   line " | toString(a+1) | ": s_H = " | toString(#onH) | ", exponents (" | toString e1 | "," | toString(tot-e1) | ")"));
  -- extra collinearities among multiple points (m >= 3) and all points
  P := apply(pts, S0 -> normK crossK(HK#(S0#0), HK#(S0#1)));
  extra := {};
  for pr in subsets(#P, 2) do (
    L0 := normK crossK(P#(pr#0), P#(pr#1));
    if any(HK, h -> normK h == L0) then continue;
    on := select(toList(0..#P-1), i -> sum(3, j -> L0#j * (P#i)#j) == 0);
    if #on >= 3 then extra = append(extra, (L0, on)));
  extra = unique extra;
  mx := if #extra == 0 then 0 else max apply(extra, x -> #(x#1));
  print("   extra lines through >= 3 intersection points: " | toString(#extra) | ", max points on one: " | toString mx);
  for x in select(extra, x -> #(x#1) >= 6) do print("      line " | toString(x#0) | " through " | toString(#(x#1)) | " points, multiplicities " | toString apply(x#1, i -> #(pts#i)));
  );
