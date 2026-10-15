-- The minimal derivation of a non-free arrangement, seen on each line as a rational self-map.
-- For theta in D_0(A) of least degree r and a line H = {alpha = 0}: theta_H = theta - (theta(alpha)/alpha) E
-- is tangent to the plane H; in a basis v1, v2 of that plane theta_H(s v1 + t v2) = F v1 + G v2, and
-- Phi_H = [F : G] : H -> H. Claim (TC_GLOBAL.md, Theorem A): every multiple point P of H is a fixed point
-- of Phi_H with local degree >= m_P - 1 (superattracting when m_P >= 3), and the Camacho-Sad index of H at P
-- equals the holomorphic fixed-point index 1/(1 - Phi'(P)).
S = QQ[x,y,z];
pen = (dp, dq, offs) -> apply(offs, v -> dq*x - dp*y - v*z);
arrs = {
 ("R17a", {z} | pen(1,0,{0,1,2,3}) | pen(0,1,{0,2,4}) | pen(1,1,{2,3,5}) | pen(4/3,1,{4}) | pen(2,1,{2,4,5,6,8})),
 ("R17b", {z} | pen(1,0,{0,1,2,3,4}) | pen(0,1,{0,1,2,4}) | pen(1,1,{4}) | pen(-1/2,1,{1}) | pen(-1,1,{-2,-1,0,1,2}))};
B = QQ[s,t];
for A in arrs do (
  L := A#1; f := product L;
  Sy := syz matrix{{diff(x,f), diff(y,f), diff(z,f)}};
  cols := entries transpose Sy; dg := apply(cols, c -> first degree c#0);
  th := cols#(position(dg, d -> d == min dg)); r := min dg;
  P0 := saturate ideal th;
  print("==== " | A#0 | "  r = " | toString r | "  Bourbaki ideal " | toString P0);
  for H in L do (
    cf := apply({x,y,z}, v -> lift(coefficient(v, H), QQ));
    KH := (th#0*cf#0 + th#1*cf#1 + th#2*cf#2) // H;
    thH := apply(3, i -> th#i - KH*({x,y,z}#i));
    -- basis of the plane alpha = 0
    kerB := gens kernel matrix{cf};   -- 3x2
    v1 := flatten entries kerB_{0}; v2 := flatten entries kerB_{1};
    subst := map(B, S, apply(3, i -> s*v1#i + t*v2#i));
    vals := apply(thH, g -> subst g);
    -- solve vals = F v1 + G v2 (exact: choose two independent coordinates)
    Vm := matrix{v1, v2}; -- 2x3
    idx := first select(subsets(3,2), p -> det submatrix(Vm, , p) != 0);
    Mi := inverse (submatrix(Vm, , idx));   -- 2x2 over QQ
    rhs := matrix{{vals#(idx#0), vals#(idx#1)}};
    FG := flatten entries (rhs * Mi);
    F := FG#0; G := FG#1;
    -- check consistency on the third coordinate
    other := first toList(set{0,1,2} - set idx);
    ok := (F*v1#other + G*v2#other) == vals#other;
    h := gcd(F, G);
    -- multiple points on H: intersections with other lines, with multiplicities
    pts := new MutableHashTable;
    for H2 in L do if H2 != H then (
        I := ideal(H, H2); pt := flatten entries gens kernel matrix{cf, apply({x,y,z}, v -> lift(coefficient(v,H2), QQ))};
        -- coordinates (s,t) of the point: solve pt = s v1 + t v2
        st := flatten entries (matrix{apply(idx, i -> pt#i)} * Mi);
        g := gcd(numerator(st#0), numerator(st#1)); -- normalize projectively
        key := if st#1 != 0 then toString(st#0/st#1) else "inf";
        if pts#?key then pts#key = (pts#key#0, pts#key#1 + 1) else pts#key = (st, 1));
    W := s*G - t*F;               -- fixed points
    Jac := det matrix{{diff(s,F), diff(t,F)}, {diff(s,G), diff(t,G)}};   -- ramification
    report := {};
    totIdx := 0;
    for key in keys pts do (
        st := pts#key#0; mu := pts#key#1;  -- mu = number of other lines through P = m_P - 1
        lf := st#1*s - st#0*t;    -- linear form vanishing at P
        ordW := 0; tmp := t*F - s*G; -- placeholder
        -- order of vanishing of (t_P F - s_P G) at P
        Wp := st#1*F - st#0*G; o := 0; while Wp != 0 and Wp % lf == 0 do (Wp = Wp // lf; o = o + 1);
        oj := 0; Jt := Jac; while Jt != 0 and Jt % lf == 0 do (Jt = Jt // lf; oj = oj + 1);
        report = append(report, (key, mu, o, oj)));
    print("  line " | toString H | ": s_H = " | toString(#keys pts) | ", deg Phi = " | toString(first degree F - first degree h) |
          ", gcd deg " | toString(first degree h) | ", consistent " | toString ok);
    print("     (point, m_P - 1, ord(F - pG), ord Jac) = " | toString sort report);
    -- residual ramification (free critical points)
    Jr := Jac; for e in report do (lf := (pts#(e#0)#0)#1*s - (pts#(e#0)#0)#0*t; for i from 1 to e#3 do Jr = Jr // lf);
    print("     free critical points: " | toString factor Jr);
  ));
