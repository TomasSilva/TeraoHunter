-- Check: for theta the minimal derivation and every line H of A,
--   (J / alpha_H)|_H  =  c * (s G - t F) * Jac(F, G),
-- where Phi_H = [F : G] is the line map of theta_H (linemaps.m2) and J = det(E, theta, theta(theta)) is the
-- inflection curve. Equivalently J'|_H = c * W * Jac / prod_P l_P^(m_P - 1): the residual inflection curve
-- meets H exactly in the excess local degrees, the extra fixed points and the free critical points of Phi_H.
S = QQ[x,y,z]; B = QQ[s,t];
pen = (dp, dq, offs) -> apply(offs, v -> dq*x - dp*y - v*z);
arrs = {
 ("R17a", {z} | pen(1,0,{0,1,2,3}) | pen(0,1,{0,2,4}) | pen(1,1,{2,3,5}) | pen(4/3,1,{4}) | pen(2,1,{2,4,5,6,8})),
 ("R17b", {z} | pen(1,0,{0,1,2,3,4}) | pen(0,1,{0,1,2,4}) | pen(1,1,{4}) | pen(-1/2,1,{1}) | pen(-1,1,{-2,-1,0,1,2})),
 ("B3 (free 3,5)", {x,y,z,x-y,x+y,x-z,x+z,y-z,y+z}),
 ("random 10 lines", {x,y,z,x+y+z,x-2*y+3*z,2*x+5*y-z,x+3*y+7*z,4*x-y+2*z,x+y-5*z,3*x+2*y+z})};
for A in arrs do (
  L := A#1; f := product L;
  Sy := syz matrix{{diff(x,f), diff(y,f), diff(z,f)}};
  cols := entries transpose Sy; dg := apply(cols, c -> first degree c#0);
  th := cols#(position(dg, d -> d == min dg)); r := min dg;
  Xf := h -> th#0*diff(x,h) + th#1*diff(y,h) + th#2*diff(z,h);
  J := det matrix{{x,y,z}, th, apply(th, Xf)};
  allok := true;
  for H in L do (
    cf := apply({x,y,z}, v -> lift(coefficient(v, H), QQ));
    KH := (th#0*cf#0 + th#1*cf#1 + th#2*cf#2) // H;
    thH := apply(3, i -> th#i - KH*({x,y,z}#i));
    kerB := gens kernel matrix{cf};
    v1 := flatten entries kerB_{0}; v2 := flatten entries kerB_{1};
    subst := map(B, S, apply(3, i -> s*v1#i + t*v2#i));
    vals := apply(thH, g -> subst g);
    Vm := matrix{v1, v2};
    idx := first select(subsets(3,2), p -> det submatrix(Vm, , p) != 0);
    Mi := inverse (submatrix(Vm, , idx));
    FG := flatten entries (matrix{{vals#(idx#0), vals#(idx#1)}} * Mi);
    F := FG#0; G := FG#1;
    W := s*G - t*F;
    Jac := det matrix{{diff(s,F), diff(t,F)}, {diff(s,G), diff(t,G)}};
    lhs := subst (J // H);
    rhs := W * Jac;
    -- proportional?
    ok := (lhs == 0 and rhs == 0) or (rhs != 0 and lhs != 0 and (
           c := leadCoefficient(lhs) / leadCoefficient(rhs); lhs == c * rhs));
    if not ok then (allok = false; print("  FAIL on line " | toString H)));
  print(A#0 | ": r = " | toString r | ", deg J = " | toString(first degree J) | ", identity holds on all " |
        toString(#L) | " lines: " | toString allok));
