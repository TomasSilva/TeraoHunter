-- the nearly free 13-line arrangement A(4,1,3) minus {x, y-z}: locate its Bourbaki (jumping) point and compute the Ziegler
-- exponents of every line (restriction multiarrangement), to check "lines through P0 <-> exponents (d1'-1, d2')".
kk = ZZ/10009; R = kk[x,y,z];
i4 = (for a from 2 to 10008 do (if (a^4 % 10009 == 1) and (a^2 % 10009 != 1) then break a)); i4 = sub(i4, kk);
A4 = {x, y, z} | apply(4, k -> x - i4^k*y) | apply(4, k -> y - i4^k*z) | apply(4, k -> z - i4^k*x);
lns = delete(y - z, delete(x, A4)); n = #lns; Q = product lns;
jac = matrix{{diff(x,Q), diff(y,Q), diff(z,Q)}}; S = mingens image syz jac;
degs = apply(numcols S, i -> first (degrees source S)#i - (n-1)); r = min degs; i0 = position(degs, d -> d == r);
th = S_{i0}; a = th_(0,0); b = th_(1,0); c = th_(2,0);
hs = for k from 0 to numcols S - 1 list (e := S_{k}; (b*e_(2,0) - c*e_(1,0)) // diff(x,Q));
I = saturate ideal hs;
<< "syz degrees " << sort degs << "; Bourbaki ideal: " << toString mingens I << " degree " << degree I << endl;
P0 = radical I;
<< "lines through P0: " << select(lns, l -> isSubset(ideal l, P0)) << endl;
-- Ziegler exponents of each line H: restrict the other lines to H with multiplicities, compute exponents of the multiarrangement
zieglerExp = (H) -> (
  others := delete(H, lns);
  -- coordinates on H: eliminate a variable with nonzero coefficient
  vars3 := {x,y,z}; k := position(vars3, v -> coefficient(v, H) != 0); vk := vars3#k; rest := delete(vk, vars3);
  invc := 1/(coefficient(vk, H));
  forms := apply(others, L -> sub(L, {vk => -(H - coefficient(vk,H)*vk) * invc}));
  U := rest#0; V := rest#1;
  pq := apply(forms, f -> (coefficient(U,f), coefficient(V,f)));
  groups := {};
  scan(#pq, i -> (j := position(groups, G -> (w := pq#(G#0); w#0*pq#i#1 - w#1*pq#i#0 == 0));
     if j === null then groups = append(groups, {i}) else groups = replace(j, append(groups#j, i), groups)));
  T := kk[symbol s, symbol t]; s := T_0; t := T_1;
  -- derivation module: (f,g) with (p s + q t)^mu | p f + q g ; compute exponents as the degrees where the kernel first appears
  N := n - 1;
  expo := {};
  for e from 0 to N do (
    mons := flatten entries basis(e, T); nm := #mons;
    rows := {};
    scan(groups, G -> (pp := pq#(G#0)#0; qq := pq#(G#0)#1; mu := #G;
      k0 := min(mu-1, e); for i from 0 to k0 do rows = append(rows, flatten {apply(mons, mo -> sub(pp*diff(s^i*t^(k0-i), mo), {s => qq, t => -pp})), apply(mons, mo -> sub(qq*diff(s^i*t^(k0-i), mo), {s => qq, t => -pp}))})));
    M := matrix rows; kdim := 2*nm - rank M;
    expo = append(expo, kdim));
  -- exponents: e1 = first e with kdim > 0; e2 = N - e1
  e1 := position(expo, k -> k > 0);
  (e1, N - e1, #groups));
scan(lns, H -> << H << " : s_H=" << (zieglerExp H)#2 << " exponents " << (zieglerExp H)_{0,1} << (if isSubset(ideal H, P0) then "  (through P0)" else "") << endl);
