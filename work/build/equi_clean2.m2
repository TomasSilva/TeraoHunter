K = toField(QQ[w]/(w^2-w+1)); use K;   -- w = e^{i pi/3}
S = K[x,y,z];
Af = apply({0, 1, w-1, w}, a -> y - a*z); Bf = apply({0, 1, w, 1-w}, b -> x - b*z);
Cf = apply({0, 1, w, w+1}, c -> x + y - c*z); Df = apply({0, 1, w-1, w}, d -> x + w*y - d*z);
forms = {z} | Af | Bf | Cf | Df;
vec = g -> {coefficient(x,g), coefficient(y,g), coefficient(z,g)};
latt = fs -> (n := #fs; vs := apply(fs, vec);
  unique apply(subsets(n,2), p -> sort unique({p#0,p#1} | select(toList(0..n-1), k -> not member(k,p) and det matrix{vs#(p#0),vs#(p#1),vs#k} == 0))));
syzd = fs -> (f := product fs; sort apply(numcols syz gens ideal(diff(x,f),diff(y,f),diff(z,f)), i -> first degree (syz gens ideal(diff(x,f),diff(y,f),diff(z,f)))_(0,i)));
-- Ziegler exponents of line a: smallest e with a derivation of degree e of the multiarrangement (clean variable names)
T = K[u,v];
ziegler = (fs, a) -> (pts := latt fs; vs := apply(fs, vec); onH := select(pts, P -> member(a,P));
  k := position({0,1,2}, j -> (vs#a)#j != 0); ij := delete(k,{0,1,2});
  forms2 := apply(onH, P -> (b := first select(P, i -> i != a); cr := (p,q) -> {p#1*q#2-p#2*q#1, p#2*q#0-p#0*q#2, p#0*q#1-p#1*q#0};
      Pt := cr(vs#a, vs#b); ((Pt#(ij#1))*u - (Pt#(ij#0))*v, #P - 1)));
  tot := sum(forms2, f -> f#1);
  first select(toList(1..tot), e -> (mons := flatten entries basis(e, T); unk := #mons;
     rows := flatten apply(forms2, fm -> (al := fm#0; m := fm#1; cs := {coefficient(u,al), coefficient(v,al)}; k0 := min(m-1, e);
        apply(k0+1, aa -> apply(2*unk, j -> (F := if j < unk then cs#0*mons#j else cs#1*mons#(j-unk); Dd := diff(u^aa*v^(k0-aa), F); sub(Dd, {u => -(cs#1), v => cs#0}))))));
     rank matrix apply(rows, r -> apply(r, c -> sub(c,K))) < 2*unk)));
report = (name, fs) -> (pts := latt fs; n := #fs;
  print("== " | name | ": " | toString n | " lines, points " | toString tally apply(pts, P -> #P) | ", b2 = " | toString sum(pts, P -> #P-1));
  print("   points per line: " | toString apply(n, a -> #select(pts, P -> member(a,P))));
  print("   Jacobian syzygy degrees: " | toString syzd fs);
  print("   Ziegler e1 per line: " | toString apply(n, a -> ziegler(fs, a))));
-- rich lines: from equi_verify (coefficient vectors in the (x,y,z) basis with the old sign convention y -> -y): recompute directly
f = product forms; J = ideal(diff(x,f),diff(y,f),diff(z,f)); Sy = syz gens J; cols = entries transpose Sy; dg = apply(cols, c -> first degree c#0); th = cols#(position(dg, d -> d == min dg));
IZ = saturate ideal th; print("   Bourbaki ideal: " | toString IZ);
pts = latt forms; vs = apply(forms, vec);
cr = (p,q) -> {p#1*q#2-p#2*q#1, p#2*q#0-p#0*q#2, p#0*q#1-p#1*q#0};
nrm = v -> (i := position(v, c -> c != 0); apply(v, c -> c/(v#i)));
coords = apply(pts, P -> nrm cr(vs#(P#0), vs#(P#1))); wt = apply(pts, P -> #P-1);
rich = {}; seen = new MutableHashTable;
for pr in subsets(#pts,2) do (i := pr#0; j := pr#1; if #(set(pts#i)*set(pts#j)) > 0 then continue; h := nrm cr(coords#i, coords#j); if seen#?h then continue; seen#h = true;
  W := sum(select(toList(0..#pts-1), k -> sum(3, s -> h#s*(coords#k)#s) == 0), k -> wt#k); if W >= 9 then rich = append(rich, (h, W)));
print("   rich lines (W >= 9): " | toString rich);
richForms = apply(rich, r -> (h := r#0; h#0*x + h#1*y + h#2*z));
U21 = forms | richForms;
report("E17 + rich lines", U21);
