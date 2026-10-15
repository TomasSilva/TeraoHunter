-- 4-line deletions of A(G26) = Hesse u dual Hesse (21 lines, over Q(w), w^2-w+1=0) with b2 = 80: window test and freeness
K = toField(QQ[w]/(w^2-w+1)); use K;
S = K[x,y,z]; T = K[u,v];
Af = apply({0, 1, w-1, w}, a -> y - a*z); Bf = apply({0, 1, w, 1-w}, b -> x - b*z);
Cf = apply({0, 1, w, w+1}, c -> x + y - c*z); Df = apply({0, 1, w-1, w}, d -> x + w*y - d*z);
rich = {x + (1-w)*y - z, x + ((w+1)/3)*y - ((w+1)/3)*z, x + (w-1)*y, x + (w+1)*y - w*z};
forms = {z} | Af | Bf | Cf | Df | rich;
n = #forms;
vec = g -> {coefficient(x,g), coefficient(y,g), coefficient(z,g)};
vs = apply(forms, vec);
pts = unique apply(subsets(n,2), p -> sort unique({p#0,p#1} | select(toList(0..n-1), k -> not member(k,p) and det matrix{vs#(p#0),vs#(p#1),vs#k} == 0)));
print("A(G26): points " | toString tally apply(pts, P -> #P) | ", b2 = " | toString sum(pts, P -> #P-1));
-- lattice of a sub-arrangement given by kept line indices
subLat = keep -> (ks := set keep; select(apply(pts, P -> select(P, i -> member(i, ks))), P -> #P >= 2));
cr = (p,q) -> {p#1*q#2-p#2*q#1, p#2*q#0-p#0*q#2, p#0*q#1-p#1*q#0};
-- Yoshinaga test on a line a of the sub-arrangement 'keep': true iff a derivation of degree 7 exists on the Ziegler restriction (=> not free)
jump7 = (keep, a) -> (lat := subLat keep; onH := select(lat, P -> member(a,P));
  k := position({0,1,2}, j -> (vs#a)#j != 0); ij := delete(k,{0,1,2});
  forms2 := apply(onH, P -> (b := first select(P, i -> i != a); Pt := cr(vs#a, vs#b); ((Pt#(ij#1))*u - (Pt#(ij#0))*v, #P - 1)));
  e := 7; mons := flatten entries basis(e, T); unk := #mons;
  rows := flatten apply(forms2, fm -> (al := fm#0; m := fm#1; cs := {coefficient(u,al), coefficient(v,al)}; k0 := min(m-1, e);
     apply(k0+1, aa -> apply(2*unk, j -> (F := if j < unk then cs#0*mons#j else cs#1*mons#(j-unk); Dd := diff(u^aa*v^(k0-aa), F); sub(Dd, {u => -(cs#1), v => cs#0}))))));
  rank matrix apply(rows, r -> apply(r, c -> sub(c,K))) < 2*unk);
cnt = 0; freeList = {}; nonfree = 0;
for del in subsets(n, 4) do (
  keep := toList(set(0..n-1) - set del);
  lat := subLat keep;
  b2 := sum(lat, P -> #P - 1);
  if b2 != 80 then continue;
  sH := apply(keep, a -> #select(lat, P -> member(a,P)));
  mx := max apply(lat, P -> #P);
  if min sH < 4 or max sH > 8 or mx > 6 then continue;
  cnt = cnt + 1;
  a := keep#0;
  nf := jump7(keep, a);
  if nf then nonfree = nonfree + 1 else (freeList = append(freeList, del); print("FREE (8,8) window arrangement: delete lines " | toString del | ", points " | toString tally apply(lat, P -> #P) | ", s_H " | toString tally sH)));
print("4-line deletions with b2 = 80 in the window: " | toString cnt | ", non-free: " | toString nonfree | ", free: " | toString(#freeList));
