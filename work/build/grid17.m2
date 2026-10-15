-- The 17-line four-pencil arrangement: z=0, y=az (a in A), x=bz (b in B), x-y=cz (c in C), x+y=dz (d in D)
S = QQ[x,y,z];
A = {0,1,2,3}; B = {2,3,4,5}; C = {0,1,2,3}; D = {3,4,5,6};
forms = {z} | apply(A, a -> y - a*z) | apply(B, b -> x - b*z) | apply(C, c -> x - y - c*z) | apply(D, d -> x + y - d*z);
n = #forms; print("lines: " | toString n);
f = product forms;
J = ideal(diff(x,f), diff(y,f), diff(z,f));
Sy = syz gens J;
print("Jacobian syzygy degrees: " | toString sort apply(numcols Sy, i -> first degree Sy_(0,i)));
-- lattice
vecs = apply(forms, g -> {coefficient(x,g), coefficient(y,g), coefficient(z,g)});
pts = unique apply(subsets(n,2), p -> sort unique({p#0,p#1} | select(toList(0..n-1), k -> not member(k,p) and det matrix{vecs#(p#0),vecs#(p#1),vecs#k} == 0)));
print("points by multiplicity: " | toString tally apply(pts, P -> #P) | ", b2 = " | toString sum(pts, P -> #P-1));
print("points per line: " | toString tally apply(n, a -> #select(pts, P -> member(a,P))));
-- external lines through unjoined pairs: weights
cr = (u,v) -> {u#1*v#2-u#2*v#1, u#2*v#0-u#0*v#2, u#0*v#1-u#1*v#0};
nrm = v -> (g := gcd v; apply(v, c -> c//g));
coords = apply(pts, P -> nrm cr(vecs#(P#0), vecs#(P#1)));
w = apply(pts, P -> #P - 1);
tal = new MutableHashTable; seen = new MutableHashTable;
for pr in subsets(#pts,2) do (
  i := pr#0; j := pr#1;
  if #(set(pts#i)*set(pts#j)) > 0 then continue;
  h := nrm cr(coords#i, coords#j);
  if seen#?h then continue; seen#h = true;
  W := sum(select(toList(0..#pts-1), k -> sum(3, t -> h#t * (coords#k)#t) == 0), k -> w#k);
  tal#W = (if tal#?W then tal#W else 0) + 1);
print("external lines by weight (free (8,8) needs all <= 8): " | toString sort pairs new HashTable from tal);
