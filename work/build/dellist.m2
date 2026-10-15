-- list the 4-line deletions with b2 = 80 (window) of A(G26) and A(6,1,3), with the lattice of the ambient arrangement
K = toField(QQ[w]/(w^2-w+1)); use K; S = K[x,y,z];
vec = g -> {coefficient(x,g), coefficient(y,g), coefficient(z,g)};
doIt = (name, forms) -> (
  n := #forms; vs := apply(forms, vec);
  pts := unique apply(subsets(n,2), p -> sort unique({p#0,p#1} | select(toList(0..n-1), k -> not member(k,p) and det matrix{vs#(p#0),vs#(p#1),vs#k} == 0)));
  fo := openOut(name | "_lattice.txt");
  for P in pts do fo << "point " << demark(" ", apply(P, toString)) << endl;
  close fo;
  fd := openOut(name | "_dels.txt"); cnt := 0;
  for del in subsets(n,4) do (
    ks := set(toList(set(0..n-1) - set del));
    lat := select(apply(pts, P -> select(P, i -> member(i, ks))), P -> #P >= 2);
    if sum(lat, P -> #P - 1) != 80 then continue;
    keep := sort toList ks;
    sH := apply(keep, a -> #select(lat, P -> member(a,P))); mx := max apply(lat, P -> #P);
    if min sH < 4 or max sH > 8 or mx > 6 then continue;
    cnt = cnt + 1; fd << demark(" ", apply(del, toString)) << endl);
  close fd; print(name | ": " | toString cnt | " deletions written"));
Af = apply({0, 1, w-1, w}, a -> y - a*z); Bf = apply({0, 1, w, 1-w}, b -> x - b*z);
Cf = apply({0, 1, w, w+1}, c -> x + y - c*z); Df = apply({0, 1, w-1, w}, d -> x + w*y - d*z);
rich = {x + (1-w)*y - z, x + ((w+1)/3)*y - ((w+1)/3)*z, x + (w-1)*y, x + (w+1)*y - w*z};
doIt("g26", {z} | Af | Bf | Cf | Df | rich);
mu6 = apply(6, k -> w^k);
doIt("fermat", {x, y, z} | apply(mu6, u -> x - u*y) | apply(mu6, u -> y - u*z) | apply(mu6, u -> z - u*x));
