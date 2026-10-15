-- Rich-line lemma check.  For an arrangement A (lines given as coefficient vectors over a field K):
--  * r = mdr(A) and a minimal Jacobian syzygy theta (a derivation of degree r killing Q);
--  * Z = zero scheme of theta (the Bourbaki scheme), of length c2(E(r)) = c2 - r(n-1) + r^2;
--  * every line l outside A through >= 2 multiple points ("diagonal"), its weight
--    W(l) = n - #(A meet l), and z(l) = length of Z meet l.
-- Claim (weight dichotomy): W(l) <= r - z(l)  or  W(l) = n - 1 - r + z(l)   [valid when 2r < n - 1].
-- Claim (leaves): when W(l) = n - 1 - r + z(l) the line l is invariant under theta.
-- Claim (count): the number of such "rich" lines is at most 3r - n.
normK = v -> (i := position(v, c -> c != 0); apply(v, c -> c / v#i));
crossK = (u,v) -> {u#1*v#2-u#2*v#1, u#2*v#0-u#0*v#2, u#0*v#1-u#1*v#0};
checkRich = (name, HK) -> (
  n := #HK; K := ring first first HK;
  S := K[x,y,z];
  forms := apply(HK, h -> h#0*x + h#1*y + h#2*z);
  f := product forms;
  M := matrix{{diff(x,f), diff(y,f), diff(z,f)}};
  Sy := syz M;
  degs := apply(numcols Sy, i -> first degree(Sy_(0,i)) );
  dmin := min degs;
  r := dmin;  -- entries of a syzygy of degree r
  cand := select(toList(0..numcols Sy - 1), i -> degs#i == dmin);
  th := Sy_(first cand);
  a := th_0; b := th_1; c := th_2;
  IZ := saturate ideal(a,b,c);
  lenZ := if IZ == ideal(1_S) then 0 else degree IZ;
  -- points of A with multiplicities
  conc := (i,j,k) -> det matrix{HK#i, HK#j, HK#k} == 0;
  pts := unique apply(subsets(n,2), p -> sort unique({p#0,p#1} | select(toList(0..n-1), k -> not member(k,p) and conc(p#0,p#1,k))));
  b2 := sum(pts, P -> #P - 1);
  c2 := b2 - n + 1;
  coords := apply(pts, P -> normK crossK(HK#(P#0), HK#(P#1)));
  print("== " | name | ": n = " | toString n | ", mdr r = " | toString r | ", #min generators = " | toString(#cand)
        | ", syzygy degrees " | toString sort degs);
  print("   points " | toString tally apply(pts, P -> #P) | ", b2 = " | toString b2 | ", c2 = " | toString c2
        | ", Bourbaki length = " | toString lenZ | " (predicted " | toString(c2 - r*(n-1) + r^2) | ")");
  if 2*r >= n-1 then print("   (2r >= n-1: the leaf claim does not apply)");
  -- diagonals
  joined := (i,j) -> any(subsets(n,1), h -> member(h#0, pts#i) and member(h#0, pts#j));
  diags := new MutableHashTable;
  for pr in subsets(#pts, 2) do (
    i := pr#0; j := pr#1;
    if #(set(pts#i) * set(pts#j)) > 0 then continue;   -- on a common line of A
    L0 := normK crossK(coords#i, coords#j);
    diags#L0 = 1);
  nrich := 0; bad := 0; leafbad := 0;
  wtally := new MutableHashTable;
  for L0 in keys diags do (
    meet := unique apply(HK, h -> normK crossK(L0, h));
    k := #meet; W := n - k;
    al := L0#0*x + L0#1*y + L0#2*z;
    zl := (J := saturate(IZ + ideal(al)); if J == ideal(1_S) then 0 else degree J);
    ok := (W <= r - zl) or (W == n - 1 - r + zl);
    key := (W, zl);
    wtally#key = (if wtally#?key then wtally#key else 0) + 1;
    if not ok then (bad = bad + 1; print("   VIOLATION: line " | toString L0 | " W = " | toString W | " z = " | toString zl));
    if W == n - 1 - r + zl and 2*r < n-1 then (
      nrich = nrich + 1;
      tv := a*diff(x,al) + b*diff(y,al) + c*diff(z,al);
      inv := (tv % ideal(al)) == 0;
      if not inv then leafbad = leafbad + 1;
      print("   rich line " | toString L0 | ": W = " | toString W | ", |A meet l| = " | toString k | ", z = " | toString zl | ", invariant: " | toString inv)));
  print("   diagonals: " | toString(#keys diags) | ", (W, z) tally " | toString sort pairs new HashTable from wtally);
  print("   dichotomy violations: " | toString bad | ";  rich lines: " | toString nrich | " (bound 3r - n = " | toString(3*r-n)
        | ");  rich lines not invariant: " | toString leafbad);
  );
