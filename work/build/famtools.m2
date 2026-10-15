-- Tools for one-parameter families of line arrangements over R = QQ[l0,l1] (parameter lambda = l1/l0).
-- A family is a list of lines, each a list {a,b,c} of elements of R (homogeneous of the same degree).
R = QQ[l0,l1];
prim = v -> (g := fold(gcd, v); if g == 0 then v else apply(v, c -> c // g));
crossv = (u,v) -> prim {u#1*v#2-u#2*v#1, u#2*v#0-u#0*v#2, u#0*v#1-u#1*v#0};
det3 = (a,b,c) -> det matrix{a,b,c};
dot3 = (a,b) -> a#0*b#0 + a#1*b#1 + a#2*b#2;
ff = (x,k) -> product(k, r -> x - r);
monic = q -> q * (1/(leadCoefficient q));
irr = F -> select(apply(toList factor F, p -> {monic(p#0), p#1}), p -> first degree(p#0) > 0);
sameLine = (u,v) -> (crossv(u,v) == {0_R,0_R,0_R});

-- generic lattice: list of point-sets (lists of line indices), and the degeneration factors
latticeOf = H -> (
  n := #H;
  D3 := new MutableHashTable;
  for t in subsets(n,3) do D3#t = det3(H#(t#0), H#(t#1), H#(t#2));
  pts := unique apply(subsets(n,2), p -> sort unique({p#0,p#1} | select(toList(0..n-1), k -> not member(k,p) and D3#(sort{p#0,p#1,k}) == 0)));
  degF := new MutableHashTable;
  for t in subsets(n,3) do if D3#t != 0 then for p in irr(D3#t) do (q := p#0; degF#q = (if degF#?q then degF#q else {}) | {t});
  -- coincident lines
  for p in subsets(n,2) do (c := crossv(H#(p#0), H#(p#1)); for e in c do if e != 0 then for q in irr e do if all(c, x -> x % (q#0) == 0) then (degF#(q#0) = (if degF#?(q#0) then degF#(q#0) else {}) | {p}));
  (pts, degF));

-- window report for the balanced class (d,d), n = 2d+1
windowReport = (H, pts) -> (
  n := #H; d := (n-1)//2;
  b2 := sum(pts, S -> #S - 1);
  mlt := tally apply(pts, S -> #S);
  sH := apply(n, a -> #select(pts, S -> member(a,S)));
  ok := (b2 == d^2 + n - 1) and (max keys mlt <= d - 2) and (min sH >= 4) and (max sH <= d) and (member(4,sH) or member(5,sH) or member(d,sH));
  print("   points " | toString mlt | ", b2 = " | toString b2 | " (need " | toString(d^2+n-1) | "), points per line " | toString tally sH | ", window: " | toString ok);
  ok);

-- Ziegler jump function on line a (balanced case): gcd over charts of det of the (n-1)x(n-1) matrix of degree-(d-1) derivations
thetaH = (H, pts, a) -> (
  n := #H; d := (n-1)//2; e := d - 1;
  onH := select(pts, S -> member(a, S));
  dets := {};
  for k from 0 to 2 do if (H#a)#k != 0 then (
    ij := delete(k, {0,1,2});
    rows := flatten apply(onH, S -> (
        b := first select(S, x -> x != a);
        P := crossv(H#a, H#b);
        Pi := P#(ij#0); Pj := P#(ij#1);
        m := #S - 1;
        apply(m, aa -> (
           Dd := u -> if u >= aa and (e-u) >= (m-1-aa) then ff(u,aa)*ff(e-u,m-1-aa)*Pi^(u-aa)*Pj^(e-u-(m-1-aa)) else 0_R;
           apply(e+1, u -> Pj*Dd(u)) | apply(e+1, u -> -Pi*Dd(u))))));
    dets = dets | {det matrix rows});
  fold(gcd, dets));

-- classify the factors of Theta_H: collision (new concurrency through H), remote, or INTERIOR
classifyTheta = (H, pts, degF, a, Th) -> (
  coll := select(keys degF, q -> any(degF#q, t -> member(a, t)));
  remo := select(keys degF, q -> not member(q, coll));
  fs := irr Th;
  typ := q -> if member(q, coll) then "collision" else if member(q, remo) then "remote" else "INTERIOR";
  for p in fs do print("     " | toString(p#0) | " ^" | toString(p#1) | "   " | typ(p#0));
  select(fs, p -> typ(p#0) == "INTERIOR"));

-- hidden lines: for every unjoined pair of multiple points, the line through them, its generic weight,
-- and the extra weight it acquires at special parameters; report weight >= wmin (interior only)
hiddenLines = (H, pts, degF, wmin) -> (
  n := #H;
  P := apply(pts, S -> crossv(H#(S#0), H#(S#1)));
  w := apply(pts, S -> #S - 1);
  degKeys := keys degF;
  isDeg := q -> member(q, degKeys);
  found := {};
  seen := new MutableHashTable;
  for pr in subsets(#pts, 2) do (
    i := pr#0; j := pr#1;
    if #(set(pts#i) * set(pts#j)) > 0 then continue;
    h := crossv(P#i, P#j);
    key := h; if seen#?key then continue; seen#key = true;
    -- generic points on h and conditional ones
    gen := 0; cond := new MutableHashTable;
    for k from 0 to #pts-1 do (
      c := dot3(h, P#k);
      if c == 0 then gen = gen + w#k
      else for p in irr c do (q := p#0; cond#q = (if cond#?q then cond#q else 0) + w#k));
    if gen >= wmin then found = append(found, (h, gen, "generic", null));
    for q in keys cond do if gen + cond#q >= wmin and not isDeg q then found = append(found, (h, gen + cond#q, "at", q)));
  found);
