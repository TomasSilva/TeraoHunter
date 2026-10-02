-- 26-plane P^3 arrangement over QQ(w), w a primitive 6th root of unity.  Saito check with exponents (1,7,9,9):
-- D(A)_k computed by exact linear algebra over the number field (substitution on each plane).
K = toField(QQ[w]/(w^2 - w + 1));
wb = 1 - w;
R = K[x_1..x_4];
e = apply(4, i -> apply(4, j -> if i == j then 1_K else 0_K));
add = (u, v) -> apply(4, i -> u#i + v#i);
scl = (c, u) -> apply(4, i -> c * u#i);
X1 = add(e#0, e#1); X2 = add(e#2, e#3);
normalize = v -> (j := position(v, c -> c != 0); apply(v, c -> c / v#j));
-- derivations theta = sum f_i d_i of degree k with theta(alpha) in (alpha) for all planes:
-- substitute a parametrisation of each plane and require theta(alpha) to vanish there.
derSpace = (planes, k) -> (
   T := K[s_1..s_3];
   mons := flatten entries basis(k, R);
   cols := flatten for i from 0 to 3 list for mu in mons list (
      flatten for P in planes list (
         al := P#0; B := P#1;          -- normal vector, 3 basis vectors of the plane
         sb := map(T, R, apply(4, c -> sum(3, j -> s_(j+1) * (B#j)#c)));
         F := al#i * sb(mu);
         tm := flatten entries basis(k, T);
         apply(tm, t -> coefficient(t, F))));
   M := transpose matrix cols;
   N := gens ker M;
   apply(numcols N, c -> (v := flatten entries N_{c}; apply(4, i -> sum(#mons, j -> v#(i*#mons + j) * mons#j)))));
for abc in toList(set{w, wb} ** set{w, wb} ** set{w, wb}) do (
  (ab, c) := toSequence abc; (a, b) := toSequence ab;
  pts := {e#0, e#1, add(e#0, scl(a, e#1)), add(X1, X2), add(X1, scl(c, X2)), e#2, e#3, add(e#2, scl(b, e#3))};
  planes := {};
  for tri in subsets(8, 3) do (
     M := matrix apply(tri, i -> pts#i);
     if rank M == 3 then planes = append(planes, normalize flatten entries gens ker M));
  planes = unique planes;
  -- plane data: normal and a basis of the plane (kernel of the normal)
  pd := apply(planes, n -> {n, apply(numcols gens ker matrix{n}, j -> flatten entries (gens ker matrix{n})_{j})});
  t0 := cpuTime();
  E := {x_1, x_2, x_3, x_4};
  D7 := derSpace(pd, 7); D9 := derSpace(pd, 9);
  -- dims of the free module S(-1)+S(-7)+S(-9)^2 in degrees 7 and 9: C(9,3)+1 = 85 and C(11,3)+C(5,3)+2 = 177
  -- random combinations: det(E, theta7, theta9, theta9') != 0 for generic choices iff A is free
  rnd := L -> (cs := apply(#L, i -> (random(-3,3) + random(-3,3)*w)); apply(4, c -> sum(#L, i -> cs#i * (L#i)#c)));
  ok := false;
  for trial from 1 to 3 do if not ok then (
     dt := det matrix{E, rnd D7, rnd D9, rnd D9};
     if dt != 0 then (ok = true; print("   Saito determinant nonzero, degree " | toString first degree dt)));
  print("(a,b,c) = " | toString(a, b, c) | ": " | toString(#planes) | " planes; dim D_7 = " | toString(#D7) | " (free: 85), dim D_9 = " | toString(#D9) | " (free: 177); FREE: " | toString ok | "  (" | toString(cpuTime() - t0) | "s)");
);
