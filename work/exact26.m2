-- 26-plane P^3 arrangement over QQ(w), w primitive 6th root of unity: planes spanned by triples of
-- e1, e2, e1+a e2, X1+X2, X1+c X2, e3, e4, e3+b e4 (X1 = e1+e2, X2 = e3+e4), for all (a,b,c) in {w, w^5}^3.
needsPackage "HyperplaneArrangements";
K = toField(QQ[w]/(w^2 - w + 1));
R = K[x_1..x_4];
wb = 1 - w;      -- = w^5, the other primitive 6th root
e = apply(4, i -> apply(4, j -> if i == j then 1_K else 0_K));
add = (u, v) -> apply(4, i -> u#i + v#i);
scl = (c, u) -> apply(4, i -> c * u#i);
X1 = add(e#0, e#1); X2 = add(e#2, e#3);
normalize = v -> (j := position(v, c -> c != 0); apply(v, c -> c / v#j));
for abc in toList(set{w, wb} ** set{w, wb} ** set{w, wb}) do (
  (ab, c) := toSequence abc; (a, b) := toSequence ab;
  pts := {e#0, e#1, add(e#0, scl(a, e#1)), add(X1, X2), add(X1, scl(c, X2)), e#2, e#3, add(e#2, scl(b, e#3))};
  planes := {};
  for tri in subsets(8, 3) do (
     M := matrix apply(tri, i -> pts#i);
     if rank M == 3 then (
        n := flatten entries gens ker M;
        v := normalize flatten entries (gens ker M);
        planes = append(planes, v)));
  planes = unique planes;
  forms := apply(planes, v -> sum(4, i -> v#i * x_(i+1)));
  A := arrangement forms;
  D := prune image der A;
  print("(a,b,c) = " | toString(a, b, c) | ": " | toString(#planes) | " planes, free: " | toString isFreeModule D | ", degrees " | toString sort flatten degrees D);
);
