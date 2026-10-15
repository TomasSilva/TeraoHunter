-- points per line and multiplicities on each line, for the F_8 hit (lattice data)
offs = mask -> select(toList(0..7), v -> (mask >> v) % 2 == 1);
R0 = ZZ/2[t]; F8 = toField(R0/ideal(t^3+t+1)); tt = sub(t, F8);
el8 = v -> sum(3, i -> if (v >> i) % 2 == 1 then tt^i else 0_F8);
lines8 = S -> {{0_F8, 0_F8, 1_F8}} | flatten apply(4, j -> apply(offs S#j, v ->
      if j == 0 then {0_F8, 1_F8, el8 v} else if j == 1 then {1_F8, 0_F8, el8 v}
      else if j == 2 then {1_F8, 1_F8, el8 v} else {1_F8, tt, el8 v}));
normalize = l -> (c := first select(l, e -> e != 0); apply(l, e -> e / c));
cr = (u, v) -> normalize {u#1*v#2 - u#2*v#1, u#2*v#0 - u#0*v#2, u#0*v#1 - u#1*v#0};
lineProfile = S -> (L := lines8 S; n := #L;
   apply(n, i -> (pts := new MutableHashTable;
      for j from 0 to n - 1 do if j != i then (P := cr(L#i, L#j); pts#P = (if pts#?P then pts#P else 0) + 1);
      sort values pts)));
for S in {{0xfe, 0xc1, 0xcc, 0xf0}, {0xfe, 0xc1, 0x3c, 0x3c}} do (
   pr := lineProfile S;
   print(toString S | ": per line (multiplicities m(X) of its points): ");
   scan(#pr, i -> print("   line " | toString i | ": s = " | toString(#(pr#i)) | "  " | toString(pr#i))));
