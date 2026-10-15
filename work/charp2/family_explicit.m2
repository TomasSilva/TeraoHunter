-- explicit generic member of a mu-flexible four-pencil lattice over F_2(mu), and its minimal syzygy
offs = mask -> select(toList(0..7), v -> (mask >> v) % 2 == 1);
explicitFamily = S -> (
   R0 := ZZ/2[t]; F8 := toField(R0/ideal(t^3+t+1)); tt := sub(t, F8);
   el8 := v -> sum(3, i -> if (v >> i) % 2 == 1 then tt^i else 0_F8);
   pen := flatten apply(4, j -> apply(offs S#j, v -> (j, v)));
   K := #pen;
   dir8 := j -> if j == 0 then {0_F8, 1_F8} else if j == 1 then {1_F8, 0_F8} else if j == 2 then {1_F8, 1_F8} else {1_F8, tt};
   L8 := {{0_F8,0_F8,1_F8}} | apply(pen, p -> dir8(p#0) | {el8(p#1)});
   C8 := select(subsets(K + 1, 3), s -> not member(0, s) and det matrix apply(s, i -> L8#i) == 0);
   A := frac(ZZ/2[mu]);
   dirA := j -> if j == 0 then {0_A, 1_A} else if j == 1 then {1_A, 0_A} else if j == 2 then {1_A, 1_A} else {1_A, mu_A};
   -- equation of a concurrent triple (i,j,k) (indices into L, 0 = H excluded): det of the 3x3 matrix, linear in offsets
   rows := apply(C8, c -> (
      apply(K, l -> (if member(l + 1, c) then (
          others := select(c, i -> i != l + 1);
          m := matrix {dirA(pen#(others#0 - 1)#0), dirA(pen#(others#1 - 1)#0)};
          -- cofactor of the offset entry of line l+1 in the 3x3 det: sign irrelevant in char 2
          det m) else 0_A))));
   M := matrix rows;
   N := gens ker M;
   T := matrix apply(K, i -> dirA(pen#i#0));
   ix := select(numcols N, j -> rank(T | N_{j}) == 3);
   v := N_{first ix};
   i0 := position(pen, p -> p#0 == 0); i1 := position(toList(0..K-1), l -> (pen#l)#0 == 0 and l != i0); i2 := position(pen, p -> p#0 == 1);
   beta := v_(i0, 0); alpha := v_(i2, 0);
   w := apply(K, l -> v_(l, 0) + (dirA(pen#l#0))#0 * alpha + (dirA(pen#l#0))#1 * beta);
   sc := w#i1; w = apply(w, e -> e / sc);
   print("lattice " | toString S | ": kernel dim " | toString(numcols N) | "; offsets over F_2(mu) by pencil (y=a, x=b, x+y=c, x+mu y=d):");
   scan(4, j -> print("   pencil " | toString(j+1) | ": " | toString apply(select(K, l -> pen#l#0 == j), l -> factor numerator w#l / factor denominator w#l)));
   -- minimal syzygy of the generic member
   P := A[x,y,z];
   f := z * product(K, l -> (dirA(pen#l#0))#0 * x + (dirA(pen#l#0))#1 * y + w#l * z);
   J := ideal(diff(x, f), diff(y, f), diff(z, f));
   Sy := gens ker matrix{{diff(x, f), diff(y, f), diff(z, f)}};
   degs := apply(numcols Sy, i -> first degree Sy_(0, i) );
   print("   syzygy generator degrees (first few): " | toString take(sort apply(numcols Sy, i -> (first degree (Sy_{i})_(0,0))), 5));
);
explicitFamily {0xfe, 0xc1, 0xcc, 0xf0};
explicitFamily {0xfe, 0xc1, 0x3c, 0x3c};
