-- For a four-pencil hit (S1..S4 bitmasks over F_8 = F_2[t]/(t^3+t+1)):
--  (1) the hit itself over F_8 (mu = t): lattice and freeness;
--  (2) the deformation to a generic mu: solve the incidence equations M(mu) s = 0 over a big field,
--      take a kernel vector outside the translations, rebuild the 19 lines, compare the lattice
--      triple by triple with the hit, and test freeness.
offs = mask -> select(toList(0..7), v -> (mask >> v) % 2 == 1);
lineList = (F, mu, el, S) -> (
   -- returns list of (pencil, coefficient vector) with offsets el(v)
   {(0, {0_F, 0_F, 1_F})} | flatten apply(4, j -> apply(offs S#j, v -> (j + 1,
      if j == 0 then {0_F, 1_F, el v} else if j == 1 then {1_F, 0_F, el v}
      else if j == 2 then {1_F, 1_F, el v} else {1_F, mu, el v}))));
conc = (L) -> (  -- set of concurrent triples (indices)
   n := #L; select(subsets(n, 3), s -> det matrix apply(s, i -> (L#i)#1) == 0));
freeness = (F, L) -> (S := F[x,y,z]; f := product(L, l -> (l#1)#0*x + (l#1)#1*y + (l#1)#2*z);
   J := ideal(diff(x, f), diff(y, f), diff(z, f)); C := res J;
   (length C, sort apply(flatten degrees C_2, d -> d - (#L - 1))));
runHit = S -> (
   R0 := ZZ/2[t]; F8 := toField(R0/ideal(t^3+t+1)); tt := sub(t, F8);
   el8 := v -> sum(3, i -> if (v >> i) % 2 == 1 then tt^i else 0_F8);
   L8 := lineList(F8, tt, el8, S);
   C8 := conc L8;
   print("hit " | toString S | ": " | toString(#L8) | " lines, " | toString(#C8) | " concurrent triples; over F_8: " | toString freeness(F8, L8));
   -- generic mu: GF(2^11)
   G := GF(2^11, Variable => w); mu := w;
   -- unknown offsets s_0..s_{K-1}, one per affine line, in the order of L8 (minus H)
   K := #L8 - 1;
   R := G[s_0..s_(K-1)];
   lin := apply(K, i -> (j := (L8#(i+1))#0; {if j == 0 then 0_R else if j == 3 then mu_R else 1_R, if j == 1 then 0_R else 1_R} ));
   -- line i: lin#i#0 * x + lin#i#1 * y + s_i z  (pencil 1: y + s z; 2: x + s z; 3: x+y+s z; 4: x + mu y + s z)
   lin = apply(K, i -> (j := (L8#(i+1))#0; if j == 1 then {0_R, 1_R} else if j == 2 then {1_R, 0_R} else if j == 3 then {1_R, 1_R} else {1_R, mu_R}));
   eqs := flatten apply(select(C8, c -> not member(0, c)), c -> {det matrix apply(c, i -> lin#(i-1) | {s_(i-1)})});
   M := transpose matrix apply(eqs, e -> apply(K, i -> coefficient(s_i, e)));
   ker1 := gens ker sub(transpose M, G);
   print("  generic mu: kernel dimension " | toString(numcols ker1));
   -- translations: alpha (x -> x + alpha): s_i += lin#i#0 * alpha ; beta (y -> y + beta): s_i += lin#i#1 * beta
   T := matrix apply(K, i -> {sub(lin#i#0, G), sub(lin#i#1, G)});
   B := ker1 | T;
   -- pick a kernel column independent of the translations, combined randomly with the others
   v := ker1 * random(G^(numcols ker1), G^1);
   Lg := {(0, {0_G, 0_G, 1_G})} | apply(K, i -> ((L8#(i+1))#0, {sub(lin#i#0, G), sub(lin#i#1, G), v_(i,0)}));
   Cg := conc Lg;
   distinct := #unique apply(Lg, l -> (m := matrix{l#1}; (1/(first select(flatten entries m, e -> e != 0))) * m)) == #Lg;
   print("  generic member: lines distinct " | toString distinct | ", same concurrent triples as the hit: " | toString(set Cg === set C8) | " (" | toString(#Cg) | ")");
   print("  generic member freeness over GF(2^11): " | toString freeness(G, Lg));
);
runHit {0xfe, 0xc1, 0xcc, 0xf0};
runHit {0xfe, 0xc1, 0x3c, 0x3c};
