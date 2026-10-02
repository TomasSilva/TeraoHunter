R = QQ[u_0,u_1,u_2,u_3,u_4,u_5];
cr = (a,b) -> matrix{{a_(0,1)*b_(0,2)-a_(0,2)*b_(0,1), a_(0,2)*b_(0,0)-a_(0,0)*b_(0,2), a_(0,0)*b_(0,1)-a_(0,1)*b_(0,0)}};
X0 = sub(matrix{{1,0,0}}, R);
X1 = sub(matrix{{0,1,0}}, R);
X2 = sub(matrix{{u_0,1,0}}, R);
X3 = sub(matrix{{u_1,1,0}}, R);
X4 = sub(matrix{{u_2,1,0}}, R);
X5 = sub(matrix{{0,0,1}}, R);
X6 = sub(matrix{{u_3,0,1}}, R);
X7 = sub(matrix{{1,1,1}}, R);
X8 = sub(matrix{{u_4,u_5,1}}, R);
L0 = cr(X0, X1);
L1 = cr(X0, X5);
L2 = cr(X0, X7);
L3 = cr(X0, X8);
L4 = cr(X1, X5);
L5 = cr(X1, X6);
L6 = cr(X1, X7);
L7 = cr(X1, X8);
L8 = cr(X2, X5);
L9 = cr(X2, X6);
L10 = cr(X2, X7);
L11 = cr(X2, X8);
L12 = cr(X3, X5);
L13 = cr(X3, X6);
L14 = cr(X3, X7);
L15 = cr(X4, X5);
L16 = cr(X4, X6);
L17 = cr(X4, X7);
L18 = cr(X4, X8);
L19 = cr(X5, X7);
L20 = cr(X5, X8);
L21 = cr(X6, X7);
L22 = cr(X6, X8);
I = ideal(det(X0 || X1 || X2),det(X0 || X1 || X3),det(X0 || X1 || X4),det(X0 || X5 || X6),det(X3 || X7 || X8),det(L0 || L19 || L22),det(L0 || L20 || L21),det(L1 || L10 || L18),det(L1 || L11 || L17),det(L2 || L4 || L18),det(L2 || L5 || L15),det(L2 || L7 || L13),det(L2 || L8 || L22),det(L2 || L16 || L20),det(L3 || L4 || L16),det(L3 || L5 || L17),det(L3 || L6 || L12),det(L3 || L9 || L19),det(L3 || L15 || L21),det(L4 || L11 || L21),det(L5 || L10 || L20),det(L6 || L11 || L13),det(L6 || L11 || L15),det(L7 || L10 || L12),det(L7 || L10 || L16),det(L8 || L13 || L17),det(L8 || L14 || L16),det(L9 || L12 || L18),det(L9 || L14 || L15));
bad = {det(X0 || X1 || X5),det(X0 || X1 || X6),det(X0 || X1 || X7),det(X0 || X1 || X8),det(X0 || X2 || X5),det(X0 || X2 || X6),det(X0 || X2 || X7),det(X0 || X2 || X8),det(X0 || X3 || X5),det(X0 || X3 || X6),det(X0 || X3 || X7),det(X0 || X3 || X8),det(X0 || X4 || X5),det(X0 || X4 || X6),det(X0 || X4 || X7),det(X0 || X4 || X8),det(X0 || X5 || X7),det(X0 || X5 || X8),det(X0 || X6 || X7),det(X0 || X6 || X8),det(X0 || X7 || X8),det(X1 || X2 || X5),det(X1 || X2 || X6),det(X1 || X2 || X7),det(X1 || X2 || X8),det(X1 || X3 || X5),det(X1 || X3 || X6),det(X1 || X3 || X7),det(X1 || X3 || X8),det(X1 || X4 || X5),det(X1 || X4 || X6),det(X1 || X4 || X7),det(X1 || X4 || X8),det(X1 || X5 || X6),det(X1 || X5 || X7),det(X1 || X5 || X8),det(X1 || X6 || X7),det(X1 || X6 || X8),det(X1 || X7 || X8),det(X2 || X3 || X5),det(X2 || X3 || X6),det(X2 || X3 || X7),det(X2 || X3 || X8),det(X2 || X4 || X5),det(X2 || X4 || X6),det(X2 || X4 || X7),det(X2 || X4 || X8),det(X2 || X5 || X6),det(X2 || X5 || X7),det(X2 || X5 || X8),det(X2 || X6 || X7),det(X2 || X6 || X8),det(X2 || X7 || X8),det(X3 || X4 || X5),det(X3 || X4 || X6),det(X3 || X4 || X7),det(X3 || X4 || X8),det(X3 || X5 || X6),det(X3 || X5 || X7),det(X3 || X5 || X8),det(X3 || X6 || X7),det(X3 || X6 || X8),det(X4 || X5 || X6),det(X4 || X5 || X7),det(X4 || X5 || X8),det(X4 || X6 || X7),det(X4 || X6 || X8),det(X4 || X7 || X8),det(X5 || X6 || X7),det(X5 || X6 || X8),det(X5 || X7 || X8),det(X6 || X7 || X8),det(L1 || L6 || L11),det(L1 || L6 || L18),det(L1 || L7 || L10),det(L1 || L7 || L17),det(L2 || L4 || L9),det(L2 || L4 || L11),det(L2 || L4 || L13),det(L2 || L4 || L16),det(L2 || L4 || L22),det(L2 || L5 || L8),det(L2 || L5 || L11),det(L2 || L5 || L12),det(L2 || L5 || L18),det(L2 || L5 || L20),det(L2 || L7 || L8),det(L2 || L7 || L9),det(L2 || L7 || L12),det(L2 || L7 || L15),det(L2 || L7 || L16),det(L2 || L8 || L13),det(L2 || L8 || L16),det(L2 || L8 || L18),det(L2 || L9 || L12),det(L2 || L9 || L15),det(L2 || L9 || L18),det(L2 || L9 || L20),det(L2 || L11 || L12),det(L2 || L11 || L13),det(L2 || L11 || L15),det(L2 || L11 || L16),det(L2 || L12 || L16),det(L2 || L12 || L18),det(L2 || L12 || L22),det(L2 || L13 || L15),det(L2 || L13 || L18),det(L2 || L13 || L20),det(L2 || L15 || L22),det(L3 || L4 || L9),det(L3 || L4 || L10),det(L3 || L4 || L13),det(L3 || L4 || L17),det(L3 || L4 || L21),det(L3 || L5 || L8),det(L3 || L5 || L10),det(L3 || L5 || L12),det(L3 || L5 || L15),det(L3 || L5 || L19),det(L3 || L6 || L8),det(L3 || L6 || L9),det(L3 || L6 || L13),det(L3 || L6 || L15),det(L3 || L6 || L16),det(L3 || L8 || L13),det(L3 || L8 || L16),det(L3 || L8 || L17),det(L3 || L8 || L21),det(L3 || L9 || L12),det(L3 || L9 || L15),det(L3 || L9 || L17),det(L3 || L10 || L12),det(L3 || L10 || L13),det(L3 || L10 || L15),det(L3 || L10 || L16),det(L3 || L12 || L16),det(L3 || L12 || L17),det(L3 || L12 || L21),det(L3 || L13 || L15),det(L3 || L13 || L17),det(L3 || L13 || L19),det(L3 || L16 || L19),det(L4 || L9 || L14),det(L4 || L9 || L17),det(L4 || L9 || L18),det(L4 || L10 || L13),det(L4 || L10 || L16),det(L4 || L10 || L18),det(L4 || L10 || L22),det(L4 || L11 || L13),det(L4 || L11 || L16),det(L4 || L11 || L17),det(L4 || L13 || L17),det(L4 || L13 || L18),det(L4 || L14 || L16),det(L4 || L17 || L22),det(L4 || L18 || L21),det(L5 || L8 || L14),det(L5 || L8 || L17),det(L5 || L8 || L18),det(L5 || L10 || L12),det(L5 || L10 || L15),det(L5 || L10 || L18),det(L5 || L11 || L12),det(L5 || L11 || L15),det(L5 || L11 || L17),det(L5 || L11 || L19),det(L5 || L12 || L17),det(L5 || L12 || L18),det(L5 || L14 || L15),det(L5 || L17 || L20),det(L5 || L18 || L19),det(L6 || L8 || L13),det(L6 || L8 || L16),det(L6 || L8 || L18),det(L6 || L8 || L22),det(L6 || L9 || L12),det(L6 || L9 || L15),det(L6 || L9 || L18),det(L6 || L9 || L20),det(L6 || L11 || L12),det(L6 || L11 || L16),det(L6 || L12 || L16),det(L6 || L12 || L18),det(L6 || L12 || L22),det(L6 || L13 || L18),det(L6 || L13 || L20),det(L6 || L15 || L22),det(L6 || L16 || L20),det(L7 || L8 || L13),det(L7 || L8 || L16),det(L7 || L8 || L17),det(L7 || L8 || L21),det(L7 || L9 || L12),det(L7 || L9 || L15),det(L7 || L9 || L17),det(L7 || L9 || L19),det(L7 || L10 || L13),det(L7 || L10 || L15),det(L7 || L12 || L17),det(L7 || L12 || L21),det(L7 || L13 || L15),det(L7 || L13 || L17),det(L7 || L13 || L19),det(L7 || L15 || L21),det(L7 || L16 || L19),det(L8 || L13 || L18),det(L8 || L17 || L22),det(L8 || L18 || L21),det(L9 || L12 || L17),det(L9 || L17 || L20),det(L9 || L18 || L19),det(L10 || L12 || L18),det(L10 || L12 || L22),det(L10 || L13 || L15),det(L10 || L13 || L18),det(L10 || L13 || L20),det(L10 || L15 || L22),det(L10 || L16 || L20),det(L11 || L12 || L16),det(L11 || L12 || L17),det(L11 || L12 || L21),det(L11 || L13 || L17),det(L11 || L13 || L19),det(L11 || L15 || L21),det(L11 || L16 || L19),det(L12 || L17 || L22),det(L12 || L18 || L21),det(L13 || L17 || L20),det(L13 || L18 || L19)};
bad = select(bad, f -> f != 0);
if any(bad, f -> f == 0) then print "WARNING zero degeneracy factor";
-- remove degenerate components first: saturate by each degeneracy factor in turn (cheap factors first)
G = I;
bad = sort(bad, f -> (first degree f, size f));
for f in bad do (if G != ideal(1_R) then G = saturate(G, f));
if G == ideal(1_R) then (print "components: 0, nondegenerate: 0  (NOT REALIZABLE in characteristic 0)"; exit 0);
cs = if G == 0 then {ideal(0_R)} else decompose G;
good = select(cs, c -> all(bad, f -> f % c != 0));
print("components: " | toString(#cs) | ", nondegenerate: " | toString(#good) | ", dims " | toString apply(good, dim));
S = R[x,y,z]; T = R[s,t]; kdeg = 10;
mons = flatten entries basis(kdeg, S);
lineData = {{L1,X0,X5},{L2,X0,X7},{L3,X0,X8},{L4,X1,X5},{L5,X1,X6},{L6,X1,X7},{L7,X1,X8},{L8,X2,X5},{L9,X2,X6},{L10,X2,X7},{L11,X2,X8},{L12,X3,X5},{L13,X3,X6},{L14,X3,X7},{L15,X4,X5},{L16,X4,X6},{L17,X4,X7},{L18,X4,X8},{L19,X5,X7},{L20,X5,X8},{L21,X6,X7},{L22,X6,X8}};
condCols = for slot from 0 to 1 list for mu in mons list (
   flatten for D in lineData list (
      al := D#0; P1 := D#1; P2 := D#2;
      sb := map(T, S, {s*P1_(0,0)+t*P2_(0,0), s*P1_(0,1)+t*P2_(0,1), s*P1_(0,2)+t*P2_(0,2)});
      F := al_(0,slot) * sb(mu);
      apply(kdeg+1, j -> lift(coefficient(s^j*t^(kdeg-j), F), R))));
M = transpose matrix flatten condCols;
print("derivation matrix " | toString(numrows M) | " x " | toString(numcols M) | " in degree " | toString kdeg);
Bprod = product bad;
for c in good do (
   A := R/c;
   if dim c == 0 then (
      F := toField A;
      r := rank sub(M, F);
      print("COMPONENT dim 0, degree " | toString(degree c) | ": rank " | toString r | " / " | toString(numcols M) | " -> " | (if r == numcols M then "FREE" else "NOT FREE"));
   ) else (
      r := rank sub(M, frac A);
      an := ann coker transpose sub(M, A);
      ok := (sub(Bprod, A) % radical an) == 0;
      print("COMPONENT dim " | toString(dim c) | ": generic rank " | toString r | " / " | toString(numcols M) | (if r == numcols M then " (generic member FREE)" else " (generic member NOT free)") | "; non-free locus inside degeneracy locus: " | toString ok);
   ));

