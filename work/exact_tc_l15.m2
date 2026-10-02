R = QQ[u_0,u_1,u_2,u_3];
cr = (a,b) -> matrix{{a_(0,1)*b_(0,2)-a_(0,2)*b_(0,1), a_(0,2)*b_(0,0)-a_(0,0)*b_(0,2), a_(0,0)*b_(0,1)-a_(0,1)*b_(0,0)}};
X0 = sub(matrix{{1,0,0}}, R);
X1 = sub(matrix{{0,1,0}}, R);
X2 = sub(matrix{{0,0,1}}, R);
X3 = sub(matrix{{1,1,1}}, R);
X4 = sub(matrix{{u_0,u_1,1}}, R);
X5 = sub(matrix{{u_2,u_3,1}}, R);
L0 = cr(X0, X1);
L1 = cr(X0, X2);
L2 = cr(X0, X3);
L3 = cr(X0, X4);
L4 = cr(X0, X5);
L5 = cr(X1, X2);
L6 = cr(X1, X3);
L7 = cr(X1, X4);
L8 = cr(X1, X5);
L9 = cr(X2, X3);
L10 = cr(X2, X4);
L11 = cr(X2, X5);
L12 = cr(X3, X4);
L13 = cr(X3, X5);
L14 = cr(X4, X5);
I = ideal(det(L0 || L9 || L14),det(L1 || L6 || L14),det(L2 || L7 || L11),det(L2 || L8 || L10),det(L3 || L5 || L13),det(L4 || L5 || L12));
bad = {det(X0 || X1 || X2),det(X0 || X1 || X3),det(X0 || X1 || X4),det(X0 || X1 || X5),det(X0 || X2 || X3),det(X0 || X2 || X4),det(X0 || X2 || X5),det(X0 || X3 || X4),det(X0 || X3 || X5),det(X0 || X4 || X5),det(X1 || X2 || X3),det(X1 || X2 || X4),det(X1 || X2 || X5),det(X1 || X3 || X4),det(X1 || X3 || X5),det(X1 || X4 || X5),det(X2 || X3 || X4),det(X2 || X3 || X5),det(X2 || X4 || X5),det(X3 || X4 || X5),det(L0 || L10 || L13),det(L0 || L11 || L12),det(L1 || L7 || L13),det(L1 || L8 || L12),det(L2 || L5 || L14),det(L3 || L6 || L11),det(L3 || L8 || L9),det(L4 || L6 || L10),det(L4 || L7 || L9)};
bad = select(bad, f -> f != 0);
if any(bad, f -> f == 0) then print "WARNING zero degeneracy factor";
cs = if I == 0 then {ideal(0_R)} else decompose I;
good = select(cs, c -> all(bad, f -> f % c != 0));
print("components: " | toString(#cs) | ", nondegenerate: " | toString(#good));
S = R[x,y,z]; T = R[s,t]; kdeg = 6;
mons = flatten entries basis(kdeg, S);
lineData = {{L1,X0,X2},{L2,X0,X3},{L3,X0,X4},{L4,X0,X5},{L5,X1,X2},{L6,X1,X3},{L7,X1,X4},{L8,X1,X5},{L9,X2,X3},{L10,X2,X4},{L11,X2,X5},{L12,X3,X4},{L13,X3,X5},{L14,X4,X5}};
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

