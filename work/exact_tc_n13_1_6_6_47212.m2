R = QQ[u_0,u_1,u_2];
cr = (a,b) -> matrix{{a_(0,1)*b_(0,2)-a_(0,2)*b_(0,1), a_(0,2)*b_(0,0)-a_(0,0)*b_(0,2), a_(0,0)*b_(0,1)-a_(0,1)*b_(0,0)}};
X0 = sub(matrix{{1,0,0}}, R);
X1 = sub(matrix{{0,1,0}}, R);
X2 = sub(matrix{{u_0,1,0}}, R);
X3 = sub(matrix{{0,0,1}}, R);
X4 = sub(matrix{{1,1,1}}, R);
X5 = sub(matrix{{u_1,u_2,1}}, R);
L0 = cr(X0, X1);
L1 = cr(X0, X3);
L2 = cr(X0, X4);
L3 = cr(X0, X5);
L4 = cr(X1, X3);
L5 = cr(X1, X4);
L6 = cr(X1, X5);
L7 = cr(X2, X3);
L8 = cr(X2, X4);
L9 = cr(X2, X5);
L10 = cr(X3, X4);
L11 = cr(X3, X5);
L12 = cr(X4, X5);
I = ideal(det(X0 || X1 || X2),det(L1 || L6 || L8),det(L2 || L4 || L9),det(L3 || L5 || L7));
bad = {det(X0 || X1 || X3),det(X0 || X1 || X4),det(X0 || X1 || X5),det(X0 || X2 || X3),det(X0 || X2 || X4),det(X0 || X2 || X5),det(X0 || X3 || X4),det(X0 || X3 || X5),det(X0 || X4 || X5),det(X1 || X2 || X3),det(X1 || X2 || X4),det(X1 || X2 || X5),det(X1 || X3 || X4),det(X1 || X3 || X5),det(X1 || X4 || X5),det(X2 || X3 || X4),det(X2 || X3 || X5),det(X2 || X4 || X5),det(X3 || X4 || X5),det(L1 || L5 || L9),det(L2 || L6 || L7),det(L3 || L4 || L8)};
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
S = R[x,y,z]; T = R[s,t]; kdeg = 5;
mons = flatten entries basis(kdeg, S);
lineData = {{L1,X0,X3},{L2,X0,X4},{L3,X0,X5},{L4,X1,X3},{L5,X1,X4},{L6,X1,X5},{L7,X2,X3},{L8,X2,X4},{L9,X2,X5},{L10,X3,X4},{L11,X3,X5},{L12,X4,X5}};
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

