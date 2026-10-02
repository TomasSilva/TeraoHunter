-- Exact non-free locus of the 15-line lattice (exps (1,7,7), 6 quintuple + 6 triple points)
-- on its full realisation space  Q5 = (d^2:d:2d-1), Q6 = (d(2d-1):d^2:2d-1).
R = QQ[d];
Q = {matrix{{1_R,0,0}}, matrix{{0_R,1,0}}, matrix{{0_R,0,1}}, matrix{{1_R,1,1}},
     matrix{{d^2, d, 2*d-1}}, matrix{{d*(2*d-1), d^2, 2*d-1}}};
cross = (u,v) -> matrix{{u_(0,1)*v_(0,2)-u_(0,2)*v_(0,1), u_(0,2)*v_(0,0)-u_(0,0)*v_(0,2), u_(0,0)*v_(0,1)-u_(0,1)*v_(0,0)}};
prs = subsets(6, 2);
L = apply(prs, p -> cross(Q#(p#0), Q#(p#1)));        -- line (i,j) through Q_i, Q_j
-- degeneracy polynomial B(d): collinear Q's and the 9 non-prescribed matchings becoming concurrent
coll = apply(subsets(6, 3), s -> det(Q#(s#0) || Q#(s#1) || Q#(s#2)));
Mpres = {{{0,1},{2,3},{4,5}}, {{0,2},{1,5},{3,4}}, {{0,5},{1,4},{2,3}}, {{0,3},{1,5},{2,4}}, {{0,4},{1,2},{3,5}}, {{0,4},{1,3},{2,5}}};
lineOf = p -> L#(position(prs, q -> q == sort p));
conc = m -> det(lineOf(m#0) || lineOf(m#1) || lineOf(m#2));
allM = unique flatten apply(subsets(6,2), s -> (rest := toList(set(0..5) - set s);
          apply(subsets(rest, 2), t -> sort {sort s, sort t, sort toList(set rest - set t)})));
others = select(allM, m -> not member(m, apply(Mpres, x -> sort apply(x, sort))));
print("prescribed concurrencies vanish identically: " | toString(all(Mpres, m -> conc m == 0)));
B = product(coll) * product(apply(others, conc));
print("degeneracy polynomial B factors: " | toString factor B);
-- derivations of degree k with theta(z) = 0 (z = line Q1Q2):  theta = f d_x + g d_y
k = 6;
S = R[x,y,z]; T = R[s,t];
mons = flatten entries basis(k, S);
cols = {};
rows = {};
condCols = for slot from 0 to 1 list for mu in mons list (
   flatten for i from 0 to #prs-1 list (
      p := prs#i;
      if p == {0,1} then continue;
      a := L#i;
      P := Q#(p#0); Qv := Q#(p#1);
      sub1 := map(T, S, {s*P_(0,0)+t*Qv_(0,0), s*P_(0,1)+t*Qv_(0,1), s*P_(0,2)+t*Qv_(0,2)});
      coef := if slot == 0 then a_(0,0) else a_(0,1);
      F := coef * sub1(mu);
      apply(k+1, j -> lift(coefficient(s^j*t^(k-j), F), R))
   ));
M = transpose matrix flatten condCols;      -- rows: conditions, columns: unknowns
print("matrix size: " | toString(numrows M) | " x " | toString(numcols M));
Kf = frac R;
print("generic rank (over QQ(d)): " | toString rank(sub(M, Kf)));
N = coker transpose M;
an = ann N;
print("annihilator of the torsion module (non-free locus): " | toString an);
g = an_0;
print("factored: " | toString factor g);
-- every non-free value must be a degeneracy value of the lattice
print("non-free locus contained in degeneracy locus: " | toString((B % radical an) == 0));
