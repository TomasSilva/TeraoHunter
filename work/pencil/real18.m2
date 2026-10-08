-- realization space of the lattice of the pasted example: three pencils through the coordinate points,
-- A: x = a_i y (i=0..5), B: z = c_j x (j=0..5), C: y = b_k z (k=2..5); triple point (i,j,k) iff k = -(i+j) mod 6, k in {2,3,4,5}
R = QQ[a_1..a_5, c_1..c_5, b_2..b_5];
av = i -> if i == 0 then 1_R else a_i; cv = j -> if j == 0 then 1_R else c_j;
eqs = flatten for i from 0 to 5 list for j from 0 to 5 list ( k := (-(i+j)) % 6; if member(k, {2,3,4,5}) then av(i)*cv(j)*b_k - 1 else continue );
<< "#equations = " << #eqs << endl;
I = ideal eqs;
<< "dim = " << dim I << "  degree = " << degree I << endl;
P = primaryDecomposition I;
<< "#components = " << #P << endl;
scan(P, q -> << "  component: dim " << dim q << " degree " << degree q << " : " << toString mingens q << endl);
