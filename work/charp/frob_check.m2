-- Frobenius bound: an arrangement defined over F_q with p not dividing n has a Jacobian syzygy of degree q
mdrOf = Q -> (R := ring Q; jac := matrix{{diff(R_0,Q),diff(R_1,Q),diff(R_2,Q)}}; S := mingens image syz jac;
   sort apply(numcols S, i -> first (degrees source S)#i - (first degree Q - 1)));
-- (1) the seed over F_7 (zeta = 3 has order 6 in F_7)
kk = ZZ/7; R = kk[x,y,z]; w = 3_kk;
L = {y, z} | apply(6, i -> x - w^i*y) | apply(6, j -> z - w^j*x) | apply({2,3,4,5}, k -> y - w^k*z);
<< "seed over F_7: n=" << #L << " syz degrees " << mdrOf product L << endl;
-- (2) random 11 lines over F_4 = GF(4): mdr should be <= 4
kk = GF(4); R = kk[x,y,z];
L = apply(11, i -> (random kk)*x + (random kk)*y + (random kk)*z);
<< "random 11 lines over F_4 (distinct: " << (#unique L == 11) << "): syz degrees " << mdrOf product L << endl;
-- (3) random 13 lines over F_5
kk = ZZ/5; R = kk[x,y,z];
L = apply(13, i -> (random kk)*x + (random kk)*y + (random kk)*z);
<< "random 13 lines over F_5: syz degrees " << mdrOf product L << endl;
-- (4) the explicit Frobenius syzygy for (3)
Q = product L; sigma = sum apply(L, l -> l^4); n = 13;
th = {x^5 - (1/n)*sigma*x, y^5 - (1/n)*sigma*y, z^5 - (1/n)*sigma*z};
<< "explicit theta' kills Q: " << (th#0*diff(x,Q) + th#1*diff(y,Q) + th#2*diff(z,Q) == 0) << endl;
-- (5) full PG(2,4): exponents (4,16)
kk = GF(4); R = kk[x,y,z]; a = kk_0;
pts = {}; els = {0_kk, 1_kk, a, a+1};
L = unique flatten {apply(els, b -> apply(els, c -> x + b*y + c*z)), apply(els, c -> y + c*z), {z}};
<< "PG(2,4): n=" << #L << " syz degrees " << mdrOf product L << endl;
