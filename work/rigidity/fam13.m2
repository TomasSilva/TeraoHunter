load "linejumps.m2";
-- Abe-Cuntz-Kawanoue-Nozawa 13-line family (arXiv:1411.3351, Sec. 5), rescaled x -> x/sqrt3 so it is over QQ.
-- lambda = l1/l0.  Lines as coefficient vectors (x, y, z).
H13 = {
 {-l0, -l0, l1+l0}, {0_R, 2*l0, l1+l0}, {l0, -l0, l1+l0},
 {l0, -l0, l1-2*l0}, {-l0, -l0, l1-2*l0}, {0_R, 2*l0, l1-2*l0},
 {0_R, 2*l0, -2*l1+l0}, {l0, -l0, -2*l1+l0}, {-l0, -l0, -2*l1+l0},
 {(l0-l1)*l0, (l1+l0)*l0, -l1^2+l1*l0-l0^2},
 {l1*l0, (l1-2*l0)*l0, -l1^2+l1*l0-l0^2},
 {-l0^2, (l0-2*l1)*l0, -l1^2+l1*l0-l0^2},
 {0_R, 0_R, 1_R}};
H13 = apply(H13, prim);
out = analyseFamily("ACKN 13", H13);
print("SUMMARY (line, deg ThetaH, #interior, #remote missed): " | toString out);
