load "linejumps.m2";
-- 15-line (1,7,7) family: the 15 lines joining Q1..Q6 (TC_ATTEMPT Theorem 2), lambda = l1/l0.
Q = {{1_R,0_R,0_R}, {0_R,1_R,0_R}, {0_R,0_R,1_R}, {1_R,1_R,1_R},
     {l1^2, l1*l0, 2*l1*l0-l0^2}, {l1*(2*l1-l0), l1^2, (2*l1-l0)*l0}};
H15 = apply(subsets(6,2), p -> prim crossv(Q#(p#0), Q#(p#1)));
out = analyseFamily("15-line (7,7)", H15);
print("SUMMARY (line, deg ThetaH, #interior, #remote missed): " | toString out);
