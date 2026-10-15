load "richlines.m2";
-- 1. the 18-line seed  yz(x^6-y^6)(z^6-x^6) prod_{b in mu6 - {1,zt}} (y - b z)   (chi-bar (t-8)(t-9), mdr 7)
K6 = toField(QQ[w]/(w^2+w+1)); use K6;
zt = 1 + w;  -- primitive 6th root of unity
mu6 = apply(6, k -> zt^k);
seed = {{0_K6,1_K6,0_K6},{0_K6,0_K6,1_K6}} | apply(6, k -> {1_K6, -(zt^k), 0_K6}) | apply(6, k -> {-(zt^k), 0_K6, 1_K6}) |
       apply(select(mu6, b -> b != 1 and b != zt), b -> {0_K6, 1_K6, -b});
checkRich("18-line seed", seed);
-- 2. Ceva(5) = A(5,5,3): free (6,8); its rich lines must be leaves of theta_1
K5 = toField(QQ[e]/(e^4+e^3+e^2+e+1)); use K5;
ceva5 = apply(5, k -> {1_K5, -(e^k), 0_K5}) | apply(5, k -> {0_K5, 1_K5, -(e^k)}) | apply(5, k -> {-(e^k), 0_K5, 1_K5});
checkRich("Ceva(5), free (6,8)", ceva5);
-- 3. Ceva(5) minus one line: nearly free (6,9)
checkRich("Ceva(5) minus a line (nearly free)", drop(ceva5, 1));
-- 4. Ceva(5) minus two lines from different pencils
checkRich("Ceva(5) minus two lines", drop(drop(ceva5, 1), {5,5}));
