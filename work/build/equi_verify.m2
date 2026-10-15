load "../rigidity/specialise.m2";
load "../contradiction/richlines.m2";
K = toField(QQ[t]/(t^2-t+1)); use K;
o = 1_K;
A = {0*o, o, t-1, t}; B = {0*o, o, t, 1-t}; C = {0*o, o, t, t+1}; D = {0*o, o, t-1, t};
E17 = {{0*o,0*o,o}} | apply(A, a -> {0*o, -o, -a}) | apply(B, b -> {o, 0*o, -b}) | apply(C, c -> {o, -o, -c}) | apply(D, d -> {o, -t, -d});
analyseAt("E17 (equianharmonic four-pencil)", E17);
checkRich("E17", E17);
-- the harmonic grid G17 for comparison (same weak combinatorics)
G17 = {{0*o,0*o,o}} | apply({0,1,2,3}, a -> {0*o, -o, -a*o}) | apply({2,3,4,5}, b -> {o, 0*o, -b*o}) | apply({0,1,2,3}, c -> {o, -o, -c*o}) | apply({3,4,5,6}, d -> {o, o, -d*o});
checkRich("G17 (harmonic grid, free)", G17);
