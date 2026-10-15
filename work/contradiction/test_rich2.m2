load "richlines.m2";
Ki = toField(QQ[jj]/(jj^2+1)); use Ki;
-- A(4,1,3) = xyz(x^4-y^4)(y^4-z^4)(z^4-x^4): free (5,9); prediction: no rich line at all (3*5-15 = 0)
fourth = {1_Ki, jj, -1_Ki, -jj};
a413 = {{1_Ki,0_Ki,0_Ki},{0_Ki,1_Ki,0_Ki},{0_Ki,0_Ki,1_Ki}} | apply(fourth, u -> {1_Ki, -u, 0_Ki}) |
       apply(fourth, u -> {0_Ki, 1_Ki, -u}) | apply(fourth, u -> {-u, 0_Ki, 1_Ki});
checkRich("A(4,1,3), free (5,9)", a413);
checkRich("A(4,1,3) minus x=0", drop(a413, 1));
checkRich("A(4,1,3) minus two pencil lines", drop(drop(a413, {3,3}), {7,7}));
-- B3 (9 lines, free (3,5)) and a deletion; and the A(3,3,3)+lines family member
b3 = {{1,0,0},{0,1,0},{0,0,1},{1,1,0},{1,-1,0},{1,0,1},{1,0,-1},{0,1,1},{0,1,-1}};
b3 = apply(b3, v -> apply(v, c -> c*1/1));
checkRich("B3, free (3,5)", b3);
checkRich("B3 minus a line", drop(b3, 1));
