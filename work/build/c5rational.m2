load "../contradiction/richlines.m2";
q = x -> x * 1/1;
-- line  dq*x - dp*y - v*z = 0 for direction (dp,dq) and offset v
pen = (dp, dq, offs) -> apply(offs, v -> {q dq, q(-dp), q(-v)});
R0 = {{q 0, q 0, q 1}} | pen(1,0,{0,1,2,3}) | pen(0,1,{0,2,4}) | pen(1,1,{2,3,5}) | pen(4/3,1,{4}) | pen(2,1,{2,4,5,6,8});
R2 = {{q 0, q 0, q 1}} | pen(1,0,{0,1,2,3,4}) | pen(0,1,{0,1,2,4}) | pen(1,1,{4}) | pen(-1/2,1,{1}) | pen(-1,1,{-2,-1,0,1,2});
checkRich("R17a (H-vector 5,4,3,3,1; t = 4/3, 2)", R0);
checkRich("R17b (H-vector 5,5,4,1,1; t = -1/2, -1)", R2);
