-- The critically fixed maps behind the rational (harmonic) jump points of 4-point lines with K = 16.
-- Points x (k1), y (k2), x - y (k3), x - mu y (k4), i.e. z = x/y = 0, inf, 1, mu.
load "../rational/jump4.m2";
pats = {({1,1,7,7},-1), ({1,3,6,6},-1), ({1,4,4,7},null), ({1,5,5,5},-1),
        ({2,2,5,7},null), ({3,3,3,7},null), ({3,3,5,5},-1), ({3,4,4,5},null)};
T = QQ[X,Y];
for p in pats do (
  k := p#0;
  rr := ratRoots jumpPoly(k, 7);
  for m0 in rr do (
    M := sub(condMatrix(k, 7), {mu => m0});
    Mq := lift(M, QQ);
    kv := flatten entries gens kernel Mq;
    if #kv == 0 then continue;
    f := sum(8, i -> kv#i * X^(7-i)*Y^i); g := sum(8, i -> kv#(8+i) * X^(7-i)*Y^i);
    W := X*g - Y*f; Jc := det matrix{{diff(X,f), diff(Y,f)}, {diff(X,g), diff(Y,g)}};
    print("k = " | toString k | "  mu = " | toString m0 | "  gcd(f,g) degree " | toString(first degree gcd(f,g)));
    print("   fixed points  W = " | toString factor W);
    print("   critical pts Jac = " | toString factor Jc)));
