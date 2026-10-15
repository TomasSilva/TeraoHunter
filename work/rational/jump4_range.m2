load "jump4.m2";
Kmin = value(scriptCommandLine#1); Kmax = value(scriptCommandLine#2);
for K from Kmin to Kmax do
 for k1 from 1 to K do for k2 from k1 to K do for k3 from k2 to K do (
   k4 := K - k1 - k2 - k3;
   if k4 >= k3 and 2 * k4 <= K then (
      k := {k1, k2, k3, k4};
      eg := genericE1 k;
      g := jumpPoly(k, eg - 1);
      fac := if g == 0 then "ALWAYS" else toString factor g;
      rr := if g == 0 then {} else select(ratRoots g, r -> r != 0 and r != 1);
      print("K " | toString K | " k " | toString k | " e1gen " | toString eg | " ratjump " | toString rr | " | " | fac);
   ));
