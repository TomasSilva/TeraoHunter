S = QQ[x,y,z];
load "thmF_example_lines.m2";
arrInfo = L -> (
  L = apply(L, l -> sub(l, S)); f := product L;
  Sy := syz matrix{{diff(x,f), diff(y,f), diff(z,f)}};
  degs := sort apply(entries transpose Sy, c -> first degree c#0);
  -- multiple points profile
  pts := {}; for i from 0 to #L-1 do for j from i+1 to #L-1 do (
     I := ideal(L#i, L#j); pts = append(pts, I));
  print("syzygy degrees of the Jacobian ideal: " | toString degs));
print "generic member (mu = 7/3):"; arrInfo gen;
print "limit member (mu0 = 2):"; arrInfo lim;
