-- Absorption check: adding a rich line l (|A meet l| = r - z + 1) lowers the Bourbaki length by z(l);
-- when it reaches 0 the enlarged arrangement is free with exponents (r, |A| - r).
syzdegs = HK -> (K := ring first first HK; S := K[x,y,z];
  f := product apply(HK, h -> h#0*x + h#1*y + h#2*z);
  Sy := syz matrix{{diff(x,f), diff(y,f), diff(z,f)}};
  sort apply(numcols Sy, i -> first degree(Sy_(0,i))));
K6 = toField(QQ[w]/(w^2+w+1)); use K6;
zt = 1 + w; mu6 = apply(6, k -> zt^k);
seed = {{0_K6,1_K6,0_K6},{0_K6,0_K6,1_K6}} | apply(6, k -> {1_K6, -(zt^k), 0_K6}) | apply(6, k -> {-(zt^k), 0_K6, 1_K6}) |
       apply(select(mu6, b -> b != 1 and b != zt), b -> {0_K6, 1_K6, -b});
print("seed (18 lines): syzygy degrees " | toString syzdegs seed);
print("seed + {x=0} (19 lines, rich line with z = 2): " | toString syzdegs(seed | {{1_K6,0_K6,0_K6}}) | "   predicted free (7,11)");
print("seed + {y=z} (19 lines, rich line with z = 1): " | toString syzdegs(seed | {{0_K6,1_K6,-1_K6}}) | "   predicted Bourbaki 1: not free");
print("seed + all 3 rich lines (21 lines): " | toString syzdegs(seed | {{1_K6,0_K6,0_K6},{0_K6,1_K6,-1_K6},{0_K6,1_K6,-zt}}) | "   predicted free (7,13) = A(6,1,3)");
K5 = toField(QQ[e]/(e^4+e^3+e^2+e+1)); use K5;
ceva5 = apply(5, k -> {1_K5, -(e^k), 0_K5}) | apply(5, k -> {0_K5, 1_K5, -(e^k)}) | apply(5, k -> {-(e^k), 0_K5, 1_K5});
nf = drop(ceva5, 1);
print("Ceva(5) minus a line (14): " | toString syzdegs nf);
print("  + {z=0} (rich, z = 1): " | toString syzdegs(nf | {{0_K5,0_K5,1_K5}}) | "   predicted free (6,8)");
print("  + {x=0} (rich, z = 0): " | toString syzdegs(nf | {{1_K5,0_K5,0_K5}}) | "   predicted Bourbaki stays 1");
