-- the degree-8 Jacobian syzygy of the generic member of family 1 (over F_2(mu))
A = frac(ZZ/2[mu]);
P = A[x,y,z];
V = {0_A, 1_A, mu^2+mu+1, mu+1, mu^2+1, mu^2, mu};
lins = {z} | apply(V, a -> y + a*z) | apply({0_A, mu, mu^2}, b -> x + b*z) | apply({1_A, mu^2+mu+1, mu^2, mu}, c -> x + y + c*z) |
       apply({mu, mu^2, mu^3+mu^2+mu, mu^3}, d -> x + mu*y + d*z);
f = product lins;
Sy = gens ker matrix{{diff(x, f), diff(y, f), diff(z, f)}};
i8 = first select(numcols Sy, i -> first degree (Sy_{i})_(0,0) == 8 or (Sy_(0,i) == 0 and first degree (Sy_{i})_(1,0) == 8));
th = flatten entries Sy_{i8};
print("degree-8 syzygy (a, b, c) with a f_x + b f_y + c f_z = 0:");
scan(th, e -> print("   " | toString factor e));
-- compare with the additive derivation of V = span(1, mu, mu^2): P_V(X,Z) = prod_{v in V+0} (X + v Z)
Vfull = {0_A, 1_A, mu, mu+1, mu^2, mu^2+1, mu^2+mu, mu^2+mu+1};
PV = (X, Z) -> product(Vfull, v -> X + v*Z);
print("P_V(x,z) = " | toString PV(x, z));
