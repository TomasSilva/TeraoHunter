R = ZZ/2[x_1..x_6];
K = {x_6,x_5,x_4,x_3,x_2,x_1,x_3+x_4+x_5,x_2+x_4+x_6,x_1+x_4+x_5,x_1+x_3+x_5,x_1+x_3+x_4+x_6,x_1+x_2+x_3+x_6,x_1+x_2+x_3+x_4};
QK = product K;
Rb = R/ideal(QK);
v = j -> matrix apply(6, i -> {(sub(R_i, Rb))^(2^j)});
M = v(0) | v(1) | v(2) | v(3) | v(4);
t0 = cpuTime();
r = v(5) % M;
print("MEMBERSHIP theta_F5 in span(E..theta_F4) mod Q_K: " | toString(r == 0) | "  (" | toString(cpuTime() - t0) | "s)");
