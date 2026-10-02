
R = ZZ/2[x_1..x_6];
L = {x_5+x_6,x_4+x_6,x_4+x_5,x_4+x_5+x_6,x_3+x_6,x_3+x_5,x_3+x_5+x_6,x_3+x_4,x_3+x_4+x_6,x_3+x_4+x_5+x_6,x_2+x_6,x_2+x_5,x_2+x_5+x_6,x_2+x_4,x_2+x_4+x_5,x_2+x_4+x_5+x_6,x_2+x_3,x_2+x_3+x_6,x_2+x_3+x_5,x_2+x_3+x_5+x_6,x_2+x_3+x_4,x_2+x_3+x_4+x_6,x_2+x_3+x_4+x_5,x_2+x_3+x_4+x_5+x_6,x_1+x_6,x_1+x_5,x_1+x_5+x_6,x_1+x_4,x_1+x_4+x_6,x_1+x_4+x_5+x_6,x_1+x_3,x_1+x_3+x_6,x_1+x_3+x_5+x_6,x_1+x_3+x_4,x_1+x_3+x_4+x_5,x_1+x_3+x_4+x_5+x_6,x_1+x_2,x_1+x_2+x_6,x_1+x_2+x_5,x_1+x_2+x_5+x_6,x_1+x_2+x_4,x_1+x_2+x_4+x_6,x_1+x_2+x_4+x_5,x_1+x_2+x_4+x_5+x_6,x_1+x_2+x_3,x_1+x_2+x_3+x_5,x_1+x_2+x_3+x_5+x_6,x_1+x_2+x_3+x_4+x_6,x_1+x_2+x_3+x_4+x_5,x_1+x_2+x_3+x_4+x_5+x_6};
n = #L;
M = matrix apply(n, H -> join(apply(6, i -> coefficient(R_i, L#H)), apply(n, j -> if j == H then L#H else 0_R)));
Phi = map(R^n, , M);
t0 = cpuTime();
S = syz(Phi, DegreeLimit => 19);
print("SYZ done " | toString(cpuTime() - t0) | "s, generators: " | toString numcols S);
G = submatrix(S, {0..5}, );
G = mingens image G;
print("GENDEGS " | toString sort flatten degrees source G);
Q = product L;
fro = apply(5, j -> transpose matrix{apply(6, i -> R_i^(2^j))});
cand = select(numcols G, c -> (degrees source G)#c == {19});
found = false;
for c in cand do (
  T = fold((a, b) -> a | b, fro) | G_{c};
  dt = det T;
  if dt != 0 then (
     print("DET degree " | toString first degree dt | ", equals Q_A: " | toString(dt == Q));
     found = true;
     break;
  );
);
if not found then print("NO degree-19 element with nonzero determinant");
