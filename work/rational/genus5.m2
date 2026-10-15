-- geometric genus of the 5-point jump curves (homogenised in P^2)
load "jump5.m2";
needsPackage "AdjointIdeal";
k = value(scriptCommandLine#1);
F = jumpCurve k;
G = product apply(select(toList factor F, p -> first degree(p#0) > 1), p -> p#0);
P = QQ[A, B, C];
d = first degree G;
Gh = sum(terms G, t -> (e := first exponents t; leadCoefficient t * A^(e#0) * B^(e#1) * C^(d - e#0 - e#1)));
I = ideal Gh;
print("k " | toString k | ": degree " | toString d | ", singular points (scheme degree) " | toString degree(I + ideal jacobian I));
g = elapsedTime geometricGenus(I);
print("k " | toString k | ": geometric genus " | toString g);
