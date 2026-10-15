load "jump5.m2";
P = QQ[A, B, C];
for k in {{5,5,4,1,1},{5,5,3,2,1},{5,5,2,2,2},{5,4,4,2,1},{5,4,3,3,1},{5,4,3,2,2},{5,3,3,3,2},{4,4,4,3,1},{4,4,4,2,2},{4,4,3,3,2},{4,3,3,3,3}} do (
   F := jumpCurve k;
   G := product apply(select(toList factor F, p -> first degree(p#0) > 1), p -> p#0);
   d := first degree G;
   Gh := sum(terms G, t -> (e := first exponents t; sub(leadCoefficient t, QQ) * A^(e#0) * B^(e#1) * C^(d - e#0 - e#1)));
   I := ideal Gh; Js := I + ideal jacobian I;
   sdeg := degree Js;
   print("k " | toString k | ": degree " | toString d | ", arithmetic genus " | toString((d-1)*(d-2)//2) | ", total Tjurina number " | toString sdeg | ", irreducible over QQ: " | toString(#(toList factor Gh) == 1 or #select(toList factor Gh, p -> first degree(p#0) > 0) == 1)));
