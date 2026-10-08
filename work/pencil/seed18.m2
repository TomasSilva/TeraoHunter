-- the pasted 18-line example: yz (x^6-y^6)(z^6-x^6) prod_{b in mu6 \ {1,w6}} (y - b z), over ZZ/p with w6 of order 6
p = 10009; kk = ZZ/p;
w6 = (for a from 2 to p-1 do (if (a^6 % p == 1) and (a^3 % p != 1) and (a^2 % p != 1) then break a));
<< "p=" << p << " w6=" << w6 << " w6^2-w6+1 mod p = " << (w6^2 - w6 + 1) % p << endl;
R = kk[x,y,z]; z6 = sub(w6, kk);
uniqIdeals = L -> ( out := {}; scan(L, p -> if not any(out, q -> q == p) then out = append(out, p)); out );
analyse = (name, lines) -> (
  Q := product lines; n := #lines;
  jac := matrix{{diff(x,Q), diff(y,Q), diff(z,Q)}};
  S := mingens image syz jac;
  degs := apply(numcols S, i -> first (degrees source S)#i - (n-1));
  r := min degs; i := position(degs, d -> d == r);
  th := S_{i}; a := th_(0,0); b := th_(1,0); c := th_(2,0);
  w1 := b*z - c*y; w2 := c*x - a*z; w3 := a*y - b*x;
  g := gcd(gcd(w1,w2),w3);
  W := matrix{{w1,w2,w3}}; J := det jacobian W;
  Zw := ideal(w1,w2,w3);
  pts := uniqIdeals apply(subsets(lines,2), p -> radical ideal p);
  ms := apply(pts, p -> number(lines, l -> isSubset(ideal l, p)));
  b2 := sum apply(ms, m -> m-1);
  chi := t -> t^2 - (n-1)*t + b2 - (n-1);
  multI := intersect pts;
  off := saturate(Zw, Q);
  Qx := diff(x,Q);
  hs := for k from 0 to numcols S - 1 list (
     e := S_{k}; a2 := e_(0,0); b2' := e_(1,0); c2 := e_(2,0);
     cr1 := b*c2 - c*b2'; cr1 // Qx );
  Inu := ideal hs; Inus := saturate Inu;
  << "== " << name << ": n=" << n << " syz degs=" << degs << " r=" << r << " j=3r-n=" << 3*r-n << " gcd(omega)=" << g
    << " | profile=" << tally ms << " b2=" << b2 << " chibar(r)=" << chi r << endl;
  << "   points per line: " << tally apply(lines, l -> number(pts, p -> isSubset(ideal l, p))) << endl;
  << "   Z(omega)=" << (if dim Zw <= 1 then degree Zw else "curve") << " (r^2+r+1=" << r^2+r+1 << ")"
    << "  off-line=" << (if dim off == 1 then degree off else 0) << endl;
  << "   Bourbaki ideal: dim R/I=" << dim Inus << " degree=" << degree Inus << endl;
  << "   Bourbaki scheme components: " << apply(primaryDecomposition Inus, q -> (degree q, radical q)) << endl;
  << "   off-line singular points: " << (if dim off == 1 then primaryDecomposition off else "none") << endl;
  << "   J/Q = " << (if J % Q == 0 then factor(J // Q) else "not divisible") << endl;
  offB := saturate(Inus, Q); << "   off-line part of Bourbaki = off-line singular scheme of omega? " << (offB == off) << endl;
  if J % Q == 0 then ( Jp := J // Q; sj := ideal(Jp, diff(x,Jp), diff(y,Jp), diff(z,Jp));
     << "   Sing(J/Q): " << apply(primaryDecomposition saturate sj, q -> (degree q, radical q)) << endl;
     << "   Bourbaki scheme inside Sing(J/Q)? " << isSubset(saturate sj, Inus) << endl; );
  -- Ziegler exponents: restriction to each line (multiarrangement); e1 = min degree of a syzygy of the restricted Jacobian... use the standard trick:
  -- e1(H) = mdr of the multiarrangement; compute via the restriction of D0(A) in the ambient: dim of D0(A)_k restricted. Simpler: Ziegler: e1+e2 = n-1 and e1 = min k with a nonzero derivation of the restricted multiarrangement.
  << endl;
);
mu6 = apply(6, i -> z6^i);
L0 = {y, z} | apply(6, i -> x - z6^i*y) | apply(6, j -> z - z6^j*x);
analyse("omit {1,w6}   (pasted example)", L0 | apply({2,3,4,5}, k -> y - z6^k*z));
analyse("omit {1,w6^2}", L0 | apply({1,3,4,5}, k -> y - z6^k*z));
analyse("omit {1,w6^3}", L0 | apply({1,2,4,5}, k -> y - z6^k*z));
