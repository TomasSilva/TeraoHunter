-- Realisation space of a 4-pencil lattice.  Input via the variables sizes, lam0, A0, B0, C0, D0 (rationals).
-- Unknown offsets a_i, b_i, c_i, d_i and the cross-ratio t; pencil directions (1:0:0),(0:1:0),(1:1:0),(t:1:0).
-- Normalisation: a_0 = 0, b_0 = 0 (translations), a_1 = 1 (scaling)  [if sizes#0 >= 2].
realSpace = (sizes, lam0, A0, B0, C0, D0) -> (
  nA := sizes#0; nB := sizes#1; nC := sizes#2; nD := sizes#3;
  vars0 := (apply(nA, i -> "a"|toString i)) | (apply(nB, i -> "b"|toString i)) | (apply(nC, i -> "c"|toString i)) | (apply(nD, i -> "d"|toString i)) | {"t"};
  Rr := QQ(monoid[apply(vars0, v -> getSymbol v)]);
  gensR := gens Rr;
  a := i -> gensR#i; b := i -> gensR#(nA+i); c := i -> gensR#(nA+nB+i); d := i -> gensR#(nA+nB+nC+i); t := gensR#(nA+nB+nC+nD);
  -- lines as coefficient vectors (q, -p, -v) with directions (p,q): (1,0),(0,1),(1,1),(t,1)
  lines := {{0_Rr,0_Rr,1_Rr}} | apply(nA, i -> {0_Rr, -1_Rr, -a i}) | apply(nB, i -> {1_Rr, 0_Rr, -b i}) | apply(nC, i -> {1_Rr, -1_Rr, -c i}) | apply(nD, i -> {1_Rr, -t, -d i});
  n := #lines;
  -- the given realisation
  vals := {0_QQ, 0_QQ, 1_QQ}; -- dummy
  subs0 := apply(nA, i -> a i => A0#i) | apply(nB, i -> b i => B0#i) | apply(nC, i -> c i => C0#i) | apply(nD, i -> d i => D0#i) | {t => lam0};
  ev := f -> sub(f, subs0);
  dets := new MutableHashTable;
  for T in subsets(n,3) do dets#T = det matrix apply(T, i -> lines#i);
  conc := select(subsets(n,3), T -> ev(dets#T) == 0);
  nonc := select(subsets(n,3), T -> ev(dets#T) != 0);
  I := ideal apply(conc, T -> dets#T);
  I = I + ideal(a 0, b 0) + (if nA >= 2 then ideal(a 1 - 1) else ideal());
  -- degeneracies: equal offsets in a pencil, t in {0,1}
  degs := apply(subsets(nA,2), p -> a(p#0) - a(p#1)) | apply(subsets(nB,2), p -> b(p#0) - b(p#1)) | apply(subsets(nC,2), p -> c(p#0) - c(p#1)) | apply(subsets(nD,2), p -> d(p#0) - d(p#1)) | {t, t-1};
  for g in degs do I = saturate(I, g);
  comps := minimalPrimes I;
  print("   incidence ideal: " | toString(#conc) | " concurrencies; components: " | toString(#comps));
  for Cc in comps do (
    extra := select(nonc, T -> (dets#T) % Cc == 0);
    tconst := (dim(Cc + ideal(t - lam0)) < dim Cc) or (dim Cc == 0);
    print("   component dim " | toString dim Cc | ", degree " | toString degree Cc | ", extra incidences " | toString(#extra)
          | (if #extra > 0 then " (coarser lattice: discard)" else "") | ", t constant: " | toString(tconst)
          | (if dim Cc > 0 and #extra == 0 then (", t on component: " | toString(eliminate(Cc, delete(t, gensR)))) else "")));
  comps);
