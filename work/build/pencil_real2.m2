-- usage: M2 --script pencil_real2.m2 <classfile>
-- Realisation space over Q of a 4-pencil lattice given by its multiple points (line indices:
-- 0 = H, then pencils A,B,C,D in order). Unknown offsets and cross-ratio t.
fn = scriptCommandLine#1;
lns = lines get fn;
sizes = apply(drop(separate(" ", lns#0), 1), value);
ptsL = apply(select(lns, l -> match("^point ", l)), l -> apply(drop(separate(" ", l), 1), value));
nA = sizes#0; nB = sizes#1; nC = sizes#2; nD = sizes#3; n = 1 + nA + nB + nC + nD;
vars0 := (apply(nA, i -> "a"|toString i)) | (apply(nB, i -> "b"|toString i)) | (apply(nC, i -> "c"|toString i)) | (apply(nD, i -> "d"|toString i)) | {"t"};
Rr = QQ(monoid[apply(vars0, v -> getSymbol v)]);
g = gens Rr;
a = i -> g#i; b = i -> g#(nA+i); c = i -> g#(nA+nB+i); d = i -> g#(nA+nB+nC+i); t = g#(nA+nB+nC+nD);
lineV = {{0_Rr,0_Rr,1_Rr}} | apply(nA, i -> {0_Rr, -1_Rr, -a i}) | apply(nB, i -> {1_Rr, 0_Rr, -b i}) | apply(nC, i -> {1_Rr, -1_Rr, -c i}) | apply(nD, i -> {1_Rr, -t, -d i});
dets = new MutableHashTable;
for T in subsets(n,3) do dets#T = det matrix apply(T, i -> lineV#i);
concT = unique flatten apply(ptsL, P -> subsets(P, 3));
I = ideal apply(concT, T -> dets#T) + ideal(a 0, b 0) + (if nA >= 2 then ideal(a 1 - 1) else ideal(0_Rr));
degs = apply(subsets(nA,2), p -> a(p#0) - a(p#1)) | apply(subsets(nB,2), p -> b(p#0) - b(p#1)) | apply(subsets(nC,2), p -> c(p#0) - c(p#1)) | apply(subsets(nD,2), p -> d(p#0) - d(p#1)) | {t, t-1};
for h in degs do I = saturate(I, h);
if I == ideal(1_Rr) then (print(fn | ": NOT realisable over Q-bar (after saturation)"); exit 0);
comps = minimalPrimes I;
nonc = select(subsets(n,3), T -> not member(T, concT));
print(fn | ": sizes " | toString sizes | ", " | toString(#concT) | " concurrent triples, " | toString(#comps) | " components");
for Cc in comps do (
  extra := select(nonc, T -> (dets#T) % Cc == 0);
  dm := dim Cc;
  tinfo := if dm == 0 then ("t-values: " | toString(eliminate(Cc, delete(t, toList g)))) else (
     Jt := eliminate(Cc, delete(t, toList g));
     if Jt == 0 then "t NON-CONSTANT on this component" else ("t constant: " | toString Jt));
  print("   component dim " | toString dm | ", degree " | toString degree Cc | ", extra incidences " | toString(#extra) | (if #extra > 0 then " -> coarser lattice, discard" else " -> genuine") | "; " | tinfo));
