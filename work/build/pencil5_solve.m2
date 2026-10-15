-- usage: M2 --script pencil5_solve.m2 <classfile>
fn = scriptCommandLine#1;
lns = lines get fn;
sz = apply(drop(separate(" ", lns#0), 1), value);
bl = apply(select(lns, l -> match("^point ", l)), l -> apply(drop(separate(" ", l), 1), value));
names = {"a","b","c","e","d"};
vars0 = flatten apply(5, k -> apply(sz#k, i -> names#k | toString i)) | {"t1", "t2"};
Rr = QQ(monoid[apply(vars0, v -> getSymbol v), MonomialOrder => Lex]);
g = gens Rr; t1 = g#(#g - 2); t2 = g#(#g - 1);
off = (k, i) -> g#(sum(take(sz, k)) + i);
dirs = {{1_Rr, 0_Rr}, {0_Rr, 1_Rr}, {1_Rr, 1_Rr}, {t1, 1_Rr}, {t2, 1_Rr}};
lineV = {{0_Rr, 0_Rr, 1_Rr}} | flatten apply(5, k -> apply(sz#k, i -> {(dirs#k)#1, -(dirs#k)#0, -off(k,i)}));
concT = unique flatten apply(bl, P -> subsets(P, 3));
I = ideal apply(concT, T -> det matrix apply(T, i -> lineV#i)) + ideal(off(0,0), off(0,1) - 1, off(1,0));
degs = flatten apply(5, k -> apply(subsets(sz#k, 2), q -> off(k, q#0) - off(k, q#1))) | {t1, t1 - 1, t2, t2 - 1, t1 - t2};
for h in degs do I = saturate(I, h);
print("dim " | toString dim I | ", degree " | toString degree I);
G = flatten entries gens gb I;
for f in G do print("   " | toString f);
