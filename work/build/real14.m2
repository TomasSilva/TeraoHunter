-- Realisation space of L(F15 minus a line): 6 points, the surviving concurrent matchings.
load "famtools.m2";
Q = {{1_R,0_R,0_R}, {0_R,1_R,0_R}, {0_R,0_R,1_R}, {1_R,1_R,1_R},
     {l1^2, l1*l0, 2*l1*l0-l0^2}, {l1*(2*l1-l0), l1^2, (2*l1-l0)*l0}};
edges = subsets(6,2);
F = apply(edges, p -> crossv(Q#(p#0), Q#(p#1)));
(pts, degF) = latticeOf F;
trip = select(pts, S -> #S == 3);
print("triple points (as edge triples): " | toString apply(trip, S -> apply(S, i -> edges#i)));
sH = apply(15, a -> #select(pts, S -> member(a,S)));
print("points per line: " | toString sH);
-- generic 6 points
S6 = QQ[a,b,c,d];
Qg = {{1,0,0},{0,1,0},{0,0,1},{1,1,1},{a,b,1},{c,d,1}};
Qg = apply(Qg, v -> apply(v, x -> x*1_S6));
cr = (u,v) -> {u#1*v#2-u#2*v#1, u#2*v#0-u#0*v#2, u#0*v#1-u#1*v#0};
Lg = apply(edges, p -> cr(Qg#(p#0), Qg#(p#1)));
conc = S -> det matrix apply(S, i -> Lg#i);
-- degeneracies: collinear triples of Q's, and concurrencies of edge-triples that are not triple points of F
coll = product apply(subsets(6,3), s -> det matrix apply(s, i -> Qg#i));
allTrip = select(subsets(15,3), S -> #unique flatten apply(S, i -> edges#i) == 6);  -- matchings only (other triples share a Q)
badTrip = select(allTrip, S -> not member(S, trip));
degen = coll * product apply(badTrip, conc);
for dropLine in {0, 5} do (
  keep := select(trip, S -> not member(dropLine, S));
  I := ideal apply(keep, conc);
  I = saturate(I, degen);
  print("== drop line " | toString(edges#dropLine) | " (" | toString(sH#dropLine) | " points): " | toString(#keep) | " concurrency conditions");
  print("   dim = " | toString dim I | ", degree = " | toString degree I);
  comps := minimalPrimes I;
  for C in comps do (
     lost := select(trip, S -> member(dropLine, S));
     print("   component dim " | toString dim C | ", degree " | toString degree C | ": lost concurrencies automatic? "
           | toString apply(lost, S -> (conc S) % C == 0))));
