R = QQ[a,b,d,e];
Q1 = matrix{{1,0,0}}; Q2 = matrix{{0,1,0}}; Q3 = matrix{{0,0,1}}; Q4 = matrix{{1,1,1}};
Q5 = matrix{{a,b,1}}; Q6 = matrix{{d,e,1}};
Q = {Q1,Q2,Q3,Q4,Q5,Q6};
cross = (u,v) -> matrix{{u_(0,1)*v_(0,2)-u_(0,2)*v_(0,1), u_(0,2)*v_(0,0)-u_(0,0)*v_(0,2), u_(0,0)*v_(0,1)-u_(0,1)*v_(0,0)}};
line = (i,j) -> cross(Q#(i-1), Q#(j-1));
conc = (m) -> det( (line(m#0#0,m#0#1)) || (line(m#1#0,m#1#1)) || (line(m#2#0,m#2#1)) );
M = {{{1,2},{3,4},{5,6}}, {{1,3},{2,6},{4,5}}, {{1,6},{2,5},{3,4}}, {{1,4},{2,6},{3,5}}, {{1,5},{2,3},{4,6}}, {{1,5},{2,4},{3,6}}};
I = ideal apply(M, conc);
print("gens:"); print toString I_*;
-- degeneracy: collinear triples of the Q's, coincident points
coll = {};
for i from 0 to 5 do for j from i+1 to 5 do for k from j+1 to 5 do coll = append(coll, det(Q#i || Q#j || Q#k));
-- the other 9 matchings must not be concurrent
allM = {};
for s in subsets({1,2,3,4,5,6}, 2) do (
  rest = toList(set{1,2,3,4,5,6} - set s);
  for t in subsets(rest, 2) do (
    u = toList(set rest - set t);
    m = sort {sort s, sort t, sort u};
    allM = append(allM, m)));
allM = unique allM;
others = select(allM, m -> not member(m, apply(M, x -> sort apply(x, sort))));
print("number of matchings: " | toString(#allM) | ", others: " | toString(#others));
bad = product(coll) * product(apply(others, conc));
J = saturate(I, bad);
print("dim of realization space (affine, frame fixed): " | toString dim J);
print("degree: " | toString degree J);
cs = decompose J;
for c in cs do (print("component: dim " | toString dim c | "  gens " | toString c_*));
