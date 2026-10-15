load "specialise.m2";
fam15 = (dd) -> (K := ring dd; o := 1_K;
  Q := {{o,0*o,0*o},{0*o,o,0*o},{0*o,0*o,o},{o,o,o},{dd^2,dd,2*dd-1},{dd*(2*dd-1),dd^2,2*dd-1}};
  apply(subsets(6,2), p -> crossK(Q#(p#0), Q#(p#1))));
analyseAt("15-line, generic d = 7/5", fam15(7/5));
Ka = toField(QQ[a]/(a^2-3*a+1)); use Ka;
analyseAt("15-line, d root of d^2-3d+1 (remote for lines 6,7,12)", fam15(a));
Kb = toField(QQ[bb]/(bb^2+bb-1)); use Kb;
analyseAt("15-line, d root of d^2+d-1 (remote for lines 6,7,12)", fam15(bb));
