load "specialise.m2";
fam13 = (lam) -> (K := ring lam; o := 1_K;
 {{-o,-o,lam+1},{0*o,2*o,lam+1},{o,-o,lam+1},{o,-o,lam-2},{-o,-o,lam-2},{0*o,2*o,lam-2},
  {0*o,2*o,-2*lam+1},{o,-o,-2*lam+1},{-o,-o,-2*lam+1},{1-lam,lam+1,-lam^2+lam-1},
  {lam,lam-2,-lam^2+lam-1},{-o,1-2*lam,-lam^2+lam-1},{0*o,0*o,o}});
analyseAt("13-line, generic lambda = 7/5", fam13(7/5));
analyseAt("13-line, lambda = -1 (remote degeneration for lines 4..13)", fam13(-1/1));
Ki = toField(QQ[jj]/(jj^2+1)); use Ki;
analyseAt("13-line, lambda = i (in W)", fam13(jj));
Kf = toField(QQ[ph]/(ph^2-ph-1)); use Kf;
analyseAt("13-line, lambda = golden ratio (in W)", fam13(ph));
analyseAt("13-line, lambda = 2", fam13(2/1));
analyseAt("13-line, lambda = 1/2", fam13(1/2));
