load "specialise.m2";
-- the 15 forms of [ACKN, Sec. 7] at t = 7/5 (lines are kernels of these forms)
R15 = t -> {{1,0,0},{1,1,0},{1,0,1},{1,1,1},{1,t,1},{0,1,0},{2,1,1},{t+1,t,1},{t+1,1,1},{2*t,t,1},{1,-t+1,1},
           {-3*t+1,t^2-3*t+1,-t},{3*t-1,t,t},{-3*t+1,-t^2,-t},{3*t-1,2*t-1,t}};
analyseAt("ACKN 15 at t = 7/5", apply(R15(7/5), v -> apply(v, c -> c * 1/1)));
