load "jump5.m2";
-- R17a: points (0:1:0) x=0 mult 3, (1:0:0) y=0 mult 4, (1:1:0) mult 3, (4/3:1:0) mult 1, (2:1:0) mult 5
FA = jumpCurve {3, 4, 3, 1, 5};
print("R17a: F(4/3, 2) = " | toString sub(FA, {a => 4/3, b => 2}) | "   control F(5/7, 3) = " | toString sub(FA, {a => 5/7, b => 3}));
-- R17b: (0:1:0) mult 4, (1:0:0) mult 5, (1:1:0) mult 1, (-1/2:1:0) mult 1, (-1:1:0) mult 5
FB = jumpCurve {4, 5, 1, 1, 5};
print("R17b: F(-1/2, -1) = " | toString sub(FB, {a => -1/2, b => -1}) | "   control F(5/7, 3) = " | toString sub(FB, {a => 5/7, b => 3}));
