load "jump5.m2";
(
k = value(scriptCommandLine#1);
F = jumpCurve k;
fac = toList factor F;
print("k " | toString k | ": " | toString(#fac) | " factors");
scan(fac, p -> (h := p#0; print("   deg " | toString(first degree h) | " (a-deg " | toString(degree(a, h)) | ", b-deg " | toString(degree(b, h)) | ") ^" | toString(p#1) | " : " | (if first degree h <= 4 then toString h else "..."))));
);
