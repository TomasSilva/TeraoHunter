-- control: the free (8,8) harmonic grid G17 of TC_BUILD.md section 3, reduced mod p
for p in {7, 11, 13} do (
   S := ZZ/p[x,y,z];
   L := {z} | apply(4, a -> y - a*z) | apply({2,3,4,5}, b -> x - b*z) | apply(4, c -> x - y - c*z) | apply({3,4,5,6}, d -> x + y - d*z);
   f := product L; J := ideal(diff(x, f), diff(y, f), diff(z, f)); C := res J;
   print("G17 mod " | toString p | ": " | toString(#L) | " lines, pd " | toString(length C) | " syzygies " | toString(sort apply(flatten degrees C_2, d -> d - 16))));
