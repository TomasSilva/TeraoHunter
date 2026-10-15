-- sanity check of ag8_4pencil.c on given hits: lattice numbers and freeness over F_8
R0 = ZZ/2[t]; K = toField(R0/ideal(t^3+t+1)); S = K[x,y,z];
el = m -> sum(3, i -> if (m >> i) % 2 == 1 then t^i else 0_K);
offs = mask -> select(toList(0..7), v -> (mask >> v) % 2 == 1);
arr = (m1, m2, m3, m4) -> {z} | apply(offs m1, v -> y + el(v)*z) | apply(offs m2, v -> x + el(v)*z) |
      apply(offs m3, v -> x + y + el(v)*z) | apply(offs m4, v -> x + t*y + el(v)*z);
test = (m1, m2, m3, m4) -> (L := arr(m1, m2, m3, m4); f := product L;
   J := ideal(diff(x, f), diff(y, f), diff(z, f)); C := res J;
   print(toString(#L) | " lines, pd " | toString(length C) | ", syzygies " | toString(sort apply(flatten degrees C_2, d -> d - (#L - 1)))));
test(0xfe, 0xfe, 0xe0, 0x80);
test(0xfe, 0xfe, 0xe0, 0x40);
