#!/usr/bin/env python3
"""Realisation spaces over F_p-bar (p = 11, 13) of the lattices found by the pencil censuses (the six
five-pencil classes and E17), with a direct freeness test on every component: the arrangement over
the residue field of the component, minimal free resolution of its Jacobian ideal (free iff pd = 2).
Same frame and parametrisation as work/tc15/exact_tc_lattice.py, with ZZ/p in place of QQ."""
import sys, json, subprocess
sys.argv = ['x', '17', '8', 'dummy.jsonl', '0', '1']
sys.path.insert(0, '../tc15'); sys.path.insert(0, '../build')
import exact_tc_lattice as E
from census5_post import lattice

M2 = '/opt/homebrew/bin/M2'
n = 17

def e17_blocks():
    line = next(l for l in open('../build/fp_runs/fp4444_p13_m4.out') if l.startswith('FOUND'))
    t = line.split(); sizes = [int(x) for x in t[2:6]]; p = int(t[7]); mu = int(t[9])
    L = [(0, 0, 1)]
    for name, s, (dp, dq) in zip("ABCD", sizes, [(1, 0), (0, 1), (1, 1), (mu, 1)]):
        i = t.index(name)
        for v in t[i+1:i+1+s]: L.append((dq % p, (-dp) % p, (-int(v)) % p))
    return [sorted(P) for P in lattice(L, p) if len(P) >= 3]

lats = [('c5-%d' % i, json.loads(l)['blocks']) for i, l in enumerate(open('../build/census5_classes.jsonl'))]
lats.append(('E17', e17_blocks()))

out = ['cr = (p, q) -> matrix{{p_(0,1)*q_(0,2)-p_(0,2)*q_(0,1), p_(0,2)*q_(0,0)-p_(0,0)*q_(0,2), p_(0,0)*q_(0,1)-p_(0,1)*q_(0,0)}};']
for p in (11, 13):
    for lab, blocks in lats:
        frame, nvar, lines, eqs, bad, pts = E.setup(blocks)
        s = ['R = ZZ/%d[%s];' % (p, ','.join('u_%d' % i for i in range(nvar)) if nvar else 'u_0')]
        s.append('P%d = matrix{{1_R,0,0}}; P%d = matrix{{0_R,1,0}}; P%d = matrix{{0_R,0,1}}; P%d = matrix{{1_R,1,1}};' % frame)
        s += lines
        if nvar == 0: eqs = eqs + ['u_0']
        s.append('I = ideal(%s);' % ','.join(eqs))
        s.append('if I == ideal(1_R) then print("LAT %s p=%d NONREAL") else (' % (lab, p))
        s.append('bad = select({%s}, f -> f != 0); zerobad = %d - #bad;' % (','.join(bad), len(bad)))
        s.append('cs = minimalPrimes I;')
        s.append('good = if zerobad > 0 then {} else select(cs, c -> all(bad, f -> f % c != 0));')
        s.append('print("LAT %s p=%d components " | toString(#cs) | " good " | toString(#good) | " dims " | toString(apply(good, dim)) | " zerobad " | toString zerobad);' % (lab, p))
        s.append('Ls = {%s};' % ','.join('P%d' % e for e in range(n)))
        s.append('''for c in good do (
   if dim c == 0 then (
      F := toField(R / c); T := F[x,y,z];
      f := product(Ls, m -> (v := flatten entries sub(m, F); v#0*x + v#1*y + v#2*z));
      J := ideal(diff(x, f), diff(y, f), diff(z, f)); C := res J;
      print("  COMP dim 0 degree " | toString(degree c) | " (residue field F_" | toString(%d) | "^" | toString(degree c) | ") pd " | toString(length C) | " syzygies " | toString(sort apply(flatten degrees C_2, d -> d - 16)) | (if length C == 2 then " FREE" else " NONFREE"));
   ) else print("  COMP dim " | toString(dim c) | " (positive-dimensional: not handled)");
));''' % p)
        out.append('\n'.join(s))
open('charp_lift.m2', 'w').write('\n'.join(out) + '\nprint "DONE";\n')
res = subprocess.run([M2, '--script', 'charp_lift.m2'], capture_output=True, text=True)
print(res.stdout); print(res.stderr[-3000:], file=sys.stderr)
