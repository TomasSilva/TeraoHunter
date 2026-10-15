#!/usr/bin/env python3
"""The 236 lattice classes of the F_7 four-pencil run (work/build/fp7/fp_4444_7_6.out; none realisable
over Q-bar). Every F_7-rational member has a Jacobian syzygy of degree 7 < 8 = d1 (Frobenius bound,
TC_CHARP.md), so it is non-free; a characteristic-7 counterexample would need a FREE realisation over
F_7-bar. Compute each realisation space over F_7-bar (minimal primes over ZZ/7, same frame as
exact_tc_lattice.py) and test freeness directly on every component (resolution of the Jacobian ideal
over the residue field; positive-dimensional components over the function field)."""
import sys, subprocess, time
from collections import OrderedDict
sys.argv = ['x', '17', '8', 'dummy.jsonl', '0', '1']
sys.path.insert(0, '../tc15'); sys.path.insert(0, '../build')
import exact_tc_lattice as E
from census5_post import lattice, canon

M2 = '/opt/homebrew/bin/M2'
n, p = 17, 7

def lines4(line):
    t = line.split(); sizes = [int(x) for x in t[2:6]]; q = int(t[7]); mu = int(t[9])
    L = [(0, 0, 1)]
    for name, s, (dp, dq) in zip("ABCD", sizes, [(1, 0), (0, 1), (1, 1), (mu, 1)]):
        i = t.index(name)
        for v in t[i+1:i+1+s]: L.append((dq % q, (-dp) % q, (-int(v)) % q))
    return q, L

classes = OrderedDict()
for line in open('../build/fp7/fp_4444_7_6.out'):
    if not line.startswith('FOUND'): continue
    q, L = lines4(line); pts = lattice(L, q)
    assert sum(len(P) - 1 for P in pts) == 80
    key = canon(len(L), pts)
    if key not in classes: classes[key] = [sorted(P) for P in pts if len(P) >= 3]
print('classes', len(classes), flush=True)

def script(items):
    out = ['cr = (p, q) -> matrix{{p_(0,1)*q_(0,2)-p_(0,2)*q_(0,1), p_(0,2)*q_(0,0)-p_(0,0)*q_(0,2), p_(0,0)*q_(0,1)-p_(0,1)*q_(0,0)}};']
    for idx, blocks in items:
        frame, nvar, lines, eqs, bad, pts = E.setup(blocks)
        s = ['R = ZZ/%d[%s];' % (p, ','.join('u_%d' % i for i in range(nvar)) if nvar else 'u_0')]
        s.append('P%d = matrix{{1_R,0,0}}; P%d = matrix{{0_R,1,0}}; P%d = matrix{{0_R,0,1}}; P%d = matrix{{1_R,1,1}};' % frame)
        s += lines
        if nvar == 0: eqs = eqs + ['u_0']
        s.append('I = ideal(%s);' % ','.join(eqs))
        s.append('if I == ideal(1_R) then print("LAT %d NONREAL") else (' % idx)
        s.append('bad = select({%s}, f -> f != 0); zerobad = %d - #bad;' % (','.join(bad), len(bad)))
        s.append('cs = minimalPrimes I;')
        s.append('good = if zerobad > 0 then {} else select(cs, c -> all(bad, f -> f % c != 0));')
        s.append('print("LAT %d components " | toString(#cs) | " good " | toString(#good) | " dims " | toString(apply(good, dim)) | " degrees " | toString(apply(good, degree)));' % idx)
        s.append('Ls = {%s};' % ','.join('P%d' % e for e in range(n)))
        s.append('''for c in good do (
   F := if dim c == 0 then toField(R / c) else frac(R / c); T := F[x,y,z];
   f := product(Ls, m -> (v := flatten entries sub(m, F); v#0*x + v#1*y + v#2*z));
   J := ideal(diff(x, f), diff(y, f), diff(z, f)); C := res J;
   print("  COMP dim " | toString(dim c) | " degree " | toString(degree c) | " pd " | toString(length C) | " syzygies " | toString(sort apply(flatten degrees C_2, d -> d - 16)) | (if length C == 2 then " FREE  !!!" else " NONFREE"));
));''')
        out.append('\n'.join(s))
    return '\n'.join(out) + '\nprint "DONE";\n'

if __name__ == '__main__':
    items = list(enumerate(classes.values())); B = 12; jobs = []
    for start in range(0, len(items), B):
        path = 'charp7_lift_%03d.m2' % start
        open(path, 'w').write(script(items[start:start + B])); jobs.append(path)
    t0 = time.time(); running = []; outs = {}
    while jobs or running:
        while jobs and len(running) < 4:
            path = jobs.pop(0); running.append((path, subprocess.Popen([M2, '--script', path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)))
        path, proc = running.pop(0); so, se = proc.communicate(); outs[path] = so
        if 'DONE' not in so: outs[path] += 'ERROR ' + path + ' ' + se[-1500:] + '\n'
        print('done', path, '%.0fs' % (time.time() - t0), flush=True)
    with open('charp7_lift.out', 'w') as f:
        for path in sorted(outs): f.write(outs[path])
