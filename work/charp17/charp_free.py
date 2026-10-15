#!/usr/bin/env python3
"""Direct freeness test over F_p of every pencil-census hit (17 lines, class (8,8)[B] window).
For each hit: f = product of its lines over ZZ/p, J = Jacobian ideal, minimal free resolution of S/J.
Free iff pd(S/J) = 2 (Saito; p does not divide 17, so f lies in J). Prints syzygy degrees."""
import sys, glob, subprocess, json
from collections import Counter, defaultdict
sys.path.insert(0, '../build')
from census5_post import parse as parse5, arrangement as arr5, lattice, canon

M2 = '/opt/homebrew/bin/M2'

def parse4(line):
    t = line.split(); sizes = [int(x) for x in t[2:6]]; p = int(t[7]); mu = int(t[9])
    offs = []
    for name, s in zip("ABCD", sizes):
        i = t.index(name); offs.append([int(x) for x in t[i+1:i+1+s]])
    L = [(0, 0, 1)]
    for (dp, dq), o in zip([(1, 0), (0, 1), (1, 1), (mu, 1)], offs):
        for v in o: L.append((dq % p, (-dp) % p, (-v) % p))
    return p, L

hits = []          # (source, class label, p, lines)
keys = {}
for f in sorted(glob.glob('../build/census5_runs/*.out')):
    for line in open(f):
        if not line.startswith('FOUND5'): continue
        S, p, l1, l2, offs = parse5(line); L = arr5(S, p, l1, l2, offs)
        key = canon(len(L), lattice(L, p)); keys.setdefault(key, len(keys))
        hits.append(('5pt', 'c5-%d' % keys[key], p, L))
for f in sorted(glob.glob('../build/fp_runs/*.out')):
    if '_p11_' not in f and '_p13_' not in f: continue
    for line in open(f):
        if not line.startswith('FOUND'): continue
        p, L = parse4(line); hits.append(('4pt', 'E17', p, L))
print('hits:', Counter((h[0], h[2]) for h in hits), flush=True)

results = []
for p in sorted(set(h[2] for h in hits)):
    mine = [(i, h) for i, h in enumerate(hits) if h[2] == p]
    s = ['S = ZZ/%d[x,y,z];' % p, 'H = {']
    s.append(',\n'.join('{%d,{%s}}' % (i, ','.join('{%d,%d,%d}' % l for l in h[3])) for i, h in mine))
    s.append('};')
    s.append('''scan(H, h -> (L := h#1; f := product(L, l -> l#0*x + l#1*y + l#2*z);
   J := ideal(diff(x, f), diff(y, f), diff(z, f)); C := res J;
   print("R " | toString(h#0) | " " | toString(length C) | " " | toString(sort apply(flatten degrees C_2, d -> d - 16)))));
print "DONE";''')
    path = 'charp_free_p%d.m2' % p
    open(path, 'w').write('\n'.join(s))
    out = subprocess.run([M2, '--script', path], capture_output=True, text=True).stdout
    assert 'DONE' in out, out[-2000:]
    for line in out.splitlines():
        if line.startswith('R '):
            _, i, pd, degs = line.split(' ', 3); results.append((int(i), int(pd), degs))

table = Counter()
for i, pd, degs in results:
    src, lab, p, L = hits[i]
    table[(lab, p, 'FREE' if pd == 2 else 'NONFREE', degs)] += 1
for k in sorted(table): print(k, table[k])
print('tested', len(results), 'of', len(hits), '; free:', sum(v for k, v in table.items() if k[2] == 'FREE'))
