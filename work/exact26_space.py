# Realisation space of the 26-plane lattice inside the 3-parameter family (a, b, c), over QQ.
# 1. read the lattice from the F_7 example (a, b, c) = (3, 3, 3): planes as point sets, rank-2 and rank-3 flats;
# 2. in Macaulay2 impose every dependency of that lattice on symbolic (a, b, c), saturate by the
#    non-coplanarity of the other point quadruples, decompose;
# 3. for each solution compare the full lattice (rank-2/rank-3 flats of the plane arrangement) with F_7's.
import sys, itertools, subprocess, json
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter'); sys.path.insert(0, '/Users/agent/Projects/TeraoHunter/work')
import numpy as np
from family26 import config, spanned
from hunter.core import Lattice
p = 7
N, sets = spanned(config(3, 3, 3, p), p)
L = Lattice(N, p)
print('F_7: planes', len(N), 'exps', L.exponents(), flush=True)
flats2 = [sorted(S) for S in L.flats[2] if len(S) >= 3]
flats3 = [sorted(S) for S in L.flats[3] if len(S) >= 4]
plane_sets = [sorted(s) for s in sets]
cop = [S for S in plane_sets if len(S) >= 4]
in_plane = lambda Q: any(set(Q) <= set(S) for S in plane_sets)
noncop = [Q for Q in itertools.combinations(range(8), 4) if not in_plane(Q)]
print('plane point-sets with >=4 points:', len(cop), ' non-coplanar quadruples:', len(noncop),
      ' rank-2 flats (>=3 planes):', len(flats2), ' rank-3 flats (>=4 planes):', len(flats3), flush=True)
m2 = ['R = QQ[a,b,c];',
      'pt = {matrix{{1,0,0,0}}, matrix{{0,1,0,0}}, matrix{{1,a,0,0}}, matrix{{1,1,1,1}}, matrix{{1,1,c,c}}, matrix{{0,0,1,0}}, matrix{{0,0,0,1}}, matrix{{0,0,1,b}}};',
      'pt = apply(pt, v -> sub(v, R));',
      # normal of the plane through 3 points: signed 3x3 minors
      'nrm = (i,j,k) -> (M := pt#i || pt#j || pt#k; matrix{apply(4, c -> (-1)^c * det submatrix(M, , toList(set(0..3) - set{c})))});']
eqs = []
for S in cop:
    for extra in S[3:]:
        eqs.append('det(pt#%d || pt#%d || pt#%d || pt#%d)' % (S[0], S[1], S[2], extra))
m2.append('planeNormals = {%s};' % ','.join('nrm(%d,%d,%d)' % tuple(S[:3]) for S in plane_sets))
for F in flats2:
    m2.append('I2_%d = minors(3, fold((x,y) -> x || y, apply({%s}, i -> planeNormals#i)));' % (len(eqs), ','.join(map(str, F))))
    eqs.append('I2_%d' % len(eqs))
for F in flats3:
    for extra in F[3:]:
        eqs.append('det(planeNormals#%d || planeNormals#%d || planeNormals#%d || planeNormals#%d)' % (F[0], F[1], F[2], extra))
gens = ['ideal(%s)' % e if not e.startswith('I2_') else e for e in eqs]
m2.append('I = %s;' % ' + '.join(gens))
m2.append('bad = {%s};' % ','.join('det(pt#%d || pt#%d || pt#%d || pt#%d)' % Q for Q in noncop))
m2.append('''G = I;
for f in sort(bad, f -> first degree f) do if G != ideal(1_R) then G = saturate(G, f);
print("saturated realisation ideal: " | toString G);
cs = decompose G;
print("components: " | toString(#cs));
for C in cs do print("  " | toString C_* | "   dim " | toString dim C | " degree " | toString degree C);''')
open('/Users/agent/Projects/TeraoHunter/work/exact26_space.m2', 'w').write('\n'.join(m2) + '\n')
out = subprocess.run(['/opt/homebrew/bin/M2', '--script', '/Users/agent/Projects/TeraoHunter/work/exact26_space.m2'], capture_output=True, text=True, timeout=3000)
print(out.stdout, out.stderr[-1500:])
