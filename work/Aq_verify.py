import sys
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter'); sys.path.insert(0, '/Users/agent/Projects/TeraoHunter/work')
import sympy as sp
from hunter.fq import GF, LineArr
from test_fq import m2_free_fq
# 1. the Saito determinant identity
a, b, c = sp.symbols('a b c')
x, y, z, qq = sp.symbols('x y z q', positive=True)
lhs = a*b**2 - a*c**2 - a**2*b + a**2*c + b*c**2 - b**2*c
print('identity det = (a-b)(b-c)(c-a):', sp.expand(lhs - (a - b)*(b - c)*(c - a)) == 0)
for q in [4, 5, 7]:
    M = sp.Matrix([[x, x**q, (y*z)**(q-2)], [y, y**q, (x*z)**(q-2)], [z, z**q, (x*y)**(q-2)]])
    d = sp.expand(M.det())
    target = sp.expand((x**(q-1) - y**(q-1))*(y**(q-1) - z**(q-1))*(z**(q-1) - x**(q-1)))
    print('q=%d: det(E, Frobenius, theta2) == prod (x^(q-1)-y^(q-1))...:' % q, sp.expand(d - target) == 0)
# 2. the arrangement A_q over F_q: lattice, chi, freeness, DF
for q in [3, 4, 5, 7, 8, 9, 11, 13, 16]:
    F = GF(q)
    lines = []
    for i, j in [(0, 1), (1, 2), (0, 2)]:
        for a_ in range(1, q):
            v = [0, 0, 0]; v[i] = 1; v[j] = int(F.neg[a_]); lines.append(tuple(v))
    A = LineArr(F, lines)
    pts_per_line = sorted(set(sum(1 for S in A.flats if i in S) for i in range(A.n)))
    mults = sorted(set(len(S) for S in A.flats))
    print(f'q={q}: n={A.n} exps={A.exponents()} free={A.free()} DF={A.is_df()} points/line={pts_per_line} point multiplicities={mults}', flush=True)
    if q in (4, 5, 7, 8):
        print('      Macaulay2 free:', m2_free_fq(F, lines), flush=True)
