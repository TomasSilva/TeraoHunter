"""Dump one dominating free family from thmF_test.py (vector (1,2,2,3), jump at mu0 = 2) as M2 input:
the generic member at mu = 7/3 and the limit member at mu0 = 2."""
import random
from fractions import Fraction as Fr
import thmF_test as T
k, mu0 = (1, 2, 2, 3), 2
K = sum(k); d = K // 2; target = d * d + 2 * d
rng = random.Random(1)
found = 0
for t in range(20000):
    f = T.build(k, None, rng.random())
    if f is None:
        continue
    g1, g2 = Fr(rng.randint(5, 50), rng.randint(51, 97)) + 3, Fr(-rng.randint(5, 50), rng.randint(51, 97)) - 5
    ok1, P1 = T.lattice(f(g1)); ok2, P2 = T.lattice(f(g2))
    if not (ok1 and ok2) or T.structure(P1) != T.structure(P2) or T.b2(P1) != T.b2(P2) or T.b2(P1) != target:
        continue
    def m2(lines):
        return '{' + ', '.join('(%s)*x + (%s)*y + (%s)*z' % (str(a), str(b), str(c)) for a, b, c in lines) + '}'
    print('gen = ' + m2(f(Fr(7, 3))) + ';')
    print('lim = ' + m2(f(Fr(mu0))) + ';')
    break
