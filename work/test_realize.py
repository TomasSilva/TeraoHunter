import sys, time
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.realize import Realizer
from hunter.ffree import free_rank3_ziegler, free_simple
for q in [5, 7, 13]:
    A3 = np.array([[1,0,0],[0,1,0],[0,0,1],[1,q-1,0],[1,0,q-1],[0,1,q-1]])
    t = time.time(); sols = Realizer(A3, q).solutions()
    print('A3 over F_%d: %d realizations (%.2fs)' % (q, len(sols), time.time() - t))
# dual Hesse over F_7 and F_13: x^3-y^3 etc with cube roots of unity
for q in [7, 13]:
    w = next(a for a in range(2, q) if pow(a, 3, q) == 1)
    roots = [1, w, w * w % q]
    H = [[1, (-r) % q, 0] for r in roots] + [[0, 1, (-r) % q] for r in roots] + [[1, 0, (-r) % q] for r in roots]
    t = time.time(); sols = Realizer(np.array(H), q).solutions()
    print('dual Hesse over F_%d: %d realizations, freeness:' % (q, len(sols)), [free_rank3_ziegler(s, q) for s in sols], '(%.2fs)' % (time.time() - t))
# braid A4 in rank 4
for q in [5, 7]:
    braid = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]] + [[1 if k==i else (q-1 if k==j else 0) for k in range(4)] for i in range(4) for j in range(i+1,4)]
    t = time.time(); sols = Realizer(np.array(braid), q).solutions()
    print('A4 braid over F_%d: %d realizations (%.2fs)' % (q, len(sols), time.time() - t), [free_simple(s, q)[0] for s in sols])
