import time, random
import numpy as np
from flint import nmod_mat
p = 2147483629
for (r,c) in [(1000,800),(2700,1800),(4000,3000)]:
    A = np.random.randint(0, p, size=(r, c), dtype=np.int64)
    # make it rank-deficient
    B = A[:, :c//2]
    A = np.concatenate([B, (B*3) % p], axis=1)
    t = time.time()
    M = nmod_mat(A.tolist(), p)
    t1 = time.time()
    rk = M.rank()
    t2 = time.time()
    X, nul = M.nullspace()
    t3 = time.time()
    print(r, c, 'conv %.2f rank %.2f null %.2f' % (t1-t, t2-t1, t3-t2), rk, nul)
