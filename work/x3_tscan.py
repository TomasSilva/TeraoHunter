import sys
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import multi_freeness
p = 1009
def X3(t):
    return np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,t % p]], dtype=np.int64) % p
for m in [(2,2,2,1,1,1),(3,3,3,1,1,1),(4,4,4,1,1,1)]:
    bad = []
    for t in range(1, p - 1):
        st, g = multi_freeness(X3(t), m, p)
        if st != 'free':
            bad.append((t, st, g))
    print(m, 'nonfree t:', bad)
