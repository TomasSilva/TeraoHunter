# Rank 3 over F_{Q}, Q = r^2: lines a x + b y + c z = 0 with a^{r+1} + b^{r+1} + c^{r+1} = 0
# (dual of the Hermitian unital).  Predicted: free, exps (1, Q, r^3 - r^2), not DF.
import sys
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
from hunter.fq import GF, LineArr

for r in [2, 3, 4]:
    Q = r * r
    F = GF(Q)

    def pw(a, e):
        out = 1
        for _ in range(e):
            out = int(F.mul[out, a])
        return out
    pts = F.points(3)
    U = [v for v in pts if int(F.add[F.add[pw(v[0], r + 1), pw(v[1], r + 1)], pw(v[2], r + 1)]) == 0]
    A = LineArr(F, U)
    print('Q=%d: %d lines (r^3+1=%d), chi %s, exps %s, free %s, DF %s' % (
        Q, len(U), r ** 3 + 1, A.char_poly(), A.exponents(), A.free(), A.is_df()), flush=True)
