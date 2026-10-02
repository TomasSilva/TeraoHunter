# Odd q: A* = L u (PG(3,q) \ (Q u T_p)), Q = y1^2 - n y2^2 + y3 y4 elliptic (n non-square),
# p = (0,0,1,0), T_p = {y4 = 0}, L = {y2 = y4 = 0}.  Check |A*| = q^3 - q^2 + q + 1, blocking, minimal.
import sys
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter/work')
from covering_check import proj_points, lines_of

for q in [3, 5, 7]:
    sq = {(x * x) % q for x in range(1, q)}
    n = next(x for x in range(2, q) if x not in sq)
    Qf = lambda y: (y[0] ** 2 - n * y[1] ** 2 + y[2] * y[3]) % q
    pts = proj_points(q, 4)
    onQ = [y for y in pts if Qf(y) == 0]
    S = {y for y in pts if (y[1] == 0 and y[3] == 0)} | {y for y in pts if y[3] != 0 and Qf(y) != 0}
    lines = lines_of(q, 4)
    blocking = all(L & S for L in lines)
    tang = {h: 0 for h in S}
    for L in lines:
        I = L & S
        if len(I) == 1:
            tang[next(iter(I))] += 1
    print('q=%d: |Q|=%d (q^2+1=%d)  |A*|=%d (q^3-q^2+q+1=%d)  blocking=%s  minimal=%s  d=%d >= q^2: %s' % (
        q, len(onQ), q * q + 1, len(S), q ** 3 - q * q + q + 1, blocking, min(tang.values()) >= 1,
        len(S) - (q * q + q + 1), len(S) - (q * q + q + 1) >= q * q))
