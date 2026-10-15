# Reduce a G26 window deletion with a 5-point line mod p and compute its census coordinates (sizes order, l1, l2)
import sys
from itertools import permutations
p = 13; w = 4                     # w^2 - w + 1 = 0 mod 13
inv = lambda a: pow(a % p, p - 2, p)
third = inv(3)
# lines as coefficient vectors (x, y, z), same order as g26_lattice.txt
L = [(0, 0, 1)]
L += [(0, 1, -a) for a in (0, 1, w - 1, w)]
L += [(1, 0, -b) for b in (0, 1, w, 1 - w)]
L += [(1, 1, -c) for c in (0, 1, w, w + 1)]
L += [(1, w, -d) for d in (0, 1, w - 1, w)]
L += [(1, 1 - w, -1), (1, (w + 1) * third, -(w + 1) * third), (1, w - 1, 0), (1, w + 1, -w)]
L = [tuple(c % p for c in v) for v in L]
def cross(u, v): return ((u[1]*v[2]-u[2]*v[1]) % p, (u[2]*v[0]-u[0]*v[2]) % p, (u[0]*v[1]-u[1]*v[0]) % p)
def norm(v):
    for c in v:
        if c: s = inv(c); return tuple(x * s % p for x in v)
dele, Hi = tuple(int(x) for x in sys.argv[1].split(',')), int(sys.argv[2])
keep = [i for i in range(21) if i not in dele]
H = L[Hi]
groups = {}
for i in keep:
    if i == Hi: continue
    groups.setdefault(norm(cross(H, L[i])), []).append(i)
pts = list(groups.items())
print("points on H:", [(len(g)) for _, g in pts])
# parametrise H by two points e, f on it
e = norm(cross(H, (1, 0, 0))) if any(cross(H, (1, 0, 0))) else norm(cross(H, (0, 1, 0)))
f = norm(cross(H, (0, 0, 1))) if norm(cross(H, (0, 0, 1))) != e else norm(cross(H, (0, 1, 0)))
def coord(P):
    # P = s e + t f ; solve via 2x2 minors; return s/t or None for t = 0 (infinity)
    for (i, j) in ((0, 1), (0, 2), (1, 2)):
        det = (e[i]*f[j] - e[j]*f[i]) % p
        if det:
            s = (P[i]*f[j] - P[j]*f[i]) * inv(det) % p
            t = (e[i]*P[j] - e[j]*P[i]) * inv(det) % p
            return None if t == 0 else s * inv(t) % p
def fmap(t, t0, t1, t2):
    # Moebius map with t0 -> inf, t1 -> 0, t2 -> 1 ; None = infinity
    if t is None: return (t2 - t0) * inv(t2 - t1) % p
    if t0 is None: return (t - t1) * inv(t2 - t1) % p
    if t1 is None: return (t2 - t0) * inv(t - t0) % p
    if t2 is None: return (t - t1) * inv(t - t0) % p
    return (t - t1) * (t2 - t0) * inv((t - t0) * (t2 - t1)) % p
sizes = sorted([len(g) for _, g in pts], reverse=True)
order = [sizes[1], sizes[2], sizes[3], sizes[4], sizes[0]]   # census order: s4 = largest
print("census sizes order:", order)
taus = [coord(P) for P, _ in pts]
results = set()
for perm in permutations(range(5)):
    if [len(pts[k][1]) for k in perm] != order: continue
    t0, t1, t2, t3, t4 = [taus[k] for k in perm]
    l1, l2 = fmap(t3, t0, t1, t2), fmap(t4, t0, t1, t2)
    results.add((l1, l2))
print("lambda pairs:", sorted(results))
