"""
Consistency test of Theorem F (TC_GLOBAL.md section 7) on small arrangements.

Families: H = line at infinity with four pencils
  B: x = c (k1 lines), A: y = a (k2), C: x - y = b (k3), D: x - mu*y = d (k4),
i.e. the 4-point vector (k1, k2, k3, k4) in the order of work/rational/jump4.m2.
A, B, C offsets are random small integers, with C lines chosen through grid points; each D line
either passes through one chosen point of the A/B/C configuration (offset d = x_P - mu*y_P, so the
lattice does not force mu) or has a random offset. Such a family lies in the component of the
realisation space that dominates the cross-ratio mu.

Prediction (Theorem F + Lemma 2): if the generic member has chi_0 = (t - d)^2 (so it is free, since
mu is not a jump value), then at the rational jump value mu0 the lattice must change
(the member degenerates), and when it hits a generic point of the limit fibre, beta = 1
(exactly one new triple point).
"""
import random, sys
from fractions import Fraction as Fr
from collections import defaultdict

def cross(l1, l2):
    a = (l1[1]*l2[2] - l1[2]*l2[1], l1[2]*l2[0] - l1[0]*l2[2], l1[0]*l2[1] - l1[1]*l2[0])
    for v in a:
        if v != 0:
            return tuple(x / v for x in a)
    return None  # same line

def lattice(lines):
    """Return (ok, points) where points maps a projective point to the frozenset of line indices."""
    n = len(lines)
    pts = defaultdict(set)
    for i in range(n):
        for j in range(i + 1, n):
            p = cross(lines[i], lines[j])
            if p is None:
                return False, None
            pts[p].update((i, j))
    return True, {p: frozenset(s) for p, s in pts.items()}

def b2(pts):
    return sum(len(s) - 1 for s in pts.values())

def structure(pts):
    return frozenset(s for s in pts.values() if len(s) >= 3)

def build(k, mu, rng_state):
    rng = random.Random(rng_state)
    k1, k2, k3, k4 = k
    R = 6
    cs = rng.sample(range(-R, R + 1), k1)          # B: x = c
    as_ = rng.sample(range(-R, R + 1), k2)         # A: y = a
    bs = set()
    tries = 0
    while len(bs) < k3 and tries < 100:
        tries += 1
        if rng.random() < 0.8:
            bs.add(rng.choice(cs) - rng.choice(as_))
        else:
            bs.add(rng.randint(-2 * R, 2 * R))
    if len(bs) < k3:
        return None
    bs = list(bs)
    base = []   # candidate base points for D lines: grid, A-C, B-C intersections
    for c in cs:
        for a in as_:
            base.append((Fr(c), Fr(a)))
    for a in as_:
        for b in bs:
            base.append((Fr(a + b), Fr(a)))
    for c in cs:
        for b in bs:
            base.append((Fr(c), Fr(c - b)))
    dspec = []
    for _ in range(k4):
        if rng.random() < 0.85:
            dspec.append(('pt', rng.choice(base)))
        else:
            dspec.append(('free', Fr(rng.randint(-3 * R, 3 * R))))
    def lines_at(m):
        L = [(Fr(0), Fr(0), Fr(1))]                                  # H
        L += [(Fr(1), Fr(0), Fr(-c)) for c in cs]
        L += [(Fr(0), Fr(1), Fr(-a)) for a in as_]
        L += [(Fr(1), Fr(-1), Fr(-b)) for b in bs]
        for kind, val in dspec:
            if kind == 'pt':
                x, y = val
                d = x - m * y
            else:
                d = val
            L.append((Fr(1), -m, -d))
        return L
    return lines_at

def run(k, mu0, trials, seed=1):
    K = sum(k); d = K // 2
    target = d * d + 2 * d
    rng = random.Random(seed)
    stats = defaultdict(int)
    seen = set()
    for t in range(trials):
        f = build(k, None, rng.random())
        if f is None:
            continue
        g1, g2 = Fr(rng.randint(5, 50), rng.randint(51, 97)) + 3, Fr(-rng.randint(5, 50), rng.randint(51, 97)) - 5
        ok1, P1 = lattice(f(g1)); ok2, P2 = lattice(f(g2))
        if not (ok1 and ok2) or structure(P1) != structure(P2) or b2(P1) != b2(P2):
            stats['non-generic sample'] += 1
            continue
        if b2(P1) != target:
            stats['chi not (t-d)^2'] += 1
            continue
        key = structure(P1)
        if key in seen:
            stats['repeat'] += 1
            continue
        seen.add(key)
        stats['balanced generic lattices'] += 1
        ok0, P0 = lattice(f(Fr(mu0)))
        if not ok0:
            stats['at mu0: two lines coincide'] += 1
            continue
        if structure(P0) == structure(P1) and b2(P0) == b2(P1):
            stats['at mu0: SAME LATTICE (would contradict Theorem F)'] += 1
            continue
        beta = b2(P1) - b2(P0)
        stats['at mu0: lattice changes, beta = %d' % beta] += 1
    return stats

if __name__ == '__main__':
    cases = [((1, 1, 3, 3), -1), ((1, 2, 2, 3), 2), ((1, 3, 4, 4), -1), ((2, 2, 3, 5), -1),
             ((3, 3, 3, 3), -1), ((1, 3, 3, 5), 2)]
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
    for k, mu0 in cases:
        st = run(k, mu0, trials)
        print(k, 'mu0 =', mu0, dict(sorted(st.items())), flush=True)
