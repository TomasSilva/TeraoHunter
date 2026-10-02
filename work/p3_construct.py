# Direct P^3 search: construction sequences of planes (through existing points/lines,
# with a few free parameters), cheap combinatorial filters first:
#   chi splits over Z>=0  ->  NOT divisionally free  ->  locally free  ->  free.
# Free non-DF members are then scanned along their parameters for non-free points
# with the same labelled lattice.
import sys, random, time, json
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import Lattice, freeness, local_freeness_l4, nullspace_mod, _minors3
from hunter.divisional import df_witnesses

p = 32003
seed = int(sys.argv[1]); rng = random.Random(seed)
TLIM = float(sys.argv[2])
out = open(f'/Users/agent/Projects/TeraoHunter/work/p3c_{seed}.jsonl', 'w')


def norm(v):
    v = [int(x) % p for x in v]
    for c in v:
        if c:
            inv = pow(c, p - 2, p)
            return tuple(x * inv % p for x in v)
    return None


def plane_through(points):
    """Normal of the plane through the given points (rows), if unique-ish
    (returns a basis of the solution space)."""
    M = np.array(points, dtype=np.int64) % p
    return nullspace_mod(M, p)


def replay(recipe, params):
    planes = [norm(v) for v in ([1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[1,1,1,1])]
    for step in recipe:
        kind = step[0]
        if kind == 'ppp':        # through three points (each = triple of plane indices)
            pts = [tuple(_minors3(np.array(planes[a]), np.array(planes[b]), np.array(planes[c]), p))
                   for (a, b, c) in step[1]]
            K = plane_through(pts)
            if K.shape[0] != 1: return None
            v = norm(K[0])
        elif kind == 'lp':       # through a line (pair of planes) and a point
            (a, b), (c, d, e) = step[1], step[2]
            pt = _minors3(np.array(planes[c]), np.array(planes[d]), np.array(planes[e]), p)
            # plane in pencil  s*P_a + P_b  through pt
            A, B = np.array(planes[a]), np.array(planes[b])
            fa, fb = int(A @ pt % p), int(B @ pt % p)
            if fa == 0 and fb == 0: return None
            v = norm((fb * A - fa * B) % p) if True else None
        elif kind == 'l1':       # through a line, free parameter
            (a, b), w = step[1], step[2]
            v = norm((np.array(planes[a]) + params[w] * np.array(planes[b])) % p)
        elif kind == 'pp1':      # through two points, free parameter (pencil around their line)
            pts = [tuple(_minors3(np.array(planes[a]), np.array(planes[b]), np.array(planes[c]), p))
                   for (a, b, c) in step[1]]
            K = plane_through(pts)
            if K.shape[0] != 2: return None
            v = norm((K[0] + params[step[2]] * K[1]) % p)
        if v is None or v in planes:
            return None
        planes.append(v)
    return np.array(planes, dtype=np.int64)


def flats_of(N):
    L = Lattice(N, p)
    lines = [sorted(S) for S in L.flats[2]]
    pts = [sorted(S) for S in L.flats[3]]
    return L, lines, pts


def random_recipe(n, npar):
    params = [rng.randrange(2, p) for _ in range(npar)]
    recipe = []
    N = replay([], params)
    par_steps = set(rng.sample(range(n - 5), npar))
    for s in range(n - 5):
        L, lines, pts = flats_of(N)
        rich_pts = [P for P in pts if len(P) >= 3]
        rich_lines = [l for l in lines if len(l) >= 2]
        if s in par_steps:
            if rng.random() < 0.5:
                l = rng.choice(rich_lines)
                step = ('l1', (l[0], l[1]), len([x for x in recipe if x[0] in ('l1', 'pp1')]))
            else:
                P1, P2 = rng.sample(rich_pts, 2)
                step = ('pp1', [tuple(P1[:3]), tuple(P2[:3])], len([x for x in recipe if x[0] in ('l1', 'pp1')]))
        else:
            if rng.random() < 0.5:
                P3 = rng.sample(rich_pts, 3)
                step = ('ppp', [tuple(P[:3]) for P in P3])
            else:
                l = rng.choice(rich_lines)
                P = rng.choice(rich_pts)
                step = ('lp', (l[0], l[1]), tuple(P[:3]))
        N2 = replay(recipe + [step], params)
        if N2 is None:
            continue
        recipe.append(step)
        N = N2
    return recipe, params, N


t0 = time.time()
stats = {'built': 0, 'chi': 0, 'nondf': 0, 'locfree': 0, 'free': 0}
while time.time() - t0 < TLIM:
    n = rng.randint(11, 18)
    npar = rng.choice([1, 1, 2])
    try:
        recipe, params, N = random_recipe(n, npar)
    except (ValueError, IndexError):
        continue
    if N.shape[0] < 10:
        continue
    try:
        L = Lattice(N, p)
    except ValueError:
        continue
    stats['built'] += 1
    if L.exponents() is None:
        continue
    stats['chi'] += 1
    if df_witnesses(N, p, L):
        continue
    stats['nondf'] += 1
    if local_freeness_l4(N, p, L):
        continue
    stats['locfree'] += 1
    fr = freeness(N, p, L)
    rec = {'recipe': recipe, 'params': params, 'n': int(N.shape[0]), 'exps': fr['exps'], 'status': fr['status']}
    if fr['status'] == 'free':
        stats['free'] += 1
        out.write(json.dumps(rec) + '\n'); out.flush()
        print('FREE NON-DF', json.dumps(rec), flush=True)
print(stats, time.time() - t0, flush=True)
