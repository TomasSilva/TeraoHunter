# Exhaustive Macaulay2 check of Lemma 1: every covering class of PG(3,3) and of PG(4,2) is free with
# exponents (1, q, ..., q^(l-2), |A| - [l-1]_q).  One M2 session per field (batched).
import sys, pickle, itertools, subprocess, time
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter/work')
import numpy as np
import pynauty
from covering_check import proj_points, lines_of


def covering_classes(q, l, reps):
    pts = proj_points(q, l)
    lines = [np.array(sorted(pts.index(p) for p in L)) for L in lines_of(q, l)]
    out = []
    for S in reps:
        mask = np.zeros(len(pts), bool)
        mask[list(S)] = True
        if len(S) > sum(q ** i for i in range(l - 1)) and all(mask[L].any() for L in lines):
            out.append([pts[i] for i in S])
    return out


def orbit_reps(q, l):
    pts = proj_points(q, l)
    NP = len(pts)
    P = np.array(pts)
    inc = (P @ P.T) % q == 0
    adj = {NP + i: list(np.nonzero(inc[i])[0]) for i in range(NP)}

    def cert(S):
        g = pynauty.Graph(2 * NP, directed=False, adjacency_dict=adj,
                          vertex_coloring=[set(range(NP)), set(NP + i for i in S), set(NP + i for i in range(NP) if i not in S)])
        return pynauty.certificate(g)
    level, reps = {cert(()): ()}, []
    for k in range(1, NP + 1):
        nxt = {}
        for S in level.values():
            for j in range(NP):
                if j not in S:
                    T = tuple(sorted(S + (j,)))
                    c = cert(T)
                    if c not in nxt:
                        nxt[c] = T
        level = nxt
        reps.extend(nxt.values())
    return reps


def m2_batch(q, l, arrs, tag):
    names = ['x_%d' % (i + 1) for i in range(l)]
    lines = ['needsPackage "HyperplaneArrangements";', 'R = ZZ/%d[%s];' % (q, ','.join(names))]
    base = sum(q ** i for i in range(l - 1))
    for k, A in enumerate(arrs):
        forms = ['+'.join('%d*%s' % (c, names[i]) for i, c in enumerate(v) if c) for v in A]
        pred = sorted([q ** i for i in range(l - 1)] + [len(A) - base])
        lines.append('A = arrangement {%s}; M = prune image der A; '
                     'print("R %d " | toString(isFreeModule M) | " " | toString(sort flatten degrees M) | " %s");'
                     % (','.join(forms), k, '{' + ', '.join(map(str, pred)) + '}'))
    path = '/Users/agent/Projects/TeraoHunter/work/covering_m2_all_%s.m2' % tag
    open(path, 'w').write('\n'.join(lines) + '\n')
    out = subprocess.run(['/opt/homebrew/bin/M2', '--script', path], capture_output=True, text=True)
    ok = bad = 0
    for line in out.stdout.splitlines():
        if line.startswith('R '):
            _, k, fr, rest = line.split(' ', 3)
            got, pred = rest.rsplit(' {', 1)
            if fr == 'true' and got.strip() == '{' + pred:
                ok += 1
            else:
                bad += 1
                print('MISMATCH', tag, k, line, flush=True)
    return ok, bad, out.stderr[-500:]


t0 = time.time()
reps33 = pickle.load(open('/Users/agent/Projects/TeraoHunter/work/free_df_par_3_4_reps.pkl', 'rb'))
c33 = covering_classes(3, 4, reps33)
print('PG(3,3): %d essential covering classes' % len(c33), flush=True)
ok, bad, err = m2_batch(3, 4, c33, 'pg33')
print('PG(3,3): %d free with predicted exponents, %d mismatches (%.0fs) %s' % (ok, bad, time.time() - t0, err if bad or ok < len(c33) else ''), flush=True)
t0 = time.time()
reps42 = orbit_reps(2, 5)
c42 = covering_classes(2, 5, reps42)
print('PG(4,2): %d orbit classes, %d essential covering classes' % (len(reps42), len(c42)), flush=True)
ok, bad, err = m2_batch(2, 5, c42, 'pg42')
print('PG(4,2): %d free with predicted exponents, %d mismatches (%.0fs) %s' % (ok, bad, time.time() - t0, err if bad or ok < len(c42) else ''), flush=True)
print('DONE', flush=True)
