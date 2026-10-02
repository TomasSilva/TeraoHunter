# Enumerate B' = {P_1..P_r through X1} + triangle {L1, L2, L3} (not through X1), with each P_i through
# at most one vertex of the triangle (distinct vertices).  L0 through X1, locally heavy (m0 >= M1).
# Free (Abe-Kuhne Thm 1.3) <=> b2_off == e1*e2 (Euler-Ziegler exps on L0: points X1 (M1), L0^Lj (m_j)).
# Want: balanced generic exps (M even, max < M/2) and a jumpable pattern (from rank-2 jump table).
import itertools, ast, re, sys
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
from hunter.multi import wakamiko

jumpable = set()
for line in open('/Users/agent/Projects/TeraoHunter/work/rank2_jump_table.out'):
    mm = re.match(r"^(\([^)]*\)) generic (\([^)]*\)) jumps", line)
    if mm:
        jumpable.add(tuple(sorted(ast.literal_eval(mm.group(1)))))
print('jumpable 4-point patterns (max mult 5):', len(jumpable))

def W(*ms):
    e = wakamiko(list(ms))
    return e[0] * e[1]

hits = []
MAXP, MAXL, MAXR = 6, 6, 4
for r in range(1, MAXR + 1):
    # assignment of P_i to vertices: each P_i -> one of {None, 'V12', 'V13', 'V23'}, distinct vertices
    for assign in itertools.product([None, 0, 1, 2], repeat=r):
        vs = [v for v in assign if v is not None]
        if len(vs) != len(set(vs)):
            continue
        if list(assign) != sorted(assign, key=lambda v: (v is None, v if v is not None else 0)):
            continue  # canonical order: assigned first, sorted
        for ps in itertools.product(range(1, MAXP + 1), repeat=r):
            M1 = sum(ps)
            for (a, b, c) in itertools.product(range(1, MAXL + 1), repeat=3):
                ms = (a, b, c)
                # vertex V12 = L1^L2 (opp. L3), V13 = L1^L3 (opp L2), V23 = L2^L3 (opp L1)
                vert = {0: (0, 1, 2), 1: (0, 2, 1), 2: (1, 2, 0)}   # (j, k, opposite l)
                b2 = 0
                covered = set()
                for pi, v in zip(ps, assign):
                    if v is None:
                        b2 += pi * (a + b + c)
                    else:
                        j, k, l = vert[v]
                        b2 += W(pi, ms[j], ms[k]) + pi * ms[l]
                        covered.add(v)
                for v in range(3):
                    if v not in covered:
                        j, k, l = vert[v]
                        b2 += ms[j] * ms[k]
                pat = (M1, a, b, c)
                M = sum(pat)
                if M % 2 or 2 * max(pat) >= M:
                    continue
                if b2 == (M // 2) ** 2:
                    jp = tuple(sorted(pat)) in jumpable
                    hits.append((r, assign, ps, ms, pat, jp))
                    if jp:
                        print('HIT', 'r', r, 'assign', assign, 'p', ps, 'L', ms, 'L0pattern', pat, 'JUMPABLE', flush=True)
print('total balanced-free hits', len(hits), 'jumpable', sum(1 for h in hits if h[-1]))
from collections import Counter
print(Counter((h[0], h[1]) for h in hits).most_common(20))
for h in hits[:30]:
    print(h)
