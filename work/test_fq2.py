import sys, random, itertools
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter'); sys.path.insert(0, '/Users/agent/Projects/TeraoHunter/work')
from hunter.fq import GF, LineArr
from test_fq import m2_free_fq
rng = random.Random(5)
for q in [4, 8, 9]:
    F = GF(q); pts = F.points(3)
    allL = list(pts)
    O = (0, 0, 1)
    through_O = [L for L in allL if F.dot(L, O) == 0]
    cases = {'PG(2,q)': allL, 'PG minus a line': allL[1:], 'lines not through O': [L for L in allL if F.dot(L, O) != 0],
             'PG minus pencil + 1': [L for L in allL if F.dot(L, O) != 0] + [through_O[0]]}
    # random deletions from PG(2,q) keeping chi split
    for t in range(40):
        k = rng.randint(3, 12)
        S = rng.sample(allL, len(allL) - k)
        try:
            A = LineArr(F, S)
        except ValueError:
            continue
        if A.exponents() is not None:
            cases[f'PG minus {k} random #{t}'] = S
    mism = 0; n = 0
    for name, S in list(cases.items())[:10]:
        A = LineArr(F, S)
        ours = A.free()
        theirs = m2_free_fq(F, S)
        n += 1; mism += ours != theirs
        print(q, name, len(S), A.exponents(), 'ours', ours, 'M2', theirs, '' if ours == theirs else 'MISMATCH', flush=True)
    print('GF(%d): %d tests, %d mismatches' % (q, n, mism), flush=True)
