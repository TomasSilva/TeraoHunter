import sys, json, glob, random
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter'); sys.path.insert(0, '/Users/agent/Projects/TeraoHunter/work')
import numpy as np
from hunter.core import Lattice, freeness, rank_mod
from hunter.lift import LiftProblem
src = open('/Users/agent/Projects/TeraoHunter/work/family_search2.py').read()
ns = {}
exec(src.split('def random_recipe')[0].replace("seed = int(sys.argv[1]); rng = random.Random(seed)", "")
        .replace("out = open(f'/Users/agent/Projects/TeraoHunter/work/family2_{seed}.jsonl', 'w')", ""), ns)
replay2, p = ns['replay2'], ns['p']
def run(rec, tries, tag, out):
    recipe, line, m = rec['recipe'], rec['line'], rec['m']
    specials = [x[0] for x in rec['nonfree_u']]
    rng = random.Random(abs(hash(tag)) & 0xffff)
    ug = next(u for u in range(7, p) if u not in specials and replay2(recipe, line, u) is not None)
    P0 = LiftProblem(replay2(recipe, line, ug), m, p)
    seen, stats = set(), {}
    for trial in range(tries):
        pat = P0.greedy(rng)
        if pat is None: continue
        key = tuple(sorted((X, tuple(t)) for X, t in pat))
        if key in seen: continue
        seen.add(key)
        dim0 = len(P0.planes) - rank_mod(P0.system(pat), p)
        c0 = P0.solve(pat, random.Random(1)); A0 = P0.arrangement(c0); L0 = Lattice(A0, p)
        if freeness(A0, p, L0)['status'] != 'free':
            stats[('generic_not_free', dim0)] = stats.get(('generic_not_free', dim0), 0) + 1
            continue
        sig0 = L0.signature()
        for us in specials:
            Ps = LiftProblem(replay2(recipe, line, us), m, p)
            if Ps.points != P0.points: continue
            dims = len(Ps.planes) - rank_mod(Ps.system(pat), p)
            verdict = 'lattice_changes'
            for s in range(6):
                cs = Ps.solve(pat, random.Random(s))
                if cs is None or not Ps._distinct(cs): continue
                try:
                    As = Ps.arrangement(cs); Ls = Lattice(As, p)
                except ValueError:
                    continue
                if Ls.signature() == sig0:
                    fs = freeness(As, p, Ls)['status']
                    verdict = 'SAME_LATTICE_' + fs
                    if fs != 'free':
                        out.write(json.dumps({'tag': tag, 'rec': rec, 'pattern': pat, 'us': us, 'seed': s, 'ug': ug}) + '\n'); out.flush()
                        print('!!! CANDIDATE COUNTEREXAMPLE', tag, us, flush=True)
                    break
            stats[(verdict, dim0, dims)] = stats.get((verdict, dim0, dims), 0) + 1
    return len(seen), stats
if __name__ == '__main__':
    fn, idx, tries = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    hits = [json.loads(l[4:]) for l in open(fn) if l.startswith('HIT')]
    rec = hits[idx]
    out = open(fn.replace('.out', f'_cand{idx}.jsonl'), 'w')
    n, stats = run(rec, tries, f'{fn}#{idx}', out)
    print('hit', fn, idx, 'n', rec['n'], 'm', rec['m'], 'special u', [x[0] for x in rec['nonfree_u']], 'patterns', n, stats, flush=True)
