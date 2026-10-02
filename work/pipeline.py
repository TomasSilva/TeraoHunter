# Full pipeline for one multiarrangement hit:
#  valid lift patterns at generic t (savings == R on generic solutions, A free),
#  replay at special t*: same labelled lattice?  A_{t*} freeness?
import sys, json, glob, random
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter'); sys.path.insert(0, '/Users/agent/Projects/TeraoHunter/work')
import numpy as np
from hunter.core import Lattice, freeness, local_freeness_l4, rank_mod
from hunter.lift import LiftProblem
from lift_hit_lib import replay, p

def run(rec, tries, tag, out):
    recipe, m = rec['recipe'], rec['m']
    specials = [x[0] for x in rec['nonfree_t']]
    rng = random.Random(hash(tag) & 0xffff)
    P0 = LiftProblem(replay(recipe, 500), m, p)
    seen, stats = set(), {}
    for trial in range(tries):
        pat = P0.greedy(rng)
        if pat is None: continue
        key = tuple(sorted((X, tuple(t)) for X, t in pat))
        if key in seen: continue
        seen.add(key)
        dim0 = len(P0.planes) - rank_mod(P0.system(pat), p)
        c0 = P0.solve(pat, random.Random(1)); A0 = P0.arrangement(c0); L0 = Lattice(A0, p)
        f0 = freeness(A0, p, L0)['status']
        if f0 != 'free':
            stats[('generic_not_free', dim0)] = stats.get(('generic_not_free', dim0), 0) + 1
            continue
        sig0 = L0.signature()
        for ts in specials:
            Ps = LiftProblem(replay(recipe, ts), m, p)
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
                        out.write(json.dumps({'tag': tag, 'rec': rec, 'pattern': pat, 'ts': ts, 'seed': s}) + '\n'); out.flush()
                        print('!!! CANDIDATE COUNTEREXAMPLE', tag, ts, flush=True)
                    break
            stats[(verdict, dim0, dims)] = stats.get((verdict, dim0, dims), 0) + 1
    return len(seen), stats

if __name__ == '__main__':
    hits = []
    for fn in sorted(glob.glob('/Users/agent/Projects/TeraoHunter/work/family_search_*.out')):
        for line in open(fn):
            if line.startswith('HIT'):
                hits.append(json.loads(line[4:]))
    k = int(sys.argv[1]); tries = int(sys.argv[2])
    out = open(f'/Users/agent/Projects/TeraoHunter/work/candidates_{k}.jsonl', 'w')
    n, stats = run(hits[k], tries, f'hit{k}', out)
    print('hit', k, 'n', hits[k]['n'], 'm', hits[k]['m'], 'patterns', n, stats, flush=True)
