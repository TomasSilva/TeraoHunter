# Run the exact char-0 Terao checker on every distinct non-DF free connecting-line lattice found over
# F_11, F_13, F_17 (one representative per weak combinatorics + exponents).
import json, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
reps = {}
for q in [11, 13, 17]:
    for line in open('/Users/agent/Projects/TeraoHunter/work/connecting_%d.jsonl' % q):
        r = json.loads(line)
        key = (r['n_lines'], tuple(r['exps']), tuple(sorted(r['t'])))
        reps.setdefault(key, (q, r['P'], r['exps'][1]))
jobs = sorted(reps.items())
print(len(jobs), 'distinct lattices', flush=True)


def run(item):
    (n, ex, t), (q, P, d1) = item
    tag = 'n%d_%s_%d' % (n, '_'.join(map(str, ex)), abs(hash(t)) % 100000)
    try:
        out = subprocess.run([sys.executable, '/Users/agent/Projects/TeraoHunter/work/exact_tc_points.py', str(q), json.dumps(P), str(d1), tag],
                             capture_output=True, text=True, timeout=5400).stdout
    except subprocess.TimeoutExpired:
        out = 'TIMEOUT'
    keep = [l for l in out.splitlines() if l.startswith(('components', 'COMPONENT', 'TIMEOUT', 'k=')) or 'rror' in l]
    return '%d lines %s %s (q=%d, k=%d):\n   %s' % (n, ex, dict(sorted({x: t.count(x) for x in set(t)}.items())), q, len(P), '\n   '.join(keep))


with ThreadPoolExecutor(6) as ex:
    for res in ex.map(run, jobs):
        print(res, flush=True)
print('ALLDONE', flush=True)
