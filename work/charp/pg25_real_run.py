#!/usr/bin/env python3
"""Convert pg25_hits.txt (one representative per invariant key) to the JSON format of pg2q_real.py, generate the M2
realization-scheme scripts over F_5-bar and run them with a per-lattice timeout, summarizing the outcome."""
import json, subprocess, sys, re, os
hits = [l for l in open('pg25_hits.txt') if ';' in l]
cands = []
for l in hits:
    key, idx = l.split(';')
    n = int(re.match(r'n(\d+)', key).group(1)); d1 = int(re.search(r' d(\d+) ', key).group(1))
    S = [int(x) for x in idx.split()]
    cands.append(dict(n=n, d1=d1, profile=key.split('|')[0], S=S))
cands.sort(key=lambda c: (c['n'], c['d1']))
json.dump(cands, open('pg2_5_cands.json','w'))
print(len(cands), 'lattices', file=sys.stderr)
subprocess.run([sys.executable, '-I', 'pg2q_real.py', '5', 'pg2_5_cands.json'], check=True, stdout=subprocess.DEVNULL)
summary = []
for i, c in enumerate(cands):
    fn = f'real_pg5_{i}.m2'
    try:
        r = subprocess.run(['M2', '--script', fn], capture_output=True, text=True, timeout=int(sys.argv[1]) if len(sys.argv) > 1 else 900)
        out = r.stdout + r.stderr
    except subprocess.TimeoutExpired:
        out = 'TIMEOUT'
    lines = [l for l in out.splitlines() if 'back trace' not in l]
    flag = 'FREE' if 'FREE' in out else ('POSITIVE' if 'POSITIVE-DIMENSIONAL' in out else ('TIMEOUT' if 'TIMEOUT' in out else ('ERROR' if 'error' in out.lower() else 'rigid/nonfree')))
    summary.append((i, c['n'], c['d1'], flag))
    with open('pg25_real_results.txt','a') as f:
        f.write(f'### lattice {i} n={c["n"]} class=({c["d1"]},{c["n"]-1-c["d1"]}) {c["profile"]} -> {flag}\n' + '\n'.join(lines) + '\n')
    print(f'lattice {i}: n={c["n"]} class=({c["d1"]},{c["n"]-1-c["d1"]}) -> {flag}', file=sys.stderr, flush=True)
from collections import Counter
print(Counter((n, d1, flag) for _, n, d1, flag in summary), file=sys.stderr)
