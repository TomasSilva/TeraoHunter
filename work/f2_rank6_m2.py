import sys, json, time
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
from hunter.m2 import m2_free
d = json.load(open('/Users/agent/Projects/TeraoHunter/work/f2_rank6_cap.json'))
t0 = time.time()
res = m2_free(d['normals'], field='ZZ/2', timeout=7000)
print('M2 GF(2) rank 6 cap complement: free =', res.get('free'), 'gens', res.get('gens'), 'betti', res.get('betti'),
      '(%.0fs)' % (time.time() - t0), flush=True)
