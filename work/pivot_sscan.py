import sys, ast, re, time
sys.path.insert(0, '/Users/agent/Projects/TeraoHunter')
import numpy as np
from hunter.core import multi_freeness, Lattice
p = int(sys.argv[2]) if len(sys.argv) > 2 else 1009
cores = {
 'A3': [[1,0,0],[0,1,0],[0,0,1],[1,-1,0],[1,0,-1],[0,1,-1]],
 'delA3': [[1,0,0],[0,1,0],[0,0,1],[1,-1,0],[1,0,-1]],
 'X3_2': [[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,2]],
 'bool': [[1,0,0],[0,1,0],[0,0,1]],
 'gen4': [[1,0,0],[0,1,0],[0,0,1],[1,1,1]],
 'B3sub7': [[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,-1,0],[1,0,1],[1,0,-1]],
}
pat = re.compile(r"^(\S+) X1= (\[[^\]]*\]) m'= (\([^)]*\)) m0= (\d+)")
cands = []
for line in open(sys.argv[1]):
    mm = pat.match(line)
    if not mm: continue
    name, S, mc, m0 = mm.group(1), ast.literal_eval(mm.group(2)), ast.literal_eval(mm.group(3)), int(mm.group(4))
    core = cores[name]
    # multiplicities on L0: point X1 has M1 = sum over lines through X1; other lines meet L0 at distinct points (generic)
    M1 = sum(mc[i] for i in S)
    others = [mc[i] for i in range(len(core)) if i not in S]
    M = M1 + sum(others)
    if max([M1] + others) * 2 > M:
        continue          # dominant point on L0: Euler exponents cannot jump
    cands.append((name, S, mc, m0, M1, others))
print('non-dominant candidates', len(cands), flush=True)
seen = set()
for name, S, mc, m0, M1, others in cands:
    key = (name, tuple(S), mc)
    core = cores[name]
    u, w = core[S[0]], core[S[1]]
    def B(s):
        return np.array(core + [[(u[i] + s * w[i]) % p for i in range(3)]], dtype=np.int64) % p
    sig = Lattice(B(123), p).signature()
    m = list(mc) + [m0]
    t0 = time.time()
    bad = []
    nsame = 0
    for s in range(1, p):
        try:
            L = Lattice(B(s), p)
        except ValueError:
            continue
        if L.signature() != sig:
            continue
        nsame += 1
        st, g = multi_freeness(B(s), m, p)
        if st != 'free':
            bad.append((s, st, tuple(g)))
    print(name, S, mc, 'm0', m0, 'L0 mults', [M1] + others, '| same-lattice s:', nsame,
          '| NONFREE at', bad[:12], ('...%d' % len(bad)) if len(bad) > 12 else '', '(%.0fs)' % (time.time() - t0), flush=True)
