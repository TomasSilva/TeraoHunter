# Exact characteristic-0 Terao check for a connecting-line lattice A(P) found over F_q:
# 1. read off the lattice (collinear point triples, concurrent line sets) from the F_q example;
# 2. build its full realisation space over QQ in Macaulay2 (4 points of P fixed as the frame, the
#    coordinate zero-pattern forced by the lattice), decompose, drop degenerate components;
# 3. on every component decide freeness exactly: free <=> D_0(A)_{d1-1} = 0 (a non-free member of a
#    lattice with chi = (t-1)(t-d1)(t-d2) has mdr <= d1 - 1), including the non-free locus inside
#    positive-dimensional components.
import sys, json, itertools, subprocess, time

q = int(sys.argv[1]); P = json.loads(sys.argv[2]); d1 = int(sys.argv[3]); tag = sys.argv[4]
pts = [v for v in itertools.product(range(q), repeat=3) if any(v) and next(c for c in v if c) == 1]
X = [pts[i] for i in P]
k = len(X)


def cross(u, v):
    return ((u[1]*v[2]-u[2]*v[1]) % q, (u[2]*v[0]-u[0]*v[2]) % q, (u[0]*v[1]-u[1]*v[0]) % q)


def dot(u, v):
    return sum(a*b for a, b in zip(u, v)) % q


def norm(v):
    c = next(x for x in v if x)
    inv = pow(c, q-2, q)
    return tuple((x*inv) % q for x in v)


# lines = maximal collinear subsets of P with >= 2 points
lines = {}
for a, b in itertools.combinations(range(k), 2):
    l = norm(cross(X[a], X[b]))
    lines.setdefault(l, set()).update((a, b))
Lk = list(lines.keys())
Lset = [frozenset(lines[l]) for l in Lk]
n = len(Lk)
# intersection points: group lines by their meeting point
flats = {}
for i, j in itertools.combinations(range(n), 2):
    p = norm(cross(Lk[i], Lk[j]))
    flats.setdefault(p, set()).update((i, j))
inP = {norm(x) for x in X}
conc = [sorted(F) for p, F in flats.items() if len(F) >= 3 and p not in inP]
coll = [sorted(S) for S in Lset if len(S) >= 3]
print('k=%d points, %d lines, %d concurrency flats off P, %d collinear point sets' % (k, n, len(conc), len(coll)), flush=True)
# frame: 4 points of P with no 3 collinear (in the lattice)
def collinear3(a, b, c):
    return any({a, b, c} <= S for S in Lset)
frame = next(f for f in itertools.combinations(range(k), 4) if not any(collinear3(*t) for t in itertools.combinations(f, 3)))
e = list(frame)
others = [i for i in range(k) if i not in e]
# zero pattern of a non-frame point in the frame coordinates: x=0 iff on line e2e3, y=0 iff on e1e3, z=0 iff on e1e2
def on_line(i, a, b):
    return any({i, a, b} <= S for S in Lset)
coords, var = {}, []
coords[e[0]] = ['1', '0', '0']; coords[e[1]] = ['0', '1', '0']; coords[e[2]] = ['0', '0', '1']; coords[e[3]] = ['1', '1', '1']
for i in others:
    zero = [on_line(i, e[1], e[2]), on_line(i, e[0], e[2]), on_line(i, e[0], e[1])]
    c = []
    last = max(j for j in range(3) if not zero[j])
    for j in range(3):
        if zero[j]:
            c.append('0')
        elif j == last:
            c.append('1')
        else:
            v = 'u_%d' % len(var); var.append(v); c.append(v)
    coords[i] = c
print('frame', e, 'unknowns', len(var), flush=True)
vec = lambda i: 'matrix{{%s}}' % ','.join(coords[i])
line_of = lambda S: sorted(S)[:2]          # two points spanning the line
# M2 script
m2 = []
m2.append('R = QQ[%s];' % ','.join(var) if var else 'R = QQ[u_0];')
m2.append('cr = (a,b) -> matrix{{a_(0,1)*b_(0,2)-a_(0,2)*b_(0,1), a_(0,2)*b_(0,0)-a_(0,0)*b_(0,2), a_(0,0)*b_(0,1)-a_(0,1)*b_(0,0)}};')
for i in range(k):
    m2.append('X%d = sub(%s, R);' % (i, vec(i)))
for li, S in enumerate(Lset):
    a, b = line_of(S)
    m2.append('L%d = cr(X%d, X%d);' % (li, a, b))
eqs = []
for S in coll:
    S = sorted(S)
    for c in S[2:]:
        eqs.append('det(X%d || X%d || X%d)' % (S[0], S[1], c))
for F in conc:
    for c in F[2:]:
        eqs.append('det(L%d || L%d || L%d)' % (F[0], F[1], c))
m2.append('I = ideal(%s);' % (','.join(eqs) if eqs else '0_R'))
# degeneracy factors: non-collinear point triples, pairwise point-disjoint non-concurrent line triples
bad = []
for t in itertools.combinations(range(k), 3):
    if not collinear3(*t):
        bad.append('det(X%d || X%d || X%d)' % t)
concset = {frozenset(F) for F in conc}
for t in itertools.combinations(range(n), 3):
    A_, B_, C_ = (Lset[i] for i in t)
    if (A_ & B_) or (A_ & C_) or (B_ & C_):
        continue
    if any(set(t) <= set(F) for F in conc):
        continue
    bad.append('det(L%d || L%d || L%d)' % t)
m2.append('bad = {%s};' % ','.join(bad))
m2.append('''bad = select(bad, f -> f != 0);
if any(bad, f -> f == 0) then print "WARNING zero degeneracy factor";
-- remove degenerate components first: saturate by each degeneracy factor in turn (cheap factors first)
G = I;
bad = sort(bad, f -> (first degree f, size f));
for f in bad do (if G != ideal(1_R) then G = saturate(G, f));
if G == ideal(1_R) then (print "components: 0, nondegenerate: 0  (NOT REALIZABLE in characteristic 0)"; exit 0);
cs = if G == 0 then {ideal(0_R)} else decompose G;
good = select(cs, c -> all(bad, f -> f % c != 0));
print("components: " | toString(#cs) | ", nondegenerate: " | toString(#good) | ", dims " | toString apply(good, dim));''')
# derivation matrix of degree d1-1 with theta(alpha_0) = 0 where line 0 is chosen as z=0 if available
kdeg = d1 - 1
z_line = next((li for li, S in enumerate(Lset) if {e[0], e[1]} <= S), None)
m2.append('S = R[x,y,z]; T = R[s,t]; kdeg = %d;' % kdeg)
m2.append('mons = flatten entries basis(kdeg, S);')
rows = []
for li, S in enumerate(Lset):
    if li == z_line:
        continue
    a, b = line_of(S)
    rows.append((li, a, b))
m2.append('lineData = {%s};' % ','.join('{L%d,X%d,X%d}' % r for r in rows))
m2.append("""condCols = for slot from 0 to 1 list for mu in mons list (
   flatten for D in lineData list (
      al := D#0; P1 := D#1; P2 := D#2;
      sb := map(T, S, {s*P1_(0,0)+t*P2_(0,0), s*P1_(0,1)+t*P2_(0,1), s*P1_(0,2)+t*P2_(0,2)});
      F := al_(0,slot) * sb(mu);
      apply(kdeg+1, j -> lift(coefficient(s^j*t^(kdeg-j), F), R))));
M = transpose matrix flatten condCols;
print("derivation matrix " | toString(numrows M) | " x " | toString(numcols M) | " in degree " | toString kdeg);
Bprod = product bad;
for c in good do (
   A := R/c;
   if dim c == 0 then (
      F := toField A;
      r := rank sub(M, F);
      print("COMPONENT dim 0, degree " | toString(degree c) | ": rank " | toString r | " / " | toString(numcols M) | " -> " | (if r == numcols M then "FREE" else "NOT FREE"));
   ) else (
      r := rank sub(M, frac A);
      an := ann coker transpose sub(M, A);
      ok := (sub(Bprod, A) % radical an) == 0;
      print("COMPONENT dim " | toString(dim c) | ": generic rank " | toString r | " / " | toString(numcols M) | (if r == numcols M then " (generic member FREE)" else " (generic member NOT free)") | "; non-free locus inside degeneracy locus: " | toString ok);
   ));
""")
open('/Users/agent/Projects/TeraoHunter/work/exact_tc_%s.m2' % tag, 'w').write('\n'.join(m2) + '\n')
t0 = time.time()
out = subprocess.run(['/opt/homebrew/bin/M2', '--script', '/Users/agent/Projects/TeraoHunter/work/exact_tc_%s.m2' % tag], capture_output=True, text=True)
print(out.stdout.strip(), out.stderr[-1500:].strip(), '(%.0fs)' % (time.time() - t0), flush=True)
