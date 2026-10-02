# Independent freeness check of the F_2 rank-6 example (50 hyperplanes) with Macaulay2:
# truncated syzygies give D(A) up to degree 19; Saito's criterion is then applied to
# E, theta_F, theta_F^2, theta_F^3, theta_F^4 (sum x_i^(2^j) d_i) and a degree-19 element.
import json, subprocess, tempfile, os, time
d = json.load(open('/Users/agent/Projects/TeraoHunter/work/f2_rank6_cap_nice.json'))
N = d['normals']
forms = ['+'.join('x_%d' % (i + 1) for i in range(6) if v[i]) for v in N]
script = r'''
R = ZZ/2[x_1..x_6];
L = {%s};
n = #L;
M = matrix apply(n, H -> join(apply(6, i -> coefficient(R_i, L#H)), apply(n, j -> if j == H then L#H else 0_R)));
Phi = map(R^n, , M);
t0 = cpuTime();
S = syz(Phi, DegreeLimit => 19);
print("SYZ done " | toString(cpuTime() - t0) | "s, generators: " | toString numcols S);
G = submatrix(S, {0..5}, );
G = mingens image G;
print("GENDEGS " | toString sort flatten degrees source G);
Q = product L;
fro = apply(5, j -> transpose matrix{apply(6, i -> R_i^(2^j))});
cand = select(numcols G, c -> (degrees source G)#c == {19});
found = false;
for c in cand do (
  T = fold((a, b) -> a | b, fro) | G_{c};
  dt = det T;
  if dt != 0 then (
     print("DET degree " | toString first degree dt | ", equals Q_A: " | toString(dt == Q));
     found = true;
     break;
  );
);
if not found then print("NO degree-19 element with nonzero determinant");
''' % ','.join(forms)
with open('/Users/agent/Projects/TeraoHunter/work/f2_rank6_saito.m2', 'w') as fh:
    fh.write(script)
t0 = time.time()
out = subprocess.run(['/opt/homebrew/bin/M2', '--script', '/Users/agent/Projects/TeraoHunter/work/f2_rank6_saito.m2'],
                     capture_output=True, text=True)
print(out.stdout, out.stderr[-2000:], '(%.0fs)' % (time.time() - t0), flush=True)
