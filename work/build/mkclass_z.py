# write a class file (incidences) from an integer 4-pencil configuration: sizes, lam, A, B, C, D
import sys
from fractions import Fraction as Fr
sys.path.insert(0, '.')
from pencil_post import lines_of, lattice
sizes = tuple(int(x) for x in sys.argv[1].split(',')); lam = Fr(sys.argv[2])
A = [Fr(x) for x in sys.argv[3].split(',')]; B = [Fr(x) for x in sys.argv[4].split(',')]
C = [Fr(x) for x in sys.argv[5].split(',')]; D = [Fr(x) for x in sys.argv[6].split(',')]
L = lines_of(sizes, lam, A, B, C, D); pts = lattice(L)
print("sizes", *sizes)
for S in pts:
    if len(S) >= 3: print("point", *sorted(S))
