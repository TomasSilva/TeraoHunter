# Three-pencil lattices at 18 lines (TC_FOLIATION.md §11)

Lines `x = a·y` (a ∈ A, through [0:0:1]), `z = c·x` (c ∈ B, through [0:1:0]), `y = b·z` (b ∈ C,
through [1:0:0]) and coordinate lines K ⊂ {x, y, z}. Triple points are exactly the triples with
abc = 1, all other intersections are double. The realization space of such a lattice is the
group of solutions of the binomial system a_i c_j b_k = 1 modulo the diagonal torus; its
components are the torsion of Z^m/(relations), computed by Smith normal form.

- `seed18.m2` / `seed18.out`: the arrangement yz(x⁶−y⁶)(z⁶−x⁶)∏_{b∈μ₆∖{1,ζ}}(y−bz) (= A(6,1,3)
  minus x, y−z, y−ζz) and its two variants, over F_10009: mdr 7, Bourbaki degree 2, J/Q = product
  of the three deleted lines.
- `real18.m2` / `real18.out`: realization space of its lattice over ℚ (6 points, only the
  Galois pair ζ₆^{±1} non-degenerate: rigid).
- `pencil_enum.py [P1 ...] [--wide]`: enumerates (A,B,C) ⊂ Z/N, Z, Z×Z/N, Z² with the required
  numbers of triple points for every multiplicity pattern of the class (8,9), dedupes by a
  signature, computes the components of the realization space and writes a generic point of
  every non-degenerate component over F_p into `pencil_entries*.json` and an M2 batch
  `pencil_mdr*.m2`.
- `pencil_iso.py entries.json`: exact isomorphism classes of the triple systems found.
- `pencil_jump.py entries.json [sigid ...]`: for every lattice, the normalized components
  (a_0 = c_0 = 1), mdr at a random point of each, and for positive-dimensional components the
  non-free locus by Yoshinaga's criterion (16-minors of the 17×16 matrix of the degree-7
  derivations of the Ziegler multirestriction on a coordinate line, saturated by the
  degeneration binomials). `jump_P1_3_QQ.m2/.out`: the one-parameter family, over ℚ.
- `pencil_mdr_P1.out`, `pencil_P1_summary.txt`, `pencil_entries_P1.json`: results for the
  pattern (5,5,5; xyz; 19): six lattices, all free.
