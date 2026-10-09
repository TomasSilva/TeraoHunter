# TeraoHunter — research notes (P^3 and above)

## Tooling (validated)
- `hunter/core.py`: arrangements over F_p (numpy + python-flint). Lattice (l=3,4), chi(A,t),
  exponents, local freeness (l=4), Hilbert function of D_H(A) (gauge D(A) = S.E (+) D_H(A)),
  Saito-determinant freeness test; multiarrangements D(A,m) via Taylor conditions,
  minimal generator degrees, Saito test.
- `hunter/m2.py`: Macaulay2 (HyperplaneArrangements) exact check over QQ.
- Cross-checked: Hilbert functions and freeness agree with M2 on A4, B4, generic, random
  arrangements (l = 3, 4) and on several multiarrangements.

## Reduction used (Yoshinaga, l = 4)
A free  <=>  A locally free along H0  and  Ziegler restriction (A^{H0}, m^{H0}) free.
Decone w.r.t. H0: affine planes a_L . x = c_{L,i} in K^3.  (A^{H0}, m) = the normal
directions B (lines of P^2) with multiplicity = number of parallel planes.
A_X for X in P(H0) (a direction v) = H0 + planes containing direction v: an affine line
arrangement (projected along v) + line at infinity; must be free (combinatorial if Terao
holds in rank 3 for that size).  Hence a P^3 counterexample needs a rank-3 MULTIarrangement
family (B_t, m) with constant lattice and t-dependent freeness, which lifts to A_t with
constant lattice.

## Finding 1: multiarrangement freeness is not combinatorial (X3 family)
X3(t) = {x, y, z, x+y, x+z, y+tz}, lattice constant for t != 0, -1 (3 triple points).
Over F_1009, multiplicities in {1..4}^6: generically free only (k,k,k,1,1,1), k=2,3,4,
exps (k+1,k+1,k+1).
- k=2: free for all t.
- k=3: NON-free exactly at t = 1 (gens 3,5,5,5).
- k=4: NON-free exactly at t^2 - t + 1 = 0 (t = -omega) (gens 4,6,6,6).

## Finding 2: the P^3 lift of X3 forces t (monodromy obstruction)
Local freeness at the triple point (x,y,x+y) with multiplicities (k,k,1) forces the single
plane x+y=d through k grid points (T = k triple lines): b_{s(i)} = d - a_i.  Same at the other
two triple points.  Combining: a_i + t a_{pi(i)} = C for a permutation pi of {1..k}.
Fixed points of pi create 6-fold points with localization X3 (non-free) => pi derangement;
distinctness of a_i => -t is a PRIMITIVE j-th root of unity for every cycle length j.
So the lattice of the lift forces t into one Galois orbit (roots of Phi_j(-t)), e.g.
k=3 (pi = 2-cycle+fixed pt is excluded...) — t is never free to move.  Since freeness is
Galois invariant, no counterexample from this family.  Lesson: local freeness at
infinity imposes Ceva-type monodromy conditions on the moduli of B.

## Rank-2 exponent jumps (points 0, inf, 1, lambda on P^1)
(3,3,3,3): (5,7) at harmonic lambda; (2,2,3,3): (4,6) at lambda^2-7lambda+1=0;
(2,2,2,4): (4,6) at lambda^2-lambda+1=0; (1,2,2,3): (3,5) at lambda=2.

## Literature anchors (survey, Oct 2026)
- No counterexample known in any rank/field; no exhaustive rank-4 verification exists.
- DiPasquale–Wakefield (arXiv:1707.03961): X3 family x^n y^n z^n (x-ay)(x+z)(y+z) free iff
  ord(a) does not divide n-1; their Cor 5.5: every rank-4 extension with this Ziegler
  restriction satisfies TC.  (Matches Finding 2: lift forces ord(-t) | n, non-free needs
  ord(-t) | n-1, coprime.)
- Abe–Kuhne (arXiv:1906.02188) Thm 1.3: L0 locally heavy (m0 >= sum of other mults at each
  multiple point on L0) => (A,m) free iff b2_off = e1 e2 (Euler–Ziegler exps on L0).
  Thm 1.6: arrangements with a locally heavy FLAG satisfy TC.
- Abe–Yoshinaga Cor 5.3: A1 free, isomorphic Ziegler restriction + same b2 => A2 free.
  So in a counterexample family EVERY Ziegler restriction must vary.

## Finding 3: rank-2 jump table (4 points, mults <= 5)
Jumps only for even M (|e1-e2| <= #points - 2 = 2).  Generic e1 = min(M//2, M-k+1, M-max).

## Finding 4: pivot construction (heavy L0 rotating through X1) fails with 4 points on L0
Exhaustive over pencil-through-X1 + triangle / concurrent triple: b2_off never equals the
balanced e1e2 for a non-dominant pattern.  Counting heuristic: freeness needs savings
~1/3 of cross pairs, forest-like incidence gives far less.

## Finding 5: cone + translates lift  A = cone(C) u {w} u {alpha_L - c w}
Locally free iff <= 1 heavy line per multiple point (generic offsets); then
A free <=> (C,m) free (validated numerically).  Lattice constant along moduli of C.
But: one heavy line => always free (ATW addition; combinatorial); two heavy lines meeting
at a double point => never free.  Database test (827 free C with moduli, 13.5k m's): 0 free.
Allowing 2 heavy lines per point (forest heavy-pair graph): 2.5k free, all of the pattern
"two heavy lines through one multiple point, exps (1+k, d1, d2)" = combinatorial; one outlier
(two pencils sharing a line) is a local cross-ratio effect (does not lift).

## Finding 6: multiarrangement families with GLOBAL t-dependent freeness are common
Random 1-parameter construction sequences (6-10 lines, F_1009), mults 1-4: 12 genuine hits
(after hardening the randomized Saito test: 2 earlier hits were false negatives).
All of type exps (k,k,k) generically -> (k-1,k+1,k+1,k+1) at special t; local exponents
unchanged (global phenomenon).  Special t: 1, 2, -2, +-1/2, -omega, t^2+7t+1=0 (golden
ratio: t = -phi^{+-4}), others.

## Finding 7: P^3 lifts exist generically, but DEGENERATE exactly at t*  (protection)
hunter/lift.py: greedy search of concurrency patterns realising the savings R_X at every
point of B (=> locally free along H), validated on generic solutions.  For 5 hits, 322
distinct valid patterns: all lifts free for generic t with constant lattice (genuine
1-parameter families of free plane arrangements in P^3), all RIGID (intercepts unique up to
affine gauge), and EVERY one acquires an extra incidence at t* (e.g. 4 planes y, z, x+z, L7
through one affine point -> non-free localization).  So TC survives: the lattice of A "sees"
the special parameter that the lattice of B does not.
Pitfall fixed: greedy patterns with solution dim 5 were artifacts (savings only at a special
solution); now re-validated on fresh generic solutions.
Next: 2-parameter families (N = non-free curve vs D = lift-degeneration curve).

## Finding 8: the protection mechanism is DIVISIONAL FREENESS (Abe, Invent. Math. 2016)
Every valid free lift found (all hits, ~60+ lifts) has planes P whose simple restriction A^P
is free with chi(A^P) | chi(A) (e.g. 9-line restrictions for exps (1,4,4,4)).  DF is
combinatorial (given rank-3 TC for A^P, <= 14 lines), so TC is GUARANTEED for these lattices,
which is why the lattice must change at t*.  Ruled out as explanations: heavy Ziegler
restrictions (none), rigid restrictions (all have moduli >= 1).
=> A P^3 counterexample must be FREE but NOT divisionally free (every plane P: A^P non-free,
   or not DF, or chi(A^P) does not divide chi(A)).  hunter/divisional.py implements the test.
DF statistics (work/df_stats_*.out): 351 generic-free multiarrangements (6-10 lines), 2361
valid free P^3 lifts -> 0 non-DF.  The lift construction lands in DF every time.

## Known free non-DF arrangements (Abe, arXiv:1502.07520)
G31 IS divisionally free (Ex. 4.6).  Non-DF free: A^0_l(r) = G(r,r,l), prod (x_i^r - x_j^r),
r >= 3, l >= 3 (Thm 5.6) -- rigid.  Abe's "valley" picture: A^k_l(r) DF iff k != 0.

## Lead: bias (gain-graphic) arrangements  B_q(G) = {x_i} u {x_i - q^g x_j : [i,j,g] in G}
(arXiv:2403.08304).  Lattice depends only on additive relations of integer gains -> constant
for ALL q not a root of unity (q = 2, 3, ... included).  Paper proves freeness equivalences
only for transcendental q (3 vertices); general question open (their Q 1.17).
A free-for-generic-q, non-DF B_q on 4 vertices that fails at some algebraic non-root-of-unity
q would be a P^3 counterexample.  Search: work/gain_search.py.

## Finding 9: the unique non-DF free line arrangement with moduli in the database
Of 5968 free line arrangements, 7 are non-DF (n = 9 dual Hesse, 15, 15, 17, 18, 27, 27); only
one has moduli: n=15, exps (1,7,7), t=(2^27, 3^6, 5^6), hash 088d4d68...  Structure: the 15
lines joining 6 points Q1..Q6 (the 5-fold points) with 6 of the 15 perfect matchings
{QiQj, QkQl, QmQn} concurrent: {12,34,56},{13,26,45},{16,25,34},{14,26,35},{15,23,46},{15,24,36}.
Realization space (M2, work/real15.m2): irreducible rational curve,
    Q1..Q4 = frame,  Q5 = (d^2 : d : 2d-1),  Q6 = (d(2d-1) : d^2 : 2d-1).
Outside rank-3 verified range (<= 14 lines).  Scan over F_1009: all 1002 same-lattice members
free; 4 values of d acquire 4 extra triple points (lattice change).  [larger primes: pending]
- p=10007 scan of the 15-line family: all 10004 same-lattice members free.

## Finding 10: bias (gain-graphic) arrangements on 4 vertices are TC-safe
work/gain_search.py: ~480k random gain graphs (gains in [-2,2], random coordinate subsets),
~14k with split chi, ~13.4k locally free; ALL locally free ones are free and ALL are DF.
Structured families (all pairs, gain sets {0,1}, {-1,0,1}, {0,1,2}, ...): free only with all
4 coordinate hyperplanes, and then DF.  (Without coordinates: not even locally free.)
Dual "group-like" configurations (torsion points on elliptic/nodal quartics in P^3) fail:
coplanar quadruples dualize to 4 planes in general position through a point (non-free
localization).
- p=32003 scan of the 15-line family: all 32000 same-lattice members free.  => TC holds for
  this (previously unverified, 15-line) lattice beyond reasonable doubt.

## Finding 11: direct P^3 construction search (work/p3_construct.py)
Construction sequences of 11-18 planes (through existing points/lines, 1-2 free parameters),
8 workers x 40 min: 183k arrangements, ~2070 with split chi, ONLY ONE non-DF (and it was not
locally free).  Split-chi lattices built by incidence constructions are essentially always
divisionally free in this range -> TC automatic.  G(3,3,4) (the known non-DF free one) is
isolated: no single-plane addition through a line keeps freeness (restriction sizes 19-k with
k <= 3 never hit {11,12,14}).

## Overall conclusion so far
A P^3 counterexample must be: free, NOT divisionally free, with moduli that change the Ziegler
restriction at EVERY plane (Abe–Yoshinaga Cor 5.3), locally free, no locally heavy flag.
Every construction tried lands in DF (TC-safe) or fails local freeness:
  X3/DW lifts (monodromy), cone+translates (forced by ATW addition), heavy pivots (counting),
  Ziegler-restriction lifts of 12 new multiarrangement families (all 2361 lifts DF),
  bias arrangements (all DF), random {-1,0,1}^4 subsets, direct incidence constructions.
The unique non-DF free line arrangement with moduli in the database (15 lines) satisfies TC.

## Suggested next steps
1. Abstract search for non-DF rank-4 lattices: enumerate rank-4 matroids/lattices with split
   chi, all rank-3 localizations free-compatible, and NO divisional flag; only then realize.
   (DF test is purely combinatorial: chi divisibility along a flag.)
2. Lifts with DF-avoiding concurrency patterns (classes with m_L = 1 have forced restriction
   sizes; choose (B,m) so no plane has |A^P| = 1 + d_i + d_j).
3. Rank 5 (P^4): Ziegler restrictions are rank-4 multiarrangements; DF may be rarer.
4. Products/lifts built on the 15-line non-DF family and on G(r,r,l) "valley" structures.

## Finding 12: which planes witness divisional freeness (work/witness_analysis.out)
(Hit indices here refer to the hit list at analysis time; the family-search workers had
appended further hits, so they differ from the numbering used in Findings 6-7.)
Two hits with valid free lifts, 8 lifts each:
- m = [3,3,1,1,2,1,1], exps (4,4,4): witnesses in classes 0 (m=3; which plane of the class
  varies with the pattern), 2 (m=1), 3 (m=1) -- in all 8 lifts.
- m = [3,4,1,1,1,4,1], exps (5,5,5): witnesses in classes 1 (m=4; plane varies), 4 (m=1),
  6 (m=1) -- in all 8 lifts.
Three other hits produced no valid free lift within 120 greedy attempts.
=> The multiplicity-1 witnesses are forced by the multiplicity data alone (restriction size
   3k - tau_L with tau_L fixed), and each heavy class contributed a witness in every lift.
   A DF-avoiding search should first filter (B,m) so that no multiplicity-1 class has
   restriction size 1 + d_i + d_j, before searching concurrency patterns.

## Finding 13: EXHAUSTIVE rank-3 Terao check over GF(5)  (work/pg2_exhaustive.py)
Exact small-field engine hunter/ffree.py (substitution conditions + GF(p^r) Saito test;
validated vs Macaulay2 over ZZ/5, ZZ/7, ZZ/11 incl. multiarrangements).
All line subsets of PG(2,5) up to PGL(3,5): 7148 orbits (max 948 at k=15,16), 7115 lattice
classes (essential), 4721 free orbits, 27 lattice classes with >1 inequivalent realization
-> 0 classes with mixed freeness.  TC HOLDS FOR EVERY LINE ARRANGEMENT OVER F_5 (all sizes
up to 31 lines; beyond the 14-line literature bound).  GF(3) sanity run: 25 classes, all
uniquely realised (ternary matroids), 0 mixed.

## Finding 14: more finite-field results
- 15-line non-DF lattice: ALL F_p-realizations enumerated (hunter/realize.py) for p <= 113:
  4..110 realizations per prime, all free.  Other 6 non-DF database lattices: 1-2 realizations
  per prime (rigid), all free.
- G(r,r,4) (non-DF free, rank 4): (3,7),(3,13),(3,19),(4,5),(4,13),(4,17),(5,11),(6,7):
  2-4 realizations (Galois-type), all free.
- Over GF(p) the Frobenius derivation sum x_i^p d_i is logarithmic for EVERY F_p-arrangement
  => d1 <= p; free with d1 = p and a full line (p+1 points) => DF.
- Random subsets of PG(2,7) (15-45 lines): ~71% free, but only ~1/10^4 non-DF; the 6 non-DF
  free ones found (15 lines (1,7,7), 16 lines (1,7,8)) are uniquely realised over GF(7).
  GF(11), GF(13): free rare, none non-DF.
- Step 2 (DF-avoiding lifts): prefilter (m=1 classes forced to DF size) passed by 13/25 hits;
  for those, 464 valid free lifts found (7 hits), ALL DF.  The lift construction is DF-bound.
- Rank 4 over GF(5), random 10-22 planes: chi splits ~1/10^4 -> needs structured inputs.

## Finding 15: connecting-line arrangements A(P) (all lines through >= 2 points of P), exhaustive
over k-point sets of PG(2,q) up to PGL (work/connecting_lines.py)
Rich source of the rare class "free, non-DF, with moduli":
- 13 lines, exps (1,6,6), points 3^3 4^3 5^3, moduli 1 (from 6 points)  [<= 14: TC known]
- 15 lines, exps (1,7,7), 3^6 5^6, moduli 1 (the database lattice; 6 points)
- 21 lines, exps (1,9,11), 3^21 6^7 (7 points), rigid
- 23 lines, exps (1,11,11), 3^26 5^4 7^4 (8 points), rigid;  24 lines (1,11,12), rigid
q=11 (k<=8), 13 (k<=7), 17 (k<=7): every realization of every such lattice free (0 mixed).
P^3 analogue (planes spanned by triples of P; localization at Q in P = connecting lines of the
projection from Q): work/spanned_planes.py, GF(5) k<=7: 32 split-chi locally free, all DF.
General-rank lattice (hunter/lattice_gen.py) validated (ranks 3,4 vs old code; braid A5).
- Connecting lines over GF(11), k=9: 129631 orbits, 15887 free, 95 non-DF free, all rigid
  (tangent moduli 0), 0 mixed.  Spanned planes in PG(3,5) up to k=9 (5093 orbits at k=9):
  all split-chi locally free ones are DF.  PG(3,7) k=7: one free non-DF = product C15 x A1
  (6 coplanar points + apex), TC reduces to rank 3.
- GF(7) random hunt: 26 free non-DF (15 lines (1,7,7), 16 lines (1,7,8)), all uniquely
  realised over GF(7).
- 15-line family: at the parameters where the 6 quintuple points (2 values) or the 6 triple
  points (4 values) lie on a conic, the lattice is preserved AND the arrangement stays free
  (no Ziegler-type conic phenomenon).  Expected codimension of the non-free locus (rank drop
  of the degree-6 derivation matrix 98x56) is 43, explaining why generic families never hit it.

## Finding 16: a rigid free NON-DF essential arrangement of 26 planes in P^3 (char 0, Q(omega))
Found by the spanned-plane enumeration over PG(3,7) (k=8).  8 points: e1, e2, e1+a e2 on l1;
e3, e4, e3+b e4 on l2 (skew); X1+X2, X1+c X2 on the transversal t through X1=e1+e2, X2=e3+e4.
A = 26 planes spanned by triples; exps (1,7,9,9); 24 triple lines + 23 quadruple lines;
NO divisional flag.  Occurs iff a, b, c are primitive 6th roots of unity (x^2-x+1=0): found
over GF(7), GF(13), GF(19), none over GF(11), GF(17).  Generic (a,b,c): different lattice,
chi does not split.  Rigid; its 2 GF(7)-realizations are both free -> TC holds for it.
(work/analyze26.py, work/family26.py, work/scan26.py)

## Infrastructure for non-prime fields (hunter/fq.py)
GF(q) table arithmetic (q = p^k <= 256), rank-3 LineArr (lattice, chi, exps, rank-3 DF test,
EXACT freeness via the Ziegler-restriction criterion), realization enumerator over F_q.
Validated: 400 random prime-field cases vs hunter/ffree (0 mismatches); 30 structured split-chi
arrangements over GF(4), GF(8), GF(9) vs Macaulay2 over GF(ZZ/p[s]/f) (0 mismatches).
Bug fixed: hunter/restrict.py hard-coded 3 coordinates (rank 4 only); now rank-generic.

## Finding 17: connecting lines over non-prime fields (work/connecting_fq.py), exhaustive
- GF(8), k<=10: k=10 -> 12870 orbits, ALL split-chi and ALL free (char-2 effect), none non-DF;
  non-DF free lattices: 13 lines (1,6,6) [6 realizations], 17 lines (1,8,8) [3],
  18 lines (1,8,9) [3] -- all realizations free.
- GF(9), k<=9: 19 distinct non-DF free lattices with 13..21 lines (e.g. seven 19-line (1,9,9)
  and five 20-line (1,9,10) types), each with 2-6 realizations -- all free.  0 mixed.

## Finding 18: P^3 and P^4 randomized point configurations (spanned hyperplanes)
- P^3 over GF(11), GF(13) (work/p3_points_random.py): ~50k arrangements, ~3940 split-chi and
  locally free -> ALL divisionally free.
- P^4 over GF(5), GF(7) (work/p4_random.py, rank-5 lattice + rank-5 DF test + exact freeness):
  1723 arrangements, 466 split-chi -> ALL divisionally free.
- GF(5) exhaustive result cross-checked: 40 random orbit representatives (15-28 lines) vs
  Macaulay2 over ZZ/5: 0 mismatches.

## Finding 19: "free => divisionally free" over F_q  (FREE_NOT_DF.md)
FALSE for every q >= 4, in every rank >= 3: A_q = prod_{i<j}(x_i^{q-1} - x_j^{q-1}) in F_q^3
(G(q-1,q-1,3)) is free with Saito basis E, Frobenius sum x_i^q d_i, and
sum (x1x2x3/x_i)^{q-2} d_i (exps (1,q,2q-4)); every line has exactly q points, so no divisional
flag unless q = 3.  Products A_q x Boolean give every rank.  Verified: sympy (determinant),
engine q = 4..16, Macaulay2 q = 4,5,7,8,9,11.
TRUE for q = 2 in rank 3 (hand proof: Frobenius d1 <= q + Ziegler restriction + counting), and
by exhaustion for q = 2 in ranks 3,4,5 (PG(4,2): 1326 essential lattices, 13 non-DF split-chi,
all non-free, Macaulay2-confirmed) and q = 3 in rank 3; q = 3 rank 4 running (111,831 classes).

## Finding 20: free but not DF over EVERY finite field -- blocking sets (FREE_NOT_DF.md)
- Covering lemma: if the normals A* meet every line of PG(V*) (every codim-2 subspace lies in
  some H in A), then A is free with exps (1, q, ..., q^{l-2}, |A| - [l-1]_q).  Proof: delete the
  missing hyperplanes one at a time from the full arrangement B (Moore basis); every restriction
  is the full arrangement of H, so Terao's deletion theorem applies.  Blocking is visible in chi:
  |A|[l-1]_q - b2(A) = [l choose 2]_q.
- Theorem: if moreover A* is a MINIMAL blocking set and d = |A| - [l-1]_q >= q^{l-2}, then A is
  not DF (a dividing chi(A^H) would force |A^H| = [l-1]_q, i.e. a full restriction).
- Examples: A_q (triangle minus vertices; q >= 4, rank 3); dual Hermitian unital over F_{r^2}
  (exps (1, r^2, r^3 - r^2); checked F4, F9, F16); F_3 rank 4: elliptic quadric Q, p in Q,
  tangent line L, A* = L u (PG(3,3) \ (Q u T_p)), 22 planes, exps (1,3,9,9) (not minimal for
  q >= 5); F_2 rank 6: complement of a 13-point complete cap of PG(5,2), 50 hyperplanes,
  exps (1,2,4,8,16,19).
- Exhaustive PG(3,3) finished: 111,802 essential classes, 6,267 split chi, 3,390 DF, 2,877
  non-DF split chi, 7 locally free, all 7 free (22 planes x2, 23 x2, 25 x3); every one is of
  blocking-set type.  So the smallest ranks of free non-DF are 3 (q >= 4), 4 (q = 3), 6 (q = 2).
- Scripts: work/find_hit_3_4.py, find_all_3_4.py, f3_*.py, covering_*.py, caps_pg52.py,
  f2_rank6_cap.py, f2_cap_basis.py, hermitian_check.py, quadric_family_check.py.
- Verification of the covering lemma: Macaulay2, ALL essential covering classes free with the
  predicted exponents -- 3,097 in PG(3,3), 38 in PG(4,2) (work/covering_m2_all.py).  The
  unfiltered PG(3,3) run (exact freeness for all 2,877 non-DF split-chi classes) found the same
  7 classes; the covering ones among the 2,877 are exactly those 7.  Algebraic certificate:
  A = B \ C free iff (if) theta_{F^(l-1)} in span(E, ..., theta_{F^(l-2)}) mod Q_C (validated on
  PG(3,3)); rank-6 F_2 runs: work/f2_rank6_membership.m2, work/f2_rank6_saito.py.
- Side check: the 7 PG(3,3) lattices are rigid in char 3 (one realisation over GF(3), GF(9),
  GF(27)) and have no realisation over GF(4,5,7,8,13,16,25): no Terao pairs from them.
- Rank-6 F_2 example: membership certificate TRUE (M2, 1951s): theta_{F^5} in span(E..theta_{F^4})
  mod Q_K, so free with exps (1,2,4,8,16,19) independently of the deletion theorem.

## Finding 21: proof attempt for Terao's conjecture (TC_ATTEMPT.md)
- No general proof.  Literature frontier (read in full): Abe 2012 (|e1-e2| <= s-2, TC if some
  s_H <= d2-d1+3), Abe 2014 (free => s_H <= d1+1 or = d2+1; TC if d2 >= 2d1-3 or d1 <= 5),
  Barakat-Kuhne 2023 (n <= 14, any char), Dimca-Kuhne-Pokora 2025 (TC if d1 in {m-1,m,n-m};
  counterexample needs d1 >= m+1), Dimca-Pokora 2026 (mdr not combinatorial for 9 lines).
- Anatomy of a P^2 counterexample: n >= 15, not DF, every line d2-d1+4 <= s_H <= d1, d1 >= 6,
  d2 <= 2d1-4, m <= d1-1; refinement: n+1-2r <= s_H <= r+1 with r = mdr(A2) in [m, d1-1];
  rigidity principle (moduli must move the points on every line).
- NEW exact results: TC holds for the 15-line (1,7,7) lattice with moduli (non-free locus
  d^7(d-1)^7(2d-1)^3(d^2-3d+1)(d^2+d-1) inside the degeneration locus; work/tc15_exact.m2);
  TC holds for the 26-plane P^3 lattice (realisation ideal (a-c, b-c, c^2-c+1): one Galois
  orbit; work/exact26_space.py); 19-line (1,9,9) F_11 lattice not char-0 realisable.
- n = 15: exps (6,8) or (7,7); (6,8) forces s_H = 6 on every line and 4 profiles; profile
  (0,25,0,3) is projectively unique = Ceva(5) (3-net + Kneser), so TC.  Others open.
- Finite fields: TC for all line arrangements over F_q, q <= 7 (Barakat-Kuhne + Yoshinaga/Abe
  Cor 6.2); covering lemma rank-3 case = Yoshinaga 2006 (attribution fixed in FREE_NOT_DF.md).

## Finding 22: Terao's conjecture beyond 14 lines -- exhaustive verification by exponent class
(work/tc15/, write-up TC_ATTEMPT.md section 10)
- Method: a counterexample lattice with n lines and exponents (d1,d2) has chi = (t-1)(t-d1)(t-d2),
  all s_H in [d2-d1+4, d1], multiplicities <= d1-1.  Enumerate such linear spaces (big blocks first,
  then triples with exact-cover branching, isomorph rejection by nauty), then decide exactly:
  realisation ideal over Q (frame-fixed dual points), components, Yoshinaga's criterion on the
  Ziegler restriction to z = 0 (rank of the (ax+by)^mu | af+bg conditions in degree d1-1).
- n = 15 (6,8): 12,721 big-block configurations, 3 candidate lattices (two 3-nets of order 5, one
  9-quadruple lattice); non-group net and 9-quad lattice NOT realisable over C; Z5-net = Ceva, one
  degree-4 component, all free.  TC HOLDS for 15 lines with exponents (6,8).  [Theorem 10.1]
- n = 15 (7,7): 99,993 big-block configurations, 6,091,145 candidate lattices (2 GB); the known
  (1,7,7) 15-line lattice is among them (index 1,140,657).  Exact realisability: fast filter
  (python-flint) 4,831,300 NONREAL / 1,259,845 deferred; pass 2 (number fields, curves) 1,245,893
  NONREAL / 13,952 deferred; Macaulay2: 13,942 NONREAL, 10 realisable, 0 mixed.  The 10: 5 rigid
  non-free, 4 rigid free (three degree-2 orbits and one rational 3^6 4^6 5^3 arrangement), and the
  known 3^6 5^6 one-parameter family (generically free, non-free locus inside the degeneracy locus).
  Finite-field cross-check (10 primes) agrees with every exact verdict.
  TC HOLDS for 15 lines with exponents (7,7) [Theorem 10.2]  =>  TC HOLDS FOR ALL 15-LINE
  ARRANGEMENTS IN CP^2 [Corollary 10.3], one beyond Barakat-Kuhne.
- Scale of larger classes (stage 1 only): n=16 (7,8): 8,884,481 big-block configurations (complete,
  4.7 h); n=17 (7,9): > 7.6M states; n=18 (7,10): > 5.7M states (killed: beyond reach of this
  pipeline without further pruning).
- Scaling analysis (TC_ATTEMPT.md section 11): profile inequalities (Hirzebruch, Bojanowski-Pokora)
  pass 99.6% of skeletons; skeleton realisability prunes only the densest (<= 8 big blocks ~100%
  realisable, dims 4-16); in-search realisability pruning: -28% nodes at 20x cost; 83% of the 6M
  (7,7) lattices have m = 5 and s_min = 5 (where the free ones are).  New lattice condition from
  Abe 2014 Cor 1.2(2): unjoined multiple points P, Q need (m_P-1)+(m_Q-1) <= d1 or = d2 (prunes
  1.3%; certifies lattices 3143, 4608 non-free).  n = 16 feasible (~400 CPU-h) with streaming +
  compiled filter; beyond needs theory (eps >= 2 free arrangements; mdr-mixed deletion lattices).
- n = 16 (7,8) [the only class]: 8,884,481 skeletons, 16,371,552 candidate lattices (6 GB); fast
  filter 13,748,942 NONREAL / 2,622,610 deferred; pass 2 2,594,921 NONREAL / 27,689 to M2; M2: 3,039
  ideal (1), 24,650 degenerate-only, 0 REALISABLE.  TC HOLDS FOR ALL 16-LINE ARRANGEMENTS IN CP^2
  [Theorem 10.4].  ~4 h on 8 cores (volume ratio lattices/skeletons 1.8, vs 61 for (7,7)).

## Finding 23: general-n proof attempt (TC_GENERAL.md)
- Window for all n collected (d1 >= 6, delta <= d1-4, delta+4 <= s_H <= d1, m <= d1-1, m <= r <= d1-1,
  unjoined-pair rule, moduli must move every line).
- Narrow classes delta = d1-4 (n = 3d1-3: (6,8), (7,10), (8,12), ...): s_H = d1 on every line,
  r = mdr(A2) = d1-1, every Ziegler restriction of A2 has exponents (d1-1, d2+1) (defect d1-3 on
  every line), Sum m_P = 3d1(d1-1), Sum m_P^2 = 12(d1-1)^2, N = d1^2-2d1+4; Melchior gives
  2N-b2 = 12-3d1 < 3, so NO real/rational counterexample in a narrow class for any n.  The minimal
  derivation theta of A2 gives omega = theta x E, a section of Omega^1(d1+1), whose restriction to
  every line of A2 is (normal) x c_H prod alpha_P: vanishes exactly at the multiple points, simply;
  c2 = d1^2-d1+1 = N + (d1-3).  Could not derive a contradiction from this.
- Reduction: for a component C with a free point, the non-free locus equals pi_H^{-1}(Sigma_H) for
  EVERY line H (Yoshinaga), so all n jump loci coincide; TC for L <=> Conjecture J for ONE line
  (jump on a line forces an extra incidence).  J proved only for the 15-line family (jump polynomial
  divides the degeneration polynomial with smaller exponents) and vacuously for rigid lattices.
- Honest status: no proof; every known tool is numerical and the numbers are consistent for all
  d1 >= 6; the verified cases are decided by non-realisability of dense linear spaces, not by the
  window.

## Finding 24: the foliation of the minimal derivation; 3 mdr(A2) >= n+2 (TC_FOLIATION.md)

- Tool: the least-degree Jacobian syzygy theta of the non-free member defines a degree-r foliation
  of P^2 (r = mdr(A2)) with all n lines invariant; saturation of omega = theta x E follows from
  s_H <= d2+1 and r <= d1-1. Its inflection curve J = det(d omega / dx) (Pereira's first extactic
  curve, degree 3r) is divisible by Q and vanishes to order >= 2 at every singular point of the
  foliation; so n <= 3r and J = Q J' with deg J' = 3r-n.
- Excess identity (Thm 4.1): (d1-r)(d2-r) = sum over off-line singular points of mu_P minus
  sum over on-line singular points of d_P, with d_P = 0 at nondegenerate points and at
  saddle-nodes along a line. Every point contributing to the right side is a singular point of J'.
- Theorem 5.1: a counterexample has 3 mdr(A2) >= n+2, hence d2 <= 2 d1 - 6, d1 >= 7, n >= 15
  (theory only).  Excluded for every n: delta = d1-4 (narrow classes) and delta = d1-5.
  Consequences: (6,7) at 14 lines proved theoretically (TC for n <= 14 now theoretical);
  (6,8) at 15 lines theoretical (was Thm 10.1 by enumeration + Kneser); (7,9) and (8,8) with r = 6
  at 17 lines excluded -- exactly what the running n=17 enumeration covers, which thereby becomes a
  consistency check; the only open 17-line class is (8,8)[B] (A2 nearly free (7,10), r = 7,
  deg J' = 4, excess 1).  (7,10) at 18 lines excluded; (8,9) survives with r = 7.
- Checked in Macaulay2 on free arrangements with n = 3r and 3r-1 (work/foliation/inflect.m2): J/Q
  constant resp. a line, zero scheme of degree r^2+r+1, no off-line points, Delta = 0.
- Not proved: TC. Open: deg J' >= 2, i.e. whether the singular points of a curve of degree
  3r-n <= d1-delta-4 can carry the excess (d1-r)(d2-r) >= delta+1; first instances j = 2
  ((7,8), (8,10)) and j = 3 ((7,7), (8,9)).
- Bourbaki link (same day): Nasrollah Nejad-Nouri (arXiv:2609.22453) have Bour(A) = chibar(mdr),
  so for A2 the excess (d1-r)(d2-r) IS the Bourbaki degree; their Bourbaki ideal contains the
  residual polar curves and the cofactor ideal (K_H), so the Bourbaki scheme is supported at the
  excess points (common zeros of the cofactors) inside Sing(J/Q). Checked on the Dimca-Pokora
  9-line pair (work/foliation/bourbaki.m2): Bourbaki scheme = excess scheme (2 nodes of J/Q;
  resp. 2 off-line points + 1 degenerate double point). Their Terao criterion needs
  n >= 3 d1 - 1, weaker than Abe inside the window; Theorem 5.1 is strictly stronger there.
- Theorem 10.1 (same day, later): deg(J/Q) = 2 is impossible, so 3 mdr(A2) >= n+3 and
  d2 <= 2 d1 - 7. Proof: (i) E|_H = O(-e1(H)) + O(-e2(H)) (Ziegler module = restricted syzygy
  bundle), so the Bourbaki points on a line are the zeros of a form of degree r - e1(H);
  (ii) Abe's gap bound gives r - e1(H) <= j/3, so for j <= 2 no Bourbaki point is on a line;
  (iii) J/Q a conic: smooth -> no excess; line pair -> Bour = 1 -> n = 13; double line ->
  every Bourbaki point is a Jordan-block singular point with eigendirection the line, the line is
  then invariant (local computation), and along an invariant double line of J the ratio
  (transverse eigenvalue)/(tangential field) solves phi' = -phi^2, so at most one Jordan point.
  Consequences: (7,8) at 16 lines theoretical (TC for 16 lines now theoretical), (8,10) and
  (9,9)[r=7] at 19, (9,12) at 22 excluded. Open: (7,7) [15, computer], (8,8)[B] [17, j=4,
  Bour 1], (8,9) [18, j=3, Bour 2].


## Finding 25: the (8,9) seed at 18 lines and the three-pencil neighbourhood (TC_FOLIATION.md §11)

- Seed (pasted by the user): Q = yz(x^6-y^6)(z^6-x^6) prod_{b in mu6 \ {1,zeta}} (y - b z), 18 lines
  = full monomial A(6,1,3) (21 lines, free (7,13)) minus x, y-z, y-zeta z. Verified
  (work/pencil/seed18.m2): profile 2^24 3^24 6 7^2, b2 = 89, chibar = (t-8)(t-9), syzygy degrees
  (7,11,12), mdr 7, j = 3, Bourbaki scheme = the two points [0:1:1], [0:zeta:1] off all lines,
  J/Q = x(y-z)(y-zeta z) = product of the deleted lines (the degree-7 derivation of A(6,1,3) has
  3*7 = 21 invariant lines, the maximum). Same for the two other variants. These are the only
  3-line deletions of A(6,1,3) in the class (8,9) (22 = 8+7+7 points must be removed).
- Lattice rigid (work/pencil/real18.m2: binomial realization scheme over QQ has 6 points, only
  the Galois pair zeta6^{+-1} non-degenerate), so no lattice-mate; Kneser's theorem forces mu6.
- Three-pencil lattices (lines x = a y, z = c x, y = b z, triple points iff abc = 1, coordinate
  lines K): realization space = binomial solution group mod torus, components = torsion of the
  relation lattice (Smith normal form), generic points explicit (work/pencil/pencil_enum.py,
  pencil_iso.py, pencil_jump.py). For n = 18, class (8,9), m <= 7: six patterns; (6,6,4;yz;24) is
  the seed; (6,5,5), (6,6,5), (6,5,6), (6,6,6) have no realizations (Kneser); (5,5,5;xyz;19) gives
  six lattices up to isomorphism (3 rigid Z/7, 2 rigid Z/8, one 1-parameter family
  x = g^i y, z = g^j x, y = g^-k z, xyz), none with two components, ALL FREE (8,9); the family is
  free at every non-degenerate g: the 16-minors of the 17x16 Yoshinaga matrix of the
  multirestriction on x over QQ[g], saturated by the degeneration binomials, give (g^45).
- Conclusion: TC holds for all three-pencil lattices with 18 lines in class (8,9); the only
  non-free ones are the seed variants, which have no lattice-mates. Weakening the cyclic
  relations inside the monomial world cannot produce a counterexample; one must leave it.

## Finding 26: Terao over fields of positive characteristic -- the Frobenius bound (TC_CHARP.md)

- An arrangement defined over F_q (q = p^k) with p not dividing n has the Jacobian syzygy
  theta' = (x^q, y^q, z^q) - (sigma/n) E, sigma = sum alpha_H^(q-1), of degree q: mdr <= q. Hence an
  F_q-rational arrangement in a lattice class with d1 > q is non-free, and a counterexample over
  F_p-bar needs only a free lattice-mate (necessarily not F_q-rational). Verified in M2.
- Searches (work/charp/): three-pencil lattices over F_q* for q = 7, 8, 9, 11 with d1 > q (36, 60,
  114, 8324 lattices): all rigid or (30 at q = 11) one-parameter families whose generic member
  keeps mdr = 11; additive three-pencil over F_5, F_7: rigid; ALL sub-arrangements of PG(2,4)
  (4 lattices, 11 lines, (5,5), rigid, syz (4,7,7)) and of PG(2,5) (203 lattices, 13-14 lines,
  (6,6)/(6,7), all rigid F_5-points, syz (5,8,8)/(5,9,10)); one-parameter families over Z: the
  smallest admissible q is always >= d1 (equality in 299 cases at D = 6), never below.
- No counterexample and no pair with different syzygy degrees found; every lattice the mechanism
  applies to is rigid.

## Finding 27: scheme-theoretic Bourbaki-in-Sing(J/Q) (TC_FOLIATION.md section 12)

- From a user-supplied proof attempt: the Bourbaki scheme of the minimal derivation is contained
  SCHEME-THEORETICALLY in the singular scheme of the residual inflection curve C = J/Q (local
  equation of C in the square of the Bourbaki ideal; proof via the Saito frame R = Euler,
  H = Hamiltonian of the local product of lines, X = fR + gH, D/q = -m(fX(g)-gX(f)+(m-1)f^2 g)
  + kappa g^3 in (f,g)^2). Verified algebraically and computationally (deletions of A(3,1,3),
  A(4,1,3); Dimca-Pokora pair; random arrangements over F_101 with non-reduced Bourbaki
  components of length 2, 3, 5): work/foliation/schemelemma.m2, nonreduced.m2.
- Consequence: C reduced => ac <= tau(C) <= (j-1)^2; for j = 3: ac <= 3, leaving (7,7)@15,
  (8,9)@18, (9,11)@21 in the reduced case. Non-reduced C: partial (C = l^2 l', off-line Bourbaki
  points on l \ l': at most one, by Lemmas 10.4(a),(b)); l ∩ l', on-line points and C = l^3 open.
- Does not touch the free lattice-mate; TC not proved. The missing step is exactly the
  Fitting-ideal saturation statement on the realization space (note's (14)-(15)).
