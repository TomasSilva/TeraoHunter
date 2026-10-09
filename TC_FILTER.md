# Combinatorial filters for the residual 17-line class (8, 8)[B]

**Status (9 Oct 2026).** Work in progress. The aim is to make the enumeration of the last open
17-line class feasible: lattices with χ = (t − 1)(t − 8)², multiplicities ≤ 6, every line with
4 to 8 points and at least one line with 4, 5 or 8 points, in which the non-free member A₂ would be
nearly free with exponents (7, 10) and the free member A₁ has exponents (8, 8)
(`TC_ATTEMPT.md` §11, Theorem 10.5). Everything here is a necessary condition on such a lattice
or on its realizations. Scripts in `work/filter/`.

## 1. The per-line invariant that separates the two members

Ziegler's restriction to a line H ∈ A is the multiarrangement of the points of H with
multiplicities m_P − 1; its exponents (e₁(H), e₂(H)), e₁ + e₂ = n − 1 = 16, depend on the
positions of the points, and by Yoshinaga's criterion A is free iff e₁(H) = 8 for one (equivalently
every) line H. In a counterexample A₁ has (8, 8) on every line, while A₂, being nearly free with
exponents (7, 10), has (7, 9) on every line not through its jumping point P₀ and (6, 10) on the
lines through P₀ (Abe–Dimca; `TC_FOLIATION.md` Lemma 10.2–10.3; checked on a 13-line nearly free
arrangement in `nf13.m2`, whose jumping point lies on no line and whose 13 lines all have (5, 7)).

The multiplicity vector μ = (m_P − 1)_{P ∈ H} of a line is lattice data (a partition of 16 into
s_H parts ≤ 5). For each such vector the locus of positions where e₁ ≤ 7 (resp. ≤ 6) is the rank
drop of an explicit 16 × 16 (resp. 16 × 14) matrix in the s_H − 3 moduli of the points
(`jumps2.m2`, `jumps4.m2`, `jumps8.m2`; the matrices are built from the Taylor coefficients
(1/k!)∂_u^k h(p, 1), which is correct in every degree).

**Lemma F1 (jump loci).**
- s_H = 4 (vectors (4,4,4,4), (5,4,4,3), (5,5,3,3), (5,5,4,2), (5,5,5,1)): e₁ ≤ 7 happens exactly
  at finitely many cross-ratios t of the four points, and e₁ ≤ 6 never. The jump sets, with the
  points at ∞, 0, 1, t and the multiplicities in that order:
  (5,5,5,1): t ∈ {−1, 2, 1/2} (harmonic) or t² − t + 1 = 0 (equianharmonic);
  (5,5,3,3): t = −1, or t⁴ − 8t³ + 7t² − 8t + 1 = 0, or t⁴ + 10t³ − 23t² + 10t + 1 = 0;
  (5,4,4,3): t = 1/2 or a root of 2t⁸ − 8t⁷ + 93t⁶ − 251t⁵ + 266t⁴ − 123t³ + 76t² − 55t + 1;
  (5,5,4,2): a root of 2t⁸ − 4t⁷ − 3t⁶ + 48t⁵ − 85t⁴ + 48t³ − 3t² − 4t + 2;
  (4,4,4,4): a root of t² − 5t + 1, t² − t + 1, t² + 3t − 3 or 3t² − 3t − 1.
- s_H = 5 (eleven vectors): e₁ ≤ 7 on an irreducible curve of degree 8 to 10 in the two
  moduli, and e₁ ≤ 6 never.
- s_H = 6: e₁ ≤ 7 on a surface and e₁ ≤ 6 at finitely many positions (for the vectors computed
  so far: (3,3,3,3,2,2), (3,3,3,3,3,1), (4,3,3,2,2,2), (4,3,3,3,2,1), (4,4,2,2,2,2); the run is
  being completed). So a line through the jumping point P₀ with 6 points sits at one of finitely
  many positions of its six points.
- s_H = 7, 8: being computed.

Consequences. (i) The jumping point P₀ of A₂ lies on no line with 4 or 5 points; the lines
through it have 6 or 7 points (this sharpens `TC_FOLIATION.md` Lemma 10.3 for j = 4). (ii) A line
with 4 points in A₂ must have its four points in one of the finitely many special cross-ratios
above; a line with 5 points in A₂ has its five points on an explicit curve.

**Lemma F2 (freeness read off a 4-point line).** Let L be a lattice of the class above with a line
H of 4 points and multiplicity vector μ. A realization of L is free if and only if the cross-ratio
of the four points of H is not in the jump set of μ. (Proof: Yoshinaga's criterion on H and
Lemma F1 for s_H = 4.) So for such lattices freeness is a single explicit algebraic condition on
one cross-ratio, and Terao's conjecture for L says: on the realization space of L the cross-ratio
of H either never meets the jump set or always does.

## 2. The line with the fewest points coordinatizes the lattice

Every other line meets H in one of its s_H points, so the arrangement is H together with s_H
pencils with collinear centres. Put the centres at [1:0:0], [0:1:0], [1:1:0], [λ_i:1:0] on
H = {z = 0}; the pencils are y = a z, x = c z, x − y = b z, x − λ_i y = d_i z, and a concurrency
of three lines from three pencils is a *linear* equation in the pencil parameters (c − a = b,
c − λa = d, a + b − λa = d, c − λ(c − b) = d, and the analogous ones with two λ's), with the
cross-ratios λ_i as the only non-linear parameters. Hence:

**Lemma F3.** The realization space of a lattice whose minimal line has s points is the solution
set of a system that is linear in the pencil parameters for fixed values of the s − 3 cross-ratios
of the points of that line, modulo a 3-dimensional affine group. For s = 4 it is a family of
linear spaces over the λ-line, and for s = 5 over the (λ₁, λ₂)-plane.

For the vector (5,5,5,1) this is a 3-net (three pencils of five lines, additive structure
c − a = b, t triple points) with the common line H and one extra line L₄ through the double point
P₄ = [λ:1:0]. Counting b₂ = 80 gives t = 10 + s_{L₄}, so t ∈ {16, 17, 18} and L₄ passes through
at least five triple or double points of the net; each pencil line carries at least three triple
points. In characteristic zero the difference set B − A of two 5-element sets has at least nine
elements (Cauchy–Davenport in a torsion-free group), so sixteen or more of the 25 pairs (a, c)
having c − a in a 5-element set C forces A, B, C to be near arithmetic progressions.

With Lemma F2 this sub-case is decidable lattice by lattice: compute the λ-fibres of the
realization space; a counterexample would be a lattice whose realization space contains both a
point with λ in the jump set and one outside it, with the same lattice at both.

**Proposition F4 (the (6,6,6,2) sub-case is empty).** No 17-line lattice of class (8, 8) with
multiplicities ≤ 6 and a line of 4 points of multiplicities 6, 6, 6, 2 is realizable over ℂ.
*Argument.* By Lemma F3 the realization is a 3-net with additive structure: A, B, C ⊂ (ℂ, +)
of size 5, the triple points being all pairs (a, c) with c − a ∈ C. The count above gives
t ∈ {16, 17, 18}, and every pencil line needs at least three triple points (a pencil line
through P₁ meets the other two pencils in ≤ 7 points, 5 + (5 − t_a) + [L₄] ≤ 7). So each a ∈ A
has at least three c with c − a ∈ C, each c at least three a, each b at least three
representations. Ordering the three sets in ℝ through a generic linear functional, the extreme
elements force the middle label b₃ to equal c₁ − a₁ = c₅ − a₅, hence A and B to have the same
diameter Δ, and a₁ + {b₃, b₄, b₅}, a₅ + {b₁, b₂, b₃} ⊂ B, c₁ − {b₁, b₂, b₃}, c₅ − {b₃, b₄, b₅} ⊂ A,
which pins A, B, C to near-arithmetic progressions; for those the triple count with all degrees
≥ 3 is at least 19. The exhaustive search `add551.py` over A, B ⊂ [0, D] (after translation),
D ≤ 16, with C any five differences, finds no structure with t ≤ 18 and all degrees ≥ 3. ∎
(computer-assisted; the ordering argument is a sketch)

The other four 4-point vectors are 4-pencil structures with pencils of sizes (5,5,4,2),
(5,5,3,3), (5,4,4,3), (4,4,4,4) and fewer forced triple points (for (5,5,4,2): t = s₁ + s₂ − 1
∈ [11, 15] with s_i the point counts of the two lines through P₄), so their additive nets may
have moduli (dim V > 3); there a counterexample is conceivable exactly when the cross-ratio λ,
a rational function on V, is non-constant and the lattice survives on a fibre λ ∈ S(μ). Note that
for (5,5,4,2) and (4,4,4,4) the jump set S(μ) contains no rational number, so a lattice whose
realization is rational and rigid is free; for (5,5,5,1), (5,5,3,3), (5,4,4,3) the rational
jump values are −1 (harmonic) and 1/2. These four sub-cases are the next computation. A first pass (`net4gen.py`, additive nets with
A, B ⊂ [0, 9] and rigid, i.e. dim V = 3) finds no complete lattice at all: for every rigid net
no choice of the lines through P₄ satisfies the window and b₂ = 80. What remains are the
non-rigid nets (dim V ≥ 4; 236, 6991, 22997 and 897890 of them in the box for the four vectors,
with as few as 7 to 12 triple points), for which the lines through P₄ impose bilinear equations
(linear in the pencil parameters, with λ as coefficient) that must be solved symbolically; since
such nets are not near-progressions, they have to be enumerated as abstract partial Latin
rectangles rather than in an integer box. This is the open part of sub-case (I).

## 3. The reduced enumeration plan

Part [B] splits by the minimal number of points on a line:
- (I) a line with 4 points: 4-pencil structures with collinear centres, one cross-ratio;
  decided by Lemma F2 and the linear realization spaces of Lemma F3.
- (II) a line with 5 points and none with 4: 5-pencil structures, two cross-ratios, with the
  jump curve of Lemma F1 playing the role of the jump set.
- (III) all lines with 6 to 8 points and at least one with 8: the window w(x) = 16 − s_x ∈ [8, 10]
  of the finished run widened by one value; it is empty if no 8-point multiplicity vector admits
  e₁ ≤ 7 (being computed).

(to be continued as the computations finish)
