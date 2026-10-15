# Terao's conjecture over fields of positive characteristic: the Frobenius bound and a search

**Status (8 Oct 2026).** Not a counterexample. This records one structural fact that is special to
characteristic p, the search it suggests, and what the search found.

## 1. The Frobenius syzygy

Let A be an arrangement of n lines in P² defined over F_q (q = p^k), with p ∤ n, Q = ∏ α_H.
Since a ↦ a^q is the identity on F_q and additive, θ_q := x^q ∂_x + y^q ∂_y + z^q ∂_z satisfies
θ_q(α_H) = α_H^q for every line, hence θ_q(Q) = Q·σ with σ = Σ_H α_H^{q−1}. Therefore

  θ′ := θ_q − (σ/n)·E,   θ′(Q) = 0,

is a Jacobian syzygy of degree q, nonzero and not a multiple of E (its x-component is
x(x^{q−1} − σ/n)). So

> **Every arrangement defined over F_q with p ∤ n has mdr ≤ q.** Consequently an F_q-rational
> arrangement whose lattice class (d₁, d₂) has d₁ > q is non-free, for any lattice-mate reason or
> none: free ⟹ mdr = d₁ (Terao's factorization holds over every field, and p ∤ n gives
> D(A) = S·E ⊕ D₀(A)).

Checked in `work/charp/frob_check.m2`: eleven random lines over F₄ have syzygy degrees (4, 5, 5),
thirteen random lines over F₅ have (5, 5, 5), and θ′ is verified to kill Q; the 18-line seed of
`TC_FOLIATION.md` §11 reduced mod 7 (ζ₆ ∈ F₇) has degrees (7, 11, 12).

**The route to a counterexample this opens.** If a lattice L has an F_q-rational realization A₂
with d₁(L) > q, then A₂ is non-free for free. A counterexample over F̄_p then needs only a
*free* realization A₁ of L over F̄_p (necessarily not F_q-rational, i.e. over F_{q^k} ∖ F_q).
Equivalently: a positive-dimensional (or multi-component) realization space over F̄_p of a
lattice with d₁ > q whose generic member is free and which contains an F_q-point with the same
lattice. Note the arithmetic: d₁ > q forces n ≥ 2q + 3 lines inside PG(2, q), a dense
sub-configuration, while a free generic member needs the lattice to be loose enough to move.

**A consistency check that is not a coincidence.** In the free one-parameter family of
`TC_FOLIATION.md` §11.2 (x = gⁱy, z = gʲx, y = g⁻ᵏz, xyz; exponents (8, 9)) the lattice survives
reduction to F_q exactly when F_q* has an element of order ≥ 7, i.e. q ≥ 8 = d₁: the F₈-member is
free with the Frobenius syzygy as its minimal one. The additive analogue (y = az, x = cz,
x − y = bz over the arithmetic progressions, with the common line z; 16 lines, class (7, 8))
survives reduction mod p exactly for p ≥ 7 = d₁. Both families miss a counterexample by one.

## 2. Searches

All F_q-rational lattices with d₁ > q in the following families were enumerated and their
realization spaces over F̄_p computed; a counterexample would be a non-rigid one with a free
generic point.

- Three-pencil lattices over F_q* (`fq_pencil.py`; multiplicative structure, Smith normal form
  of the relation lattice, torsion prime to p): q = 7 (36 lattices), 8 (60), 9 (114): all rigid,
  every realization component is Galois-conjugate to the F_q-point; M2 confirms mdr ≤ q on all
  of them (q = 7: syzygy degrees (7, 10, 10) or (7, 11, 12)). q = 11: 8324 lattices (25 to 28
  lines, classes (12, 12) to (13, 13)); 8294 are rigid with all components F₁₁-rational; the other
  30 (25 lines, class (12, 12), e.g. A = F₁₁* ∖ {1}, B = F₁₁*, C = the non-squares, K = {x}) have a
  one-dimensional realization space over F̄₁₁, and its generic member, over GF(11⁴), still has
  syzygy degrees (11, 14, 14): the Frobenius syzygy deforms along the family and nothing is free.
- Additive three-pencil lattices over F_p (`additive.py`; collinear centres, linear realization
  space): p = 5 (one lattice, 13 lines, class (6, 6)), p = 7 (two lattices, 17 lines, class
  (8, 8)): all rigid.
- All sub-arrangements of PG(2, 4) (`pg2q_enum.py`, exhaustive over the 2²¹ subsets with p ∤ n):
  only n = 11, class (5, 5), four lattices up to the invariants used; their realization schemes
  over F̄₂ (`pg2q_real.py`, `real_pg4_*.m2`) are single Galois orbits of the F₄-point, with syzygy
  degrees (4, 7, 7) (nearly free), so no lattice-mate at all.
- All sub-arrangements of PG(2, 5) (`pg25scan.c`, exhaustive): see §3.
- One-parameter families over Z with an admissible prime power q < d₁ (`spread.py`): see §3.

## 3. Results of the larger scans

**All sub-arrangements of PG(2, 5)** (`pg25scan.c`, every n-subset of the 31 lines with 5 ∤ n,
about 1.5·10⁹ subsets). Integer-rooted χ̄ with d₁ > 5 occurs only for n = 13 (class (6, 6),
96 million subsets) and n = 14 (class (6, 7), 17 million); nothing for n ≥ 16. Up to the invariant
(point profile, sorted multiplicity tuples of the points on each line) these are 203 lattices
(161 + 42; `pg25_hits.txt`). Their realization schemes over F̄₅ (`pg2q_real.py`,
`pg25_real_results.txt.gz`) are all zero-dimensional and F₅-rational: 202 consist of the single
F₅-point, one (13 lines, profile 2¹² 3¹² 5³) of two F₅-points. Syzygy degrees (5, 8, 8) for every
13-line lattice (nearly free) and (5, 9, 10) for every 14-line one. No lattice-mate exists at all.

**One-parameter three-pencil families over Z** (`spread.py`, `spread_diag.py`; A, B ⊂ [0, D],
C the sums or differences carrying most triple points, every K; D ≤ 7). For every family with
integer-rooted χ̄ the smallest prime power q admitting a reduction with the same lattice
satisfies q ≥ d₁, with equality for 299 of the families at D = 6 and strict inequality for the
rest; no family admits q < d₁. (At q = d₁ the F_q-member may well be free, with the Frobenius
syzygy as its minimal one, as in §1.) So in this world the Frobenius bound is sharp but never
violated, exactly as Terao's conjecture demands.

**Summary.** Over F̄₂ (all of PG(2, 4)), F̄₅ (all of PG(2, 5)), and for every three-pencil or
additive lattice over F₇, F₈, F₉, F₁₁ with d₁ > q, the F_q-rational arrangement that the
Frobenius syzygy makes non-free has no free lattice-mate, because its lattice is rigid or its
only deformations keep a degree-q syzygy. The mechanism is sound (it is a clean proof of
non-freeness) but every lattice it applies to in these ranges is rigid. The remaining room is
q = 7 with general (non-pencil) sub-arrangements of PG(2, 7), which cannot be enumerated
exhaustively (C(57, 17) ≈ 10¹⁴). A random sample (`pg27sample.c`, 3·10⁷ subsets of 17 lines)
contains 4.2 million subsets in class (8, 8), i.e. with d₁ = 8 > 7, in 1.66 million invariant
classes; the 25 most structured of them (fewest double points; `pg27_hits.txt`,
`pg27_real_results.txt.gz`) all have realization scheme over F̄₇ equal to the single F₇-point,
with syzygy degrees (7, 10, 10). So the sample behaves like PG(2, 4) and PG(2, 5): the
automatically non-free F₇-arrangements are rigid.

## 4. The lattices of the 17-line pencil censuses (15 Oct 2026)

The censuses of TC_BUILD.md §5.2 and §5.5 work over F_p only as a route to characteristic 0. Their
exact lifting decides freeness over ℚ̄. Here the same hits and lattices are decided in
characteristic p (`work/charp17/`). Freeness is tested by the minimal free resolution of the
Jacobian ideal over the field in question: free iff pd(S/J) = 2.

- **Every hit is non-free over its own field** (`charp_free.py`, `charp_free.out`). All 1,328
  F_p-configurations have syzygy degrees (7, 10, 10), i.e. are nearly free (7, 10):
  - 712 hits at p = 11 and 592 at p = 13 with a 5-point line;
  - the 24 E₁₇ hits at p = 13.

  This includes the 480 hits of the two characteristic-11 lattices (classes 3 and 5), which have no
  realisation in characteristic 0 or 13.
- **Realisation spaces over F̄₁₁ and F̄₁₃** (`charp_lift.py`, `charp_lift.out`) of the seven census
  lattices, E₁₇ and the six five-pencil classes, are all zero-dimensional:
  - in characteristic 11 the ℚ(√−3)-lattices (E₁₇, classes 1 and 2) are each a single point over
    F₁₂₁, and the other four lattices are each a single F₁₁-point;
  - in characteristic 13 the ℚ(√−3)-lattices split into two F₁₃-points each, and classes 3 and 5
    have no realisation.

  Every point is nearly free (7, 10). So Terao's conjecture holds for these lattices in
  characteristics 11 and 13, as well as 0.
- **The 236 lattices of the F₇ run** (`charp7_lift.py`, `charp7_lift.out`). These come from the
  four-pencil (4, 4, 4, 4) run at cross-ratio −1; none is realisable over ℚ̄. They are class (8, 8)
  with d₁ = 8 > 7, so by §1 every F₇-rational member is non-free. A characteristic-7 counterexample
  would therefore need a free realisation over F̄₇ ∖ F₇. There is none:
  - each realisation space over F̄₇ is the single F₇-point;
  - for 10 of the lattices there is also a second component, but it lies in a forbidden
    collinearity;
  - every F₇-point is nearly free with syzygies (7, 10, 10). Its degree-7 syzygy is, up to a scalar,
    the Frobenius syzygy θ′.
- **Control** (`g17_control.m2`). The free harmonic grid G₁₇ (TC_BUILD.md §3), reduced mod 11 and
  mod 13, is free with syzygies (8, 8) and pd 2.

So no census lattice is a counterexample in characteristic 0, 7, 11 or 13. No free realisation of
any of them exists over ℚ̄, F̄₇, F̄₁₁ or F̄₁₃.

## 5. A one-line Frobenius criterion (15 Oct 2026)

**Proposition 5.1.** Let A be an arrangement of n lines over a field of characteristic p, with χ̄(A, t)
having integer roots d₁ ≤ d₂. Suppose some line H ∈ A has the following properties:
- the points of A on H are projectively F_q-rational (q = pᵉ), i.e. in a suitable coordinate on H
  they all lie in P¹(F_q);
- every Ziegler multiplicity m(X) = (number of lines through X) − 1 on H is at most q.

If q < d₁, then A is not free.

*Proof.* D(A) = S·θ_E ⊕ D_H(A) in every characteristic, so if A is free with exponents (1, d₁, d₂)
then Ziegler's multirestriction (A^H, m) is free with exponents (d₁, d₂). The proof via Saito's
criterion has no characteristic assumption, and Terao's factorization is also characteristic-free.
In coordinates (x, y) on H where every point is F_q-rational, θ_q = x^q∂_x + y^q∂_y satisfies
θ_q(α_X) = α_X^q, which is divisible by α_X^{m(X)}. So θ_q ∈ D(A^H, m) has degree q, the smaller
exponent of the multirestriction is at most q < d₁, and this is a contradiction. ∎

Remarks on the proposition:
- The global bound of §1 is the special case where every line is F_q-rational.
- The new point is that a bad lattice-mate only needs ONE line in Frobenius position. Its other
  coordinates can be transcendental.
- A line with at most 3 points is always F_p-rational. So in characteristic p, a line with ≤ 3
  points, all multiplicities ≤ q and q < d₁ makes the lattice non-free in every realisation. This
  is a combinatorial statement, consistent with TC.

## 6. Why every family so far missed by one

**Lemma 6.1.** For q ≤ 7, any four distinct points of P¹(F_q) are harmonic or equianharmonic.
- F_q ∖ {0, 1} has q − 2 elements.
- These are covered by the harmonic orbit {−1, 2, ½} and the roots of μ² − μ + 1:
  - q = 3: {2};
  - q = 4: {ω, ω²};
  - q = 5: {2, 3, 4};
  - q = 7: {2, 4, 6} ∪ {3, 5}.

For q = 8, the six elements of F₈ ∖ F₂ form one generic anharmonic orbit; neither special value
lies in F₈. For q = 9, 11, 13, F_q has 6 generic values.

So for q ≤ 7 every F_q-position of a moving cross-ratio is the reduction of a characteristic-0
special position, and there TC-consistent families already degenerate:
- ACKN-13: the degeneracy set {0, 1, ∞, −1, 2, ½, λ² − λ + 1 = 0} covers P¹(F₄), P¹(F₅) and P¹(F₇);
- ACKN-15: {0, ½, 1, ∞, λ² − 3λ + 1, λ² + λ − 1} covers F₄ and F₅;
- the (8, 9) three-pencil family: the elements of F_q* have order ≤ q − 1 ≤ 6.

This is the "miss by one" of §1 and §3. The Frobenius mechanism (§1 or §5) can only produce a
counterexample when q ≥ 8, hence d₁ ≥ 9 and n ≥ 19. The first fields worth searching are F₈, F₉
and F₁₁.

## 7. Arrangements inside PG(2, q)

**Lemma 7.1.** Let A ⊂ PG(2, q) contain a line H, and let A ∖ H be K = n − 1 affine lines covering U
of the q² affine points. Then:
- b₂ = (q + 1)K − U, so χ₀(A, 0) = qK − U;
- χ̄ splits iff W = q² − U = u·v with u + v = K − 2q, and then the exponents are (q + u, q + v).

So d₁ > q iff the affine lines leave uv > 0 points uncovered with u ≥ 1. Covering every affine
point gives exactly d₁ = q, which is why the Frobenius bound is so often attained but not
exceeded. (For the 13- and 14-line PG(2, 5) lattices of §3 this gives W = 1 and W = 2.)

## 8. Exhaustive four-pencil searches over F₈ and F₉ (`work/charp2/`)

An arrangement with a 4-point line H is a four-pencil arrangement:
- every other line passes through one of the 4 points of H;
- for fixed cross-ratio μ the offsets satisfy linear equations M(μ)s = 0, one for each extra line
  through an affine point;
- so the realisation space fibres linearly over the μ-line.

With Lemma 7.1, every four-pencil sub-arrangement of PG(2, q) with d₁ > q can be enumerated.
Every such arrangement is non-free by §1. For each one the question is whether its lattice has a
realisation with μ ∉ F_q, which could be free by §5.

**q = 8 (characteristic 2).** `ag8_4pencil.c`, `ag8_flex.c`, `ag8_orbits.c`.
- The search covers all 2³² offset choices. Directions are (1:0), (0:1), (1:1), (t:1), with
  t³ = t + 1; all 4-subsets of P¹(F₈) are equivalent.
- 250,006,624 configurations have d₁ ≥ 9:
  - K = 18: 143,661,952, exponents (9, 9);
  - K = 19: 100,259,712, exponents (9, 10);
  - K = 20: 5,750,080 with (9, 11) and 317,856 with (10, 10);
  - K = 21: 17,024, exponents (9, 12).
- They form 46,698 orbits under the 5,376 semilinear affine maps that preserve the four
  directions.
- Rigid lattices: 46,696 orbits have rank M(μ) = K − 2 at generic μ and at every μ in F₄, F₁₆,
  F₃₂, F₆₄ ∖ F₈, F₁₂₈ and F₂₅₆. Their only realisations have μ ∈ F₈ (μ of degree ≥ 9 over F₂ was
  not tested).
- **Flexible lattices.** 2 orbits (10,752 configurations, `flex_classes.txt`) deform to generic μ:
  - Both are 19-line (9, 9) lattices. Each is an *additive lift*: F₈ = F₂[t]/(t³ + t + 1) is
    replaced by the F₂-span of 1, μ, μ².
  - In one of them the horizontal offsets are {0, 1, μ, μ + 1, μ², μ² + 1, μ² + μ + 1} and the
    fourth pencil is x + μy = d with d ∈ {μ, μ², μ³ + μ² + μ, μ³} (`family_explicit.out`).
  - At generic μ the generic member has exactly the lattice of the F₈ hit: every concurrent triple
    agrees (126 and 122 triples respectively).
  - **The generic member is still nearly free, with syzygies (8, 11, 11)** (`flexcheck.out`). The
    degree-8 syzygy of the F₈ member deforms along the whole family. Lower semicontinuity of mdr
    then makes every member non-free, so TC holds for both lattices.
- No counterexample.

**q = 9 (characteristic 3).** `ag9_4pencil.c`.
- The search covers all 2³⁶ offset choices, with directions (1:0), (0:1), (1:1), (i:1) and i² = −1.
- 3,607,477,992 configurations have d₁ ≥ 10. They form 697,905 orbits under 64 · 81 semilinear
  affine maps.
- None is flexible: rank at a generic μ ∈ F₂₁₈₇ is K − 2.
- None has a realisation with μ in F₂₇, F₈₁ ∖ F₉, F₂₄₃ or F₇₂₉ ∖ F₉.
- No counterexample.

**Why the flexible ones fail.**
- The F_q member of a family through a configuration with d₁ > q is dense: by Lemma 7.1, K ≥ 2q + 2
  lines cover all but uv affine points.
- So every incidence it has must hold identically in μ. That forces F_p-linear structure in the
  offsets.
- Both kinds of flexible family met so far keep the degree-q syzygy along the family: the additive
  lifts here and the multiplicative three-pencil families at q = 11 (§2).

## 9. Where a positive-characteristic counterexample could still be

By §5 to §8 it needs one of the following:
- **A bad member that is not globally F_q-rational.** For four-pencil lattices this means a
  realisation space of dimension ≥ 2 whose F_q-fibre has only degenerate F_q-points, so it is
  invisible to the PG(2, q) census.
- **Lines with 5 or more points.** Five-pencil lattices over F₈ or F₉ need two cross-ratios in F_q.
  This is about 2⁴⁰ configurations, not run.
- **Prime fields q = p ≥ 11.** There are no additive lifts there, so flexibility would have to be of
  characteristic-0 type, with d₁ ≥ 12 and n ≥ 25. The four-pencil census over F₁₁ is about 2⁴⁴
  configurations, not run.
- **A free flexible family of characteristic-0 type with d₁ ≥ 9**, reduced mod 2 or 3. The only
  known flexible free families, ACKN-13 and ACKN-15, have d₁ = 6 and 7.
