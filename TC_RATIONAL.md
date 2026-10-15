# Terao's conjecture over ℚ

**Status (15 Oct 2026).** Not proved. This note does four things:
- it makes precise how TC over ℚ sits among the other versions;
- it isolates the one feature of ℚ that the other fields lack: only *rational* special positions matter;
- it turns that feature into theorems for lines with four points;
- it maps out the rational special positions for lines with five points.

Scripts and data are in `work/rational/`.

Notation:
- A lattice L of n lines has χ(L, t) = (t − 1)(t − d₁)(t − d₂), with d₁ ≤ d₂. If χ does not split, no
  realisation is free and there is nothing to prove.
- For a line H, its Ziegler multiplicities are k_X = m_X − 1 for the points X ∈ H, with
  K = Σ k_X = n − 1.
- e₁(H) ≤ e₂(H) are the exponents of the multirestriction (A^H, k).

## 1. Where TC over ℚ sits

TC(F) says: two arrangements defined over F with the same lattice are both free or both
non-free. Freeness is unchanged by field extension, which gives two chains of implications:
- TC(ℂ) ⟹ TC(ℚ̄) ⟹ TC(F) for every number field F ⟹ TC(ℚ);
- TC(F_p) for infinitely many p ⟹ TC(ℚ) (Proposition 1.1).

None of the converses is known, so TC(ℚ) is the weakest version.

**Proposition 1.1 (spreading out).** Let A and B be arrangements over ℚ with lattice L, A free
and B not free. Then there is a finite set S of primes such that, for every p ∉ S, the
reductions A_p and B_p (in fixed integral models) have lattice L, A_p is free and B_p is not free.
In particular TC(F_p) fails for all p ∉ S.

*Proof.*
- **Lattice.** The lattice is preserved mod p unless p divides one of the finitely many nonzero
  3 × 3 determinants of non-concurrent triples.
- **A_p is free.** Choose a Saito basis θ_E, θ₁, θ₂ of D(A) with primitive integral coefficients.
  Then det = c·Q with c ∈ ℤ ∖ 0, and for p ∤ c Saito's criterion over F_p makes A_p free.
- **B_p is not free.** By du Plessis–Wall, a non-free B has mdr(B) = r < d₁: the bound
  τ ≤ (n−1)² − r(n−1−r), with equality iff free, applies with τ = (n−1)² − d₁d₂ combinatorial.
  Take a primitive integral θ ∈ D₀(B) of degree r. For p ∤ n its reduction is a nonzero derivation
  killing Q_p, hence not a multiple of θ_E. A free B_p would have only multiples of θ_E in degrees
  < d₁. ∎

So TC over ℚ would follow from TC over F_p for infinitely many p. The proposition cannot be
used computationally: S depends on the heights of A and B, so no finite list of primes is safe.
This is why the censuses of TC_BUILD.md and TC_CHARP.md can never prove TC(ℚ) by themselves.

**Why ℚ is easier.** For large p, the F_p-points meet every geometrically irreducible piece of the
non-free locus that is defined over F_p (Lang–Weil). That includes pieces with no rational points
at all. TC over ℚ only sees the rational points of the non-free locus.

So suppose a free family had all its non-free members at irrational special positions, say
equianharmonic or in ℚ(√21). That would be a counterexample over ℂ but not over ℚ. The rest of
this note exploits exactly this gap.

## 2. Freeness is read off one line

**Yoshinaga's criterion (rank 3, characteristic 0).** A realisation is free iff e₁(H) = d₁ for
one (equivalently every) line H. Otherwise e₁(H) < d₁ on every line, since e₁e₂ ≤ d₁d₂ always.

**Abe's bound.** For a balanced multiplicity (2·max k ≤ K) on s points of P¹, |e₂ − e₁| ≤ s − 2.
The generic value is ⌊K/2⌋.

Let J(k) ⊂ M_{0,s} be the *jump locus*: the positions of the s points where e₁ drops below its
generic value. A non-free realisation has every line in its jump locus. This gives:

**Lemma 2.1.** If some line H of L has a jump locus J(k_H) with no non-degenerate rational point,
then all realisations of L over ℚ have the same freeness.

Lines with ≤ 3 points have no moduli, and freeness through them is combinatorial (classical). The
first interesting case is s = 4, where J(k) is a finite set of cross-ratios μ, the zeros of an
integer polynomial.

## 3. Four-point lines: the rational jump points are harmonic

`jump4.m2` computes, for each multiplicity vector k = (k₁, k₂, k₃, k₄) on the points x, y, x − y,
x − μy:
- the matrix of the conditions on a degree-e derivation;
- the gcd of its maximal minors in ℚ[μ] at e = e₁^gen − 1, whose roots are J(k).

It was run for **all 1,005 balanced vectors with K ≤ 30** (`jump4_all.out`).

**Validation.**
- J(4,4,4,4) = roots of μ² − μ + 1 (equianharmonic: E₁₇, nearly free) and of μ² − 5μ + 1,
  μ² + 3μ − 3, 3μ² − 3μ − 1 (one anharmonic orbit in ℚ(√21)).
- This agrees with TC_FILTER.md Lemma F1, as do the other 17-line vectors.

**Theorem 3.1 (computational, K ≤ 30).** For every balanced 4-point multiplicity vector:
- (a) If K is odd, J(k) = ∅: the exponents are ((K−1)/2, (K+1)/2) at every non-degenerate
  cross-ratio. This also follows from Abe's bound.
- (b) If K is even, every rational point of J(k) lies in the harmonic orbit {−1, 2, ½}.
- (c) Rational jump points occur for exactly 73 vectors (`jump4_rational_patterns.txt`): 2, 5, 8,
  13, 18, 25 at K = 8, 12, 16, 20, 24, 28, none at K = 10, 14, 18, 22, 26, and two at K = 30,
  namely (2, 8, 9, 11) and (4, 6, 7, 13).

All the other jump points are irrational, e.g. μ² − μ + 1, μ² − 7μ + 1, μ² + 5μ + 1,
3μ² − 5μ + 3, and irreducible factors up to degree 24.

For K = 16 the eight vectors with a rational (harmonic) jump are:
- (1,1,7,7), (1,3,6,6), (1,4,4,7), (1,5,5,5);
- (2,2,5,7), (3,3,3,7), (3,3,5,5), (3,4,4,5).

**Theorem 3.2 (TC over ℚ with a 4-point line).** Let L be a lattice with n ≤ 31 lines, χ split,
and a line H with exactly four points.
- (a) If n is even, or the multiplicity on H is unbalanced, then all realisations of L in
  characteristic 0 have the same freeness.
- (b) If n ∈ {11, 15, 19, 23, 27}, all realisations of L over ℚ have the same freeness. So TC(ℚ)
  holds for L.
- (c) For the other odd n ≤ 31, TC(ℚ) can fail for L only if all of the following hold:
  - d₁ = d₂ = (n − 1)/2;
  - the multiplicity vector of H is one of the 73 above;
  - the non-free rational member has H at the jumping harmonic position.

  The same then holds on **every** 4-point line of L.

*Proof.*
- **(a), unbalanced H.** The exponents on H are (K − k_max, k_max) at every position.
- **(a), n even.** By Theorem 3.1(a), e₁(H) = (K − 1)/2 at every realisation, so by Yoshinaga
  either every realisation is free or none is.
- **(b) and (c), K = n − 1 even.** By Abe's bound e₁(H) ∈ {K/2 − 1, K/2}. Since e₁e₂ ≤ d₁d₂, this
  forces d₁ ∈ {K/2 − 1, K/2}.
  - If d₁ = K/2 − 1: every realisation needs a jump, has e₁ = d₁, and is free.
  - If d₁ = K/2: a realisation is non-free iff its cross-ratio on H lies in J(k).

  Over ℚ, such a cross-ratio is rational, so by Theorem 3.1(b) and (c) it is harmonic and k is on
  the list. Non-freeness shows on every line, which gives the last sentence. ∎

Part (b) is genuinely arithmetic. Over ℂ the same lattices may have non-free realisations at
irrational jump points; vector (2,2,3,3) already has the jump points μ² − 7μ + 1 = 0. Over ℚ these
positions do not exist.

## 4. Five-point lines (the 17-line window vectors)

For s = 5 and K even, J(k) is a plane curve F_k(a, b) = 0 in the two cross-ratios (`jump5.m2`).
For each of the eleven vectors of the (8, 8) window, the non-degenerate part:
- is irreducible over ℚ, of degree 8 to 11;
- has total Tjurina number 12 to 34, against an arithmetic genus of 21 to 45 (`sing5.out`).

At every plane-curve singularity δ ≤ τ. This follows from μ = 2δ − r + 1, μ ≥ (r − 1)², and the
known lower bound τ > 3μ/4 for plane curve germs; for a node, δ = τ = 1. So the geometric genus is
at least 3, assuming the curve is irreducible over ℚ̄. If it is not, its rational points lie in
the finite intersection of the conjugate components. **So each curve has finitely many rational points (Faltings).**

**Validation (`validate.out`).** The nearly free rational arrangements R₁₇ᵃ and R₁₇ᵇ lie exactly on
their curves, at cross-ratios (4/3, 2) and (−½, −1); random controls do not.

**Rational points found** (`jump5_ratpts_H40.out`): every rational point with a of height ≤ 40, b
unrestricted. Raising the bound from 12 to 40 found nothing new.

| vector | points |
|---|---|
| (5,5,4,1,1) | the six points with a, b ∈ {−1, ½, 2} |
| (5,5,3,2,1) | (½, 2/3), (2, 3/2) |
| (5,5,2,2,2) | twelve points with a, b ∈ {±2, ±1, ±½} |
| (5,4,4,2,1) | (2, −2), (2, 2/3) |
| (5,4,3,3,1) | (½, −1), (½, ⅓), (2, −2), (2, 2/3) |
| (5,4,3,2,2) | (−1, 3), (3, −1) |
| (4,4,4,3,1) | (−1, ⅓), (−1, 3), (½, −½), (½, 3/2), (2, −2), (2, 2/3) |
| (5,3,3,3,2), (4,4,4,2,2), (4,4,3,3,2), (4,3,3,3,3) | none |

All are small-height "harmonic-net" configurations. The curves do have real points, so there is
no real obstruction (`real5.m2`), and proving the list complete would need Chabauty-type methods.

## 5. Consequences for 17 lines over ℚ

TC over ℂ is verified for n ≤ 16, so TC(ℚ) holds there too. At n = 17 only the class (8, 8)[B] is
open; in it the non-free member A₂ is nearly free (7, 10). Over ℚ:
- **Settled lattices.** Every [B] lattice with a 4-point line of vector (4,4,4,4) or (5,5,4,2) has
  all its rational realisations free (Theorem 3.2 with Lemma 2.1), so TC(ℚ) holds for it. Unlike
  the census of TC_BUILD.md §5, this does not need good reduction at 11 or 13.
- **4-point lines of a rational A₂.** Each one sits at a harmonic position:
  - (5,5,5,1): t ∈ {−1, 2, ½};
  - (5,5,3,3): t = −1;
  - (5,4,4,3): t = ½ (in the ordering of TC_FILTER.md Lemma F1).
- **5-point lines of a rational A₂.** Each lies on one of finitely many rational points of its
  curve. If the list in §4 is complete, A₂ has no 5-point line of vector (5,3,3,3,2), (4,4,4,2,2),
  (4,4,3,3,2) or (4,3,3,3,3). Such lattices would then have only free rational realisations.
- **Jumping point.** As over ℂ, the jumping point P₀ of A₂ is rational and lies on no 4- or 5-point
  line (TC_FILTER.md). Over ℚ it also lies on no 6-point line (`TC_GLOBAL.md` §5, Step 2), so every
  line of A₂ through P₀ has exactly 7 points.

## 6. Why this is not a proof of TC(ℚ)

- **Rational special positions exist and are realised.** R₁₇ᵃ and R₁₇ᵇ are rational, non-free and
  rigid, and sit at rational jump points of their 5-point lines. Harmonic positions are rational.
  So arithmetic alone cannot exclude a non-free rational member. It only confines it to finitely
  many special configurations on its small lines.
- **Lines with ≥ 6 points.** The jump loci are hypersurfaces of dimension ≥ 2 in M_{0,s}
  (TC_FILTER.md Lemma F1), and nothing confines their rational points.
- **The harmonic case.** Theorem 3.2(c) for a 4-point line needs exactly what TC over ℂ needs: that a
  family whose generic member is free never meets the harmonic position non-degenerately, or that
  no other component of the realisation space does.

What the method does give: **over ℚ, the "special incidence" of a counterexample can only be
harmonic on 4-point lines, and one of finitely many explicit configurations on 5-point lines.**
Over ℂ it may also be equianharmonic, golden or one of many other irrational positions.

## 7. Files (`work/rational/`)

- `jump4.m2`, `jump4_run.m2`, `jump4_range.m2`: the 4-point jump loci.
  - Output: `jump4_all.out`, all 1,005 vectors, K ≤ 30, with factorisations.
  - The 73 rational-jump vectors: `jump4_rational_patterns.txt`.
- `jump5.m2`, `jump5_cli.m2`, `jump5_ratpts.m2`: the 5-point curves and their rational points.
  - Outputs: `jump5_ratpts_H40.out`, `sing5.m2`/`sing5.out`, `real5.m2`.
- `validate.m2` / `validate.out`: the checks against R₁₇ᵃ and R₁₇ᵇ.
