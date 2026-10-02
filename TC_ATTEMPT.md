# Terao's conjecture: what a counterexample must satisfy, and what can be proved

**Bottom line.** I could not prove Terao's conjecture (TC) in P² or P³. It remains open, and
nothing here closes it. What I can do:

1. Pin down what a P² counterexample would have to look like. Most conditions are in the
   literature; one refinement (Theorem 1(7)) and one principle (Theorem 1(8)) I prove here (§1).
2. Prove TC exactly, on the full complex realisation space, for the most dangerous lattices my
   searches found:
   - a free, non-divisionally-free 15-line lattice with moduli, beyond the 14-line range where
     TC is verified (§2);
   - a free, non-DF 26-plane arrangement in P³ (§5).
3. Narrow the smallest open size, 15 lines, to two exponent types, and prove TC for one of the
   four profiles that remain in the first type (§3).
4. Prove TC over finite fields in two cases (§4):
   - all line arrangements over F_q with q ≤ 7 (a direct combination of known theorems);
   - "covering" lattices in every rank (rank ≥ 4 is new).
5. Reduce P³ to P² plus one statement about multirestrictions (§5).
6. Say exactly where every argument stops (§6).

Notation (P²). L is a rank-3 lattice with n lines. A₁ is a free realisation with exponents
(1, d₁, d₂), d₁ ≤ d₂, δ = d₂ − d₁. A₂ is a non-free realisation of the same L. Both are over ℂ
(any field of characteristic 0). For a line H, s_H = |A^H| is the number of multiple points on H,
and m is the largest point multiplicity.

## 1. Anatomy of a counterexample in P²

**Theorem 1.** If (A₁, A₂) is a counterexample, then:

1. n ≥ 15. [Barakat–Kühne 2023: TC holds for ≤ 14 lines, any characteristic.]
2. L is not divisionally free: s_H ∉ {d₁ + 1, d₂ + 1} for every H. [Abe 2016: DF lattices
   satisfy TC; in rank 3, DF ⟺ some s_H − 1 is a root of χ.]
3. s_H ≤ d₁ for every line H. [Abe 2014, Cor. 1.2: a free line arrangement has s_H ≤ d₁ + 1
   or s_H = d₂ + 1; then apply (2).]
4. No point has multiplicity > (n + 1)/2, and s_H ≥ δ + 4 for every line H. [Abe 2012:
   Prop. 4.2 handles an unbalanced Ziegler restriction, i.e. such a heavy point; for balanced
   restrictions, Thm. 0.1 gives |e₁ − e₂| ≤ s − 2 and Cor. 0.4 the bound on s_H.]
5. Hence δ + 4 ≤ d₁, i.e. d₂ ≤ 2d₁ − 4, and d₁ ≥ 6. [Abe 2014, Cor. 5.4 and 5.5.]
6. m ≤ d₁ − 1. If d₁ = m + 1, then A₂ is plus-one generated with exponents
   (m, n − m, 2n − 3m − 3). [Dimca–Kühne–Pokora 2025: TC holds when d₁ ∈ {m − 1, m} or
   d₁ = n − m; their Thm. 2.10 gives the d₁ = m + 1 case.]
7. **(refinement)** A₂ has a logarithmic derivation of degree r = mdr(A₂) with m ≤ r ≤ d₁ − 1.
   Every line satisfies
     n + 1 − 2r ≤ s_H ≤ r + 1,
   and every Ziegler restriction of A₂ has exponents (e₁, e₂) with s_H − 1 ≤ e₁ ≤ r < d₁.
8. **(rigidity principle)** Suppose that for some line H, the points of A^H with their
   multiplicities are the same up to PGL(2) and field automorphisms for all realisations of L —
   for example s_H ≤ 3, or the moduli of L do not move the points of H. Then TC holds for L. So in
   a counterexample the moduli must move the points on *every* line.

*Proof of (7).* Write mdr for the least degree of a derivation that is not a multiple of E.
1. *m ≤ r ≤ d₁ − 1.* This is in Dimca–Kühne–Pokora. It comes from the du Plessis–Wall bounds:
   τ(A) = Σ_P (m_P − 1)² is combinatorial, and among arrangements with this τ the free ones are
   exactly those with mdr = d₁, while a non-free one must have mdr < d₁.
2. *e₁ ≤ r.* Fix H. Write D₀^H(A₂) for the derivations with θ(α_H) = 0; it complements S·E, so
   its least degree is r. An element of degree r restricts to a nonzero element of the Ziegler
   multirestriction (A₂^H, μ), where μ_P = m_P − 1. Indeed, the kernel of restriction is
   α_H·D₀^H(A₂), which lives in degrees ≥ r + 1. So e₁ ≤ r.
3. *s_H ≤ e₁ + 1.* For any 2-multiarrangement with s points and multiplicities μ, either
   s ≤ e₁ + 1 or s = e₂ + 1:
   - if the lowest basis element θ₁ is not a multiple of E_H, then det(θ₁, E_H) is a nonzero form
     of degree e₁ + 1 vanishing at all s points;
   - otherwise θ₁ = gE_H with α_P^{μ_P − 1} | g. Then e₁ − 1 ≥ Σ_P (μ_P − 1) = (n − 1) − s, so
     e₂ ≤ s − 1, while θ₂ (not a multiple of E_H) gives s ≤ e₂ + 1.
   The case s_H = e₂ + 1 > d₁ contradicts (3). So s_H ≤ e₁ + 1 ≤ r + 1.
4. *s_H ≥ n + 1 − 2r.* Abe's bound e₂ − e₁ ≤ s_H − 2 gives s_H ≥ n + 1 − 2e₁ ≥ n + 1 − 2r. ∎

*Proof of (8).*
- Yoshinaga's criterion: A is free iff the Ziegler multirestriction onto H has exponents
  (d₁, d₂).
- If the multirestrictions of A₁ and σ(A₂) are isomorphic, where σ is a field automorphism, they
  have the same exponents, so σ(A₂) is free. Freeness is Galois-invariant, so A₂ is free.
- This is Abe–Yoshinaga, Cor. 5.3, in rank 3. ∎

**The smallest admissible exponents** (from (1), (4), (5)):

| n | possible (d₁, d₂) | allowed s_H |
|---|---|---|
| 15 | (6, 8), (7, 7) | {6} and [4, 7] |
| 16 | (7, 8) | [5, 7] |
| 17 | (7, 9), (8, 8) | [6, 7] and [4, 8] |
| 18 | (7, 10), (8, 9) | {7} and [5, 8] |
| 19 | (8, 10), (9, 9) | [6, 8] and [4, 9] |
| 20 | (8, 11), (9, 10) | [7, 8] and [5, 9] |

## 2. Exact proofs of TC for the dangerous lattices

**Method** (as in Barakat–Kühne):
1. Write the complete realisation space over ℚ̄ as an ideal: fix a projective frame, impose the
   incidences, saturate away every degeneration (extra collinearity or concurrency), and decompose
   over ℚ.
2. On each component, decide freeness exactly. By Theorem 1(7), a realisation of L is free iff it
   has no derivation of degree d₁ − 1 apart from multiples of E. This is a full-rank condition on
   a matrix M over the coordinate ring.
3. The non-free locus is the support of coker(Mᵀ). Compare it with the degeneration locus.

Scripts: `work/tc15_exact.m2`; the general checker is `work/exact_tc_points.py`.

**Theorem 2 (new).** TC holds for the 15-line lattice L₁₅: free with exponents (1, 7, 7), not DF,
with 6 quintuple, 6 triple and 27 double points. This is the "database lattice" of Findings 9 and
15 in `NOTES.md`.

*Proof.*
1. Up to PGL(3) the realisations are the 15 lines joining Q₁, …, Q₆ with Q₁…Q₄ the standard
   frame, Q₅ = (d² : d : 2d − 1) and Q₆ = (d(2d − 1) : d² : 2d − 1). Here d ranges over ℂ minus the
   roots of
     B(d) = d²³ (d − 1)²³ (2d − 1)⁷ (d² − 3d + 1)⁴ (d² + d − 1)⁴,
   which are exactly the parameters where the lattice degenerates. The realisation ideal has one
   non-degenerate component; it is a rational curve.
2. With θ(z) = 0, the degree-6 derivations form the kernel of a 98 × 56 matrix over ℚ[d] of
   generic rank 56, so the generic member is free.
3. The annihilator of coker(Mᵀ) is generated by
     d⁷ (d − 1)⁷ (2d − 1)³ (d² − 3d + 1)(d² + d − 1).
4. Every root of it is a root of B. So every realisation of L₁₅ is free. ∎

This lattice has d₁ = m + 2 (ε = 2 in Dimca–Kühne–Pokora's notation), the class they identify
as having the most room for counterexamples. The same computation for the 13-line lattice with
moduli (exponents (1, 6, 6)) also passes; it is within the 14-line range anyway.

**Other candidates** (`work/exact_tc_driver.py`). I ran the same exact check on every distinct
free, non-DF lattice from my connecting-line searches over F₁₁, F₁₃ and F₁₇ (56 lattices,
13–28 lines):
- **13-line (1, 6, 6) lattice with moduli:** passes (one rational curve, non-free locus inside
  the degeneration locus). It is within the 14-line range anyway.
- **15-line lattice:** independently confirms Theorem 2.
- **One 19-line (1, 9, 9) lattice from F₁₁:** its saturated realisation ideal is (1), so it has
  no realisation in characteristic 0.
- **The rest are still running** (slowly; the machine was under memory pressure). Results
  append to `work/exact_tc_driver.out`.

Most of these lattices have d₁ = q (the field size). That is the signature of the Frobenius
derivation, so I expect them to be finite-field artefacts. The exception is the 21-line
(1, 9, 11) lattice on 7 points. It occurs over F₁₁ but not over F₁₃ or F₁₇, which is the
pattern of a realisation over ℚ(√−7). It is the one to watch.

## 3. Fifteen lines

**Proposition 3.** A 15-line counterexample has exponents (6, 8) or (7, 7).

**The (6, 8) case.**
- Every line carries exactly 6 multiple points, and m ≤ 5.
- Counting (Σ_P C(m_P, 2) = 105, Σ_P (m_P − 1) = 62, Σ_H s_H = 90) leaves four point-multiplicity
  profiles (n₂, n₃, n₄, n₅):
    (3, 16, 9, 0), (2, 19, 6, 1), (1, 22, 3, 2), (0, 25, 0, 3).

**Proposition 4.** A line arrangement with profile (0, 25, 0, 3) is projectively equivalent to the
Ceva arrangement (x⁵ − y⁵)(y⁵ − z⁵)(z⁵ − x⁵) = 0. So TC holds for it.

*Proof.*
1. Every line has one quintuple point and five triple points. So the lines split into three
   pencils of 5 lines through P₁, P₂, P₃, and no line contains two of the P_i.
2. A line of pencil 1 meets the 10 lines of pencils 2 and 3 at its 5 triple points. So each
   intersection of lines from different pencils is a triple point: the arrangement is a 3-net of
   order 5.
3. *If P₁, P₂, P₃ are collinear,* send their line to infinity. The pencils become parallel
   classes x = aᵢ, y = bⱼ, x + y = c_k with aᵢ + bⱼ ∈ {c_k}. So |A + B| = 5 for 5-element sets
   A, B ⊂ ℂ, contradicting |A + B| ≥ |A| + |B| − 1 in a torsion-free group.
4. *So the P_i are the coordinate points.* The pencils are x = aᵢy, y = bⱼz, z = c_kx with
   aᵢ, bⱼ, c_k ∈ ℂ*, and the triple points force aᵢbⱼ ∈ {c_k⁻¹}. So |A·B| = 5 = |A| = |B|.
5. By Kneser's theorem, A and B are cosets of the subgroup μ₅ ⊂ ℂ*.
6. Rescaling coordinates turns A, B, C into μ₅. That is the Ceva arrangement, which is free
   (the reflection arrangement of G(5, 5, 3), exponents (1, 6, 8)). ∎

**Still open.** The profiles (3, 16, 9, 0), (2, 19, 6, 1) and (1, 22, 3, 2), and the (7, 7) case,
are the smallest cases these constraints leave open. For the two profiles with m = 5 = d₁ − 1,
Dimca–Kühne–Pokora add one condition: a non-free realisation would be plus-one generated with
exponents (5, 10, 12). Settling 15 lines would need a Barakat–Kühne-type enumeration of these
profiles.

## 4. Finite fields

**Proposition 5.** TC holds for every line arrangement over F_q with q ≤ 7.

*Proof.* Let n be the number of lines.
- *n ≤ 14:* Barakat–Kühne, valid in every characteristic.
- *n ≥ 2q:* Abe 2014, Cor. 6.2, building on Yoshinaga 2006. If χ(A, q) = 0 then A is free, and
  if n ≥ 2q and A is free then χ(A, q) = 0. So freeness is equivalent to χ(A, q) = 0, which is
  combinatorial.
- For q ≤ 7 we have 2q ≤ 14, so the two ranges cover every n. ∎

For q = 8 the open range is n = 15; for q = 9 it is 15 ≤ n ≤ 17. My exhaustive PG(2, 5) search
(Finding 13) agrees with this. As far as I know this combination is not stated anywhere, but it
is immediate from published results.

**Proposition 6 (all ranks).** Over F_q, TC holds for every *covering* lattice. Here covering means
that every codimension-2 subspace lies in some hyperplane, equivalently
|A|·[ℓ−1]_q − b₂ = [ℓ choose 2]_q, which is combinatorial.
- Such arrangements are free with exponents (1, q, …, q^{ℓ−2}, |A| − [ℓ−1]_q), by Lemma 1 of
  `FREE_NOT_DF.md`.
- In rank 3 this is Yoshinaga's theorem.
- Over F₂ and F₃, TC holds trivially, because binary and ternary matroids are uniquely
  representable.

## 5. P³

**Proposition 7 (reduction).** Let (A₁, A₂) be a rank-4 counterexample. Then one of two things
holds:
- **(a)** Some rank-3 localisation of L is itself a TC counterexample, so the problem drops to P².
- **(b)** A₂ is locally free in codimension 3. Then for every plane H the Ziegler
  multirestriction (A₂^H, m^H) is a non-free rank-3 multiarrangement, while (A₁^H, m^H) is free
  [Abe–Yoshinaga, Thm. 5.1(3)]. The two share their simple lattice and multiplicities, and also
  every local exponent, since those come from the free localisations. They are never isomorphic
  [Abe–Yoshinaga, Cor. 5.3].

In addition:
- L is not DF [Abe 2016].
- L has no locally heavy flag [Abe–Kühne].

So TC in P³ follows from TC in P² plus a statement about multirestrictions. That statement is:
two locally free rank-4 realisations whose Ziegler multirestrictions share all combinatorial and
local data have the same multirestriction freeness.

**What my experiments say about that statement.**
- It cannot be proved from the multiarrangement alone. Freeness of rank-3 multiarrangements is
  not combinatorial: in the X3 family it is non-free exactly at t = 1, or at t² − t + 1 = 0
  (Finding 1).
- In every free lift I built (more than 2,300), the rank-4 lattice changed exactly where the
  multirestriction became non-free, and the lift was DF (Findings 2, 7, 8).

**Corollary 8.** If some plane H has a restriction lattice that is projectively unique up to field
automorphisms (with its multiplicities), then TC holds for L.

Known free non-DF rank-4 examples:
- G(r, r, 4): over the finite fields I tried, its realisations were few, Galois-type and all free
  (Finding 14). I have not computed its complex realisation space.
- **The 26-plane arrangement** (Finding 16; exponents (1, 7, 9, 9), not DF). TC holds for it
  (`work/exact26_space.py`), because its complex realisations form one Galois orbit:
  1. Its lattice determines the 8 points that span the planes, and their collinearity and
     coplanarity structure. So every realisation is, up to PGL(4), a member of the 3-parameter
     family of Finding 16.
  2. Imposing all 8 + 47 + 124 dependencies of the lattice and saturating away the 54 forbidden
     coplanarities gives the realisation ideal (a − c, b − c, c² − c + 1).
  3. So the only complex realisations are A(ω, ω, ω) and its complex conjugate.

  Its freeness was checked over F₇, F₁₃ and F₁₉. TC for it does not depend on that.
- The seven F₃ examples are rigid and exist only in characteristic 3 (`FREE_NOT_DF.md`).

## 6. Where the proof stops

**In P².**
- Every constraint in Theorem 1 is numerical: points per line, exponents, multiplicities. They
  rule out most lattices, but they contradict nothing in general.
- A counterexample is a lattice inside the "window" of Theorem 1(3)–(7). On its realisation space
  the logarithmic bundle splits generically but jumps on a nonempty closed set, to an unstable
  bundle with a section of degree r < d₁, while the lattice stays the same.
- In every case I could compute exactly (§2), the jump locus is contained in the degeneration
  locus. The bundle "cannot jump without an extra incidence."
- A proof would need a reason why a jump of the logarithmic bundle forces an extra concurrency.
  For a lattice with irreducible realisation space, that is equivalent to TC for the lattice.
  I see no mechanism that gives it in general.
- Dimca's Conjecture 1.2 (mdr < (n − 2)/2 is combinatorial) would imply TC. Its stronger form
  fails for a 9-line pair (Dimca–Pokora 2026).

**In P³.** The missing piece is the multirestriction statement in Proposition 7. The experiments
show it holds whenever I could test it, through "protection": the lift's lattice records the
special parameter that the restriction's lattice cannot see. But I have no proof of that
mechanism.

## Sources

- Abe, [Exponents of 2-multiarrangements and freeness of 3-arrangements](https://arxiv.org/abs/1005.5276) (2010/2012).
- Abe, [Roots of characteristic polynomials and intersection points of line arrangements](https://arxiv.org/abs/1302.3822), J. Singularities 8 (2014).
- Abe, [Divisionally free arrangements of hyperplanes](https://arxiv.org/abs/1502.07520), Invent. Math. 204 (2016).
- Abe–Yoshinaga, [Free arrangements and coefficients of characteristic polynomials](https://arxiv.org/abs/1109.0668), Math. Z. (2013).
- Barakat–Kühne, [Computing the nonfree locus of the moduli space of arrangements and Terao's freeness conjecture](https://arxiv.org/abs/2112.13065), Math. Comp. 92 (2023).
- Dimca–Kühne–Pokora, [Free line arrangements with low maximal multiplicity](https://arxiv.org/abs/2505.01733).
- Dimca–Pokora, [A nine-line counterexample to a conjecture on the minimal degree of Jacobian relations](https://arxiv.org/abs/2607.01985).
- Yoshinaga, On the freeness of 3-arrangements, Bull. LMS 37 (2005); Free arrangements over finite field, Proc. Japan Acad. 82 (2006).
