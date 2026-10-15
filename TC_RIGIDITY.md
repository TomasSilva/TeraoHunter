# TC_RIGIDITY.md — rigidity against special incidence

Started 2026-10-14. The question: can Terao's conjecture in P² be proved by the dichotomy
"a free lattice is either too rigid to carry the special position that would make another
realisation non-free, or it is flexible, and then the conditions for a non-free realisation cannot
be met"?

Notation as in TC_ATTEMPT.md: L is a lattice of n lines with χ(t) = (t − 1)(t − d₁)(t − d₂),
d₁ ≤ d₂, δ = d₂ − d₁, so χ₀(L, 0) = b₂ − n + 1 = d₁d₂.

## 0. Verdict

- **The dichotomy is Yoshinaga's problem.** Yoshinaga [Y12, p. 20, (11)] asked whether every free
  arrangement is inductively free or rigid ("rigid": the moduli space is finite and all its points
  are Galois conjugate). A yes would imply TC, because inductive freeness is combinatorial and
  freeness is Galois-invariant.
- **Its literal form is false.** Abe, Cuntz, Kawanoue and Nozawa [ACKN] give a 13-line (6, 6)
  and a 15-line (7, 7) lattice. Both are free, are not recursively free for a generic parameter, and have a moduli space
  that is a rational curve. So a flexible free lattice need not be inductively free, and the flexible
  case has to be handled by locating where non-freeness can sit.
- **Refined form.** TC for L is equivalent to two statements.
  - (F) Freeness is constant on each connected component of M(L): no "Ziegler phenomenon" inside a
    family.
  - (R) Different connected components of M(L) have the same freeness. This is automatic when they
    are Galois conjugate.
  Yoshinaga's rigid case is the case where (F) is empty and (R) holds by Galois conjugacy. (F) and
  (R) together are exactly TC, so the dichotomy does not make the problem easier by itself. It
  localises it: (F) concerns special positions inside one family, (R) separate families.
- **Proved here** (§2), three lemmas about a degeneration from free to non-free with constant lattice:
  1. *Line saturation.* Non-freeness is decided by any single line, so every line must jump at once.
  2. *Boundary jump.* Every incidence that coarsens the lattice is a forced jump of every line it does
     not touch, of order at least ⌈√β⌉ (β = the drop of b₂). This explains half of the observation
     in TC_ATTEMPT §9 that the non-free locus "lands on the extra concurrencies": they are always in
     it. The other half, that nothing else is, is (F).
  3. *Collapse.* At an interior non-free point, the Saito determinant of the free basis vanishes to
     order ≥ 2. So the family would be tangent to the jump divisor of every line.
- **Computed** (§3), for both known flexible families and every line:
  - each line's jump function factors entirely into degeneration factors;
  - every remote degeneration occurs, at order exactly ⌈√β⌉;
  - at those points the arrangement is again free, with the exponents its new χ requires;
  - the special incidences inside the 13-line family (the set W of [ACKN], where four
    intersection points become collinear) leave every line at (6, 6).
- **Census** (§4). Up to 16 lines, every realisable lattice in the window of TC_ATTEMPT §10 has a
  moduli space that is a single Galois orbit of points or a single rational curve (the 15-line
  family). So (R) holds trivially there, and (F) holds by §3.
- **Not proved.** (R) in general: nothing prevents a moduli space with two components that are not
  Galois conjugate. (F) in general: nothing excludes a tangential interior jump. TC is not proved.

## 1. Setting

- M = M(L) is the moduli space of realisations of L modulo PGL(3), a quasi-affine variety over ℚ.
- F ⊆ M is the free locus. It is Zariski open [Yuz93]. N = M ∖ F is closed, and both are defined
  over ℚ and Galois-stable.
- For a line H of L with s_H points, Conf_H is the space of s_H labelled points on P¹ modulo
  PGL(2), carrying the multiplicities m^H(P) = m_P − 1. Let π_H : M → Conf_H be the map that takes
  a realisation to the configuration of points it induces on H.
- J_H ⊆ Conf_H is the jump locus, where the multiarrangement has e₁ < d₁. It is closed, by upper
  semicontinuity of dim D(·)_{d₁−1}.
- **Two types of counterexample.**
  - *Type I:* a free and a non-free realisation lie in different connected components of M.
  - *Type II:* some connected component meets both F and N. Since F is open and N closed, curve
    selection then gives a germ γ : (Δ, 0) → M with γ(t) free for t ≠ 0 and γ(0) non-free. This is
    the "Ziegler phenomenon": a special position inside a family with fixed combinatorics.

## 2. Three lemmas

**Lemma 1 (line saturation).** N = π_H⁻¹(J_H) for every line H of L.

*Proof.*
1. Yoshinaga [Yos05]: for every realisation A and every H ∈ A, e₁(H)e₂(H) ≤ χ₀(L, 0) = d₁d₂,
   with equality iff A is free.
2. Since e₁ + e₂ = n − 1 = d₁ + d₂ and e₁ ≤ e₂, equality holds iff e₁ = d₁.
3. The exponents of (A^H, m^H) depend only on π_H(A). ∎

*Consequences.*
- (a) A type-II germ must push the configuration of every line into its jump locus at the same
  parameter.
- (b) In a type-I counterexample, π_H(C₁) ∩ π_H(C₂) = ∅ for every line H, where C₁ is a free
  component and C₂ a non-free one. Every line separates the two components.
- (c) If for some line H the image π_H(M) is a single Galois orbit (a "rigid line"), TC holds for L.

**Lemma 2 (boundary jump).** Let A_t, t ∈ Δ, be a family of n distinct lines, with lattice L for
t ≠ 0. Suppose the limit A₀ is still n distinct lines but has a strictly coarser lattice L₀ (extra
concurrencies). Put β = b₂(L) − b₂(L₀). Let H be a line on which no two points of A_t come together.
Then:
- β ≥ 1;
- the Ziegler restriction of A₀ to H has the same multiplicities as that of A_t;
- its exponents are (d₁ − j, d₂ + j) with j(j + δ) ≥ β.

In particular π_H(A₀) ∈ J_H. If d₁ = d₂, the jump function Θ_H of §3 vanishes at t = 0 to order
at least j ≥ ⌈√β⌉.

*Proof.*
1. *β ≥ 1.* A point R of L₀ is the limit of points P₁, …, P_k of L. Every pair of lines through R
   meets in some P_i, so C(m_R, 2) = Σ C(m_{P_i}, 2) and m_{P_i} < m_R when k ≥ 2. Then
   Σ (m_{P_i} − 1) = Σ 2C(m_{P_i}, 2)/m_{P_i} > 2C(m_R, 2)/m_R = m_R − 1. Since the lattice
   coarsens, some R has k ≥ 2, so b₂ drops.
2. *Multiplicities on H are unchanged.* A line ℓ that is not generically through a point P of H,
   but passes through P in the limit, makes ℓ ∩ H collide with P. No collision on H means the
   multiplicities on H stay the same.
3. *The jump.* Apply Yoshinaga's inequality to A₀ and H:
   e₁e₂ ≤ χ₀(L₀, 0) = d₁d₂ − β. With e₁ = d₁ − j and e₂ = d₂ + j this reads j(j + δ) ≥ β.
4. *The order.* A square matrix depending on t has determinant vanishing to order at least the
   dimension of its kernel at t = 0. Here that kernel is D(A₀^H, m^H)_{d−1}, of dimension j. ∎

*Reading.* Every incidence that changes the lattice is a jump point of every line it does not
touch. The boundary of M is a forced sink for jumps.

**Lemma 3 (collapse).** Let γ be a type-II germ, R = ℂ{t} and A_R the family. Use
D₀ = {θ : θ(α_{H₀}) = 0} for a fixed line H₀. Choose θ₁, θ₂ ∈ D₀(A_R) of degrees d₁, d₂ as follows:
- θ₁ generates D₀(A_R)_{d₁};
- θ₂ generates D₀(A_R)_{d₂} modulo R[x]θ₁;
- if δ = 0, (θ₁, θ₂) is an R-basis of D₀(A_R)_d.

These degree parts are free R-modules, and reduction mod t is injective. Write
det(E, θ₁, θ₂) = c(t)·Q_t. Then:
- (a) θ₁(0) = fη and θ₂(0) = hη, where η ∈ D₀(A₀) spans a saturated line subsheaf O(−r′) ⊂ E₀
  (so η has isolated zeros), and f, h are forms with deg f ≥ 1. In particular V(f, h) ≠ ∅ in P².
- (b) ord_t c ≥ 2.
- (c) For every line H, the determinant of the restricted basis is c up to a unit. Hence, if
  d₁ = d₂ and the jump on H at t = 0 is j = 1, the jump function Θ_H vanishes to order ≥ 2 along γ:
  the curve π_H ∘ γ meets the jump divisor J_H with multiplicity ≥ 2, a tangency where J_H is smooth.

*Proof.*
1. *(a).* A₀ is not free, so c(0) = 0 (Saito's criterion), and θ₁(0), θ₂(0) are proportional over
   Frac(S).
   - Let λ ⊂ E₀ be the saturation of the subsheaf generated by θ₁(0). It is reflexive of rank one on
     P², so λ ≅ O(−r′), generated by some η ∈ D₀(A₀)_{r′}. Proportionality puts θ₂(0) in λ too, so
     θ₁(0) = fη and θ₂(0) = hη with f, h ∈ S.
   - If f were constant, θ₂(0) would lie in S_δ·θ₁(0), contradicting the choice of θ₂. When δ = 0
     it would contradict the linear independence of θ₁(0) and θ₂(0).
   - So deg f ≥ 1 and deg h = deg f + δ ≥ 1, and two forms of positive degree on P² have a common
     zero.
2. *Setting up (b).* Let E_R be the family of logarithmic bundles. It is flat: it is the syzygy
   bundle of the Jacobian ideals, whose colength τ is combinatorial. Its fibres are locally free, so
   E_R is locally free. The map F_R = O(−d₁) ⊕ O(−d₂) → E_R given by (θ₁, θ₂) is injective, and its
   cokernel Q lives on t = 0 with Fitt₀(Q) = (c).
3. *(b).* Suppose ord c = 1. Then tQ = 0. The Tor sequence of 0 → F_R → E_R → Q → 0 against R/t
   gives ker(F₀ → E₀) ≅ Q ≅ coker(F₀ → E₀). The kernel sits inside the bundle F₀, so it is
   torsion-free. The cokernel contains λ/(f, h)λ ≅ O_{V(f,h)}(−r′), which is nonzero torsion.
   Contradiction.
4. *(c).* For a line H, subtract multiples of E so that θᵢ(α_H) = 0. This does not change
   det(E, θ₁, θ₂). Then, with α_H = z and θᵢ = (aᵢ, bᵢ, 0):
   - a₁b₂ − a₂b₁ = c·Q/z;
   - restricting to H gives c·∏_P α_P^{m_P}, so the restricted determinant is c up to a unit;
   - a j = 1 jump on H is pulled back from the one-dimensional versal deformation of
     O(−d + 1) ⊕ O(−d − 1), where both Θ_H and that determinant vanish simply. ∎

*Reading.*
- On a single line, a free-to-non-free jump can be transversal; the versal family shows it.
- On P², it cannot: the free basis collapses onto one derivation η, with a common factor f, and
  the family must be tangent to the jump divisor of every line at once.
- At the remote boundary points of §3, the order is exactly ⌈√β⌉, which is transversal when β = 1.
- (b) is a statement about any family of rank-2 bundles on P² whose generic member splits. No
  arrangement-specific input is used beyond flatness.

## 3. The two flexible families, line by line

**Script.** `work/rigidity/linejumps.m2`. For each line H it computes the 2d × 2d matrix of
degree-(d − 1) derivations of the Ziegler restriction, over ℚ[l₀, l₁]. Θ_H is the gcd of its
determinants in the three coordinate charts. Each irreducible factor of Θ_H is classified as:
- *collision:* a new concurrency through H;
- *remote:* a new concurrency away from H;
- *interior:* neither, which would be a TC violation.

`specialise.m2` evaluates the family at a given algebraic parameter and reports the lattice, the
exponents on every line, and the lines outside A through ≥ 3 intersection points.

**13 lines, exponents (6, 6)** [ACKN §5, rescaled x → x/√3 so the family is over ℚ]: 21 double,
3 triple, 3 quadruple and 3 quintuple points (`fam13.m2`, `spec13.m2`). The parameter is λ = l₁/l₀.

- Degeneration parameters:
  - λ ∈ {0, 1, ∞} and λ² − λ + 1 = 0: lines coincide.
  - λ ∈ {−1, 2, ½}: one new triple point (the origin lies on H₂, H₆ or H₇); β = 1.

| lines | multiplicities on H | deg Θ_H | collision factors (order) | remote factors (order) |
|---|---|---:|---|---|
| 1–3 | 1,1,1,2,3,4 | 20 | λ=0 (5), ∞ (5), λ=1 (3), λ²−λ+1 (2), λ=−1 (1) | λ=2 (1), λ=½ (1) |
| 4–6 | 1,1,1,2,3,4 | 20 | λ=0 (3), ∞ (5), λ=1 (5), λ²−λ+1 (2), λ=2 (1) | λ=−1 (1), λ=½ (1) |
| 7–9 | 1,1,1,2,3,4 | 20 | λ=0 (5), ∞ (3), λ=1 (5), λ²−λ+1 (2), λ=½ (1) | λ=−1 (1), λ=2 (1) |
| 10–12 | 1,1,1,1,4,4 | 18 | λ=0 (3), ∞ (3), λ=1 (3), λ²−λ+1 (3) | λ=−1, 2, ½ (1 each) |
| 13 | 1,1,1,3,3,3 | 18 | λ=0 (3), ∞ (3), λ=1 (3), λ²−λ+1 (3) | λ=−1, 2, ½ (1 each) |

For H₁₃, the line at infinity, Θ_H = λ³(λ − 1)³(λ + 1)(λ − 2)(2λ − 1)(λ² − λ + 1)³. This matches
the 3 × 3 reduction in [ACKN §5].

- **No interior factor on any line.** Every remote degeneration divides every Θ_H at order
  ⌈√1⌉ = 1, as Lemma 2 requires.
- At λ ∈ {−1, 2, ½}: b₂ = 47 and χ₀ = 35, and all 13 lines have exponents (5, 7). So the boundary
  arrangement is free again (inductively free by [ACKN]).
- At λ = i and at λ = φ, the golden ratio (both in W): the lattice is unchanged (b₂ = 48), all 13
  lines are (6, 6), and new lines through 4 intersection points appear (generic λ: at most 3).
  These special incidences change the lattice of A ∪ H′ but not of A, and they do not move any
  exponent.

**15 lines, exponents (7, 7)** (TC_ATTEMPT Theorem 2; the same weak combinatorics as [ACKN §7],
and the only realisable (7, 7) window lattice with moduli, so the same lattice): 27 double, 6 triple
and 6 quintuple points (`fam15.m2`, `spec15.m2`).

- Degeneration parameters:
  - λ ∈ {0, ½, 1, ∞}: the six points Q_i become special.
  - λ² − 3λ + 1 = 0 and λ² + λ − 1 = 0: four new triple points; β = 4.

| lines | s_H | multiplicities | deg Θ_H | remote factors |
|---|---:|---|---:|---|
| 1–5, 8–11, 13–15 | 7 | 1,1,1,1,2,4,4 | 28–57 | none: the quadratics are collisions |
| 6, 7, 12 | 6 | 1,1,2,2,4,4 | 38, 42, 42 | λ²−3λ+1 and λ²+λ−1, order **2** each |

- **No interior factor on any line.** The remote order is 2 = ⌈√4⌉.
- At both quadratic parameters: b₂ = 59 and χ₀ = 45, and all 15 lines have exponents (5, 9). So the
  boundary arrangement is free.

**What the examples show.** In both families:
- the jump of every line is used up exactly by the boundary;
- incidences that keep the lattice (W) do not touch the derivation module;
- incidences that change it (Lemma 2) move every untouched line, by the least amount Yoshinaga's
  inequality allows, and land on a free arrangement again.

This is the user's mechanism, "a free arrangement cannot carry a freeness-breaking special
incidence without changing its combinatorics", in the only two flexible cases known.

## 4. The rigid side, and the census up to 16 lines

- **Rigid with one Galois orbit:** TC holds (Galois invariance). This is the whole content of
  "too rigid".
- **Several components that are not Galois conjugate:** not covered. Mnëv universality gives
  rank-3 matroids whose realisation spaces have several such components. I know no reason why a
  free lattice could not be one of them. Lemma 1(b) says every line would have to separate a free
  component from a non-free one.
- **Census.** The exhaustive window enumerations of TC_ATTEMPT §10 give:

| n | class | realisable window lattices | rigid, one orbit | with moduli |
|---|---|---:|---:|---:|
| 15 | (6, 8) | 1 (ℤ₅ net, Ceva over ℚ(ζ₅)) | 1 | 0 |
| 15 | (7, 7) | 10 (5 non-free, 5 free) | 9 | 1 (the 15-line family, §3) |
| 16 | (7, 8) | 0 | – | – |

  So up to 16 lines the refined dichotomy holds exactly in the window. Every realisable lattice is
  rigid with one Galois orbit, or is the one flexible family, where §3 shows the boundary absorbs
  every jump. Below 15 lines TC is known [BK23]. For at most 12 lines Yoshinaga's dichotomy even
  holds literally: every free arrangement is inductively free or one of three exceptional lattices
  realised over quadratic fields [ACKN, Theorem 1.1].

## 5. What is missing

- **(F) in general.**
  - What is needed: on a flexible component C, Θ_H ∘ π_H has no zero inside C. Equivalently,
    Θ_H restricted to M(L) is a product of incidence brackets [i j k], a unit on M(L).
  - By invariant theory, Θ_H is a polynomial in the brackets [H i j] alone. The remote factors
    [a b c] with H ∉ {a, b, c} enter only through the Grassmann–Plücker relations of L.
  - Lemma 2 accounts for the boundary factors. Lemma 3 says an interior zero would be a tangency of
    order ≥ 2 on every line simultaneously.
  - Nothing here excludes such a tangency. A proof would need a global count: the degree of Θ_H on a
    compactified component equals the forced boundary contributions. I have a formula for neither
    side. The degrees in §3 (18–57) are matched by the boundary in each case, but I cannot predict
    them.
- **(R) in general.** No argument beyond Galois conjugacy and Lemma 1(b).
- **TC is not proved.** What is new is the reduction to (R) and (F), the three lemmas, and the
  line-by-line verification of the mechanism in the two flexible families.

## 6. Files

- `work/rigidity/linejumps.m2`: the per-line jump analysis.
- `work/rigidity/fam13.m2`, `fam15.m2`: the two families, with outputs `fam13.out` and `fam15.out`
  (about 8 s and 22 s).
- `work/rigidity/specialise.m2`, `spec13.m2`, `spec15.m2`: specialisations at the boundary and at W,
  with outputs `spec13.out` and `spec15.out`.
- `work/rigidity/ackn15check.m2`: the published 15-line family has the same weak combinatorics.

## References

- [ACKN] T. Abe, M. Cuntz, H. Kawanoue, T. Nozawa, *Non-recursive freeness and non-rigidity*,
  arXiv:1411.3351.
- [BK23] M. Barakat, L. Kühne, *Computing the nonfree locus of the moduli space of arrangements and
  Terao's freeness conjecture*, Math. Comp. 92 (2023), 1413–1452.
- [Y12] M. Yoshinaga, *Freeness of hyperplane arrangements and related topics*, arXiv:1212.3523.
- [Yos05] M. Yoshinaga, *On the freeness of 3-arrangements*, Bull. London Math. Soc. 37 (2005),
  126–134.
- [Yuz93] S. Yuzvinsky, *Free and locally free arrangements with a given intersection lattice*,
  Proc. Amer. Math. Soc. 118 (1993), 745–752.
