# Terao's conjecture in P² for all n: an attempt at a proof from the computed patterns

**Status (8 Oct 2026).** This is not a proof of Terao's conjecture (TC). It is the result of
trying to turn the patterns from the exhaustive 15-, 16- and (running) 17-line verifications
into a general argument. What came out:

1. a single clean statement, **Conjecture J** (§4), that implies TC in P² for every n, together
   with a proof of that implication (Theorem 4.2) and a proof that J is *equivalent* to TC
   (Theorem 4.3) — so the reformulation loses nothing and isolates exactly what is missing;
2. new structural facts about a counterexample that hold for every n (§2, §3): the
   non-free member has minimal Jacobian degree r = d₁ − 1 in the narrow classes, all its Ziegler
   restrictions have exponents (d₁ − 1, d₂ + 1), and on every line its lowest logarithmic
   derivation is pinned to the multiple points (θ₁ ∧ E = c·∏ α_P); a congruence-free profile
   identity for the narrow classes; the unjoined-pair rule in its sharp form;
3. a precise account of where the known tools stop (§5): every tool bounds *numbers* (exponents,
   points per line, multiplicities, Tjurina numbers) and the numbers are consistent for all
   d₁ ≥ 6, so a proof has to use *positions*, which is what J is about.

Nothing here contradicts the possibility of a counterexample with n ≥ 18 lines. Sections 2–3 are
proved; §4 proves the reduction; §5 is discussion.

Notation as in `TC_ATTEMPT.md`: L a lattice of n lines; A₁ a free realisation with exponents
(1, d₁, d₂), d₁ ≤ d₂, δ = d₂ − d₁, n = d₁ + d₂ + 1; A₂ a non-free realisation of L; for a line H,
s_H = |A^H| its number of multiple points, μ_P = m_P − 1 the Ziegler multiplicities (Σ_{P∈H} μ_P
= n − 1), (e₁(H), e₂(H)) the exponents of the Ziegler multirestriction (A₂^H, μ); r = mdr(A₂) the
least degree of a logarithmic derivation of A₂ that is not a multiple of the Euler derivation E;
m the largest multiplicity; b₂ = Σ_P (m_P − 1) = d₁d₂ + d₁ + d₂; N the number of multiple points.

## 1. The window, for every n

Collecting `TC_ATTEMPT.md` Theorem 1 and §11 (all proved there or in the cited papers), a
counterexample (A₁, A₂) satisfies, for every n:

| | constraint | source |
|---|---|---|
| (W1) | d₁ ≥ 6 and δ ≤ d₁ − 4, so n = 2d₁ + δ + 1 ∈ [2d₁ + 1, 3d₁ − 3] | Abe 2014 Cor. 5.4, 5.5 |
| (W2) | δ + 4 ≤ s_H ≤ d₁ for every line H | Abe 2012, Abe 2014, Abe 2016 |
| (W3) | m ≤ d₁ − 1 | Dimca–Kühne–Pokora |
| (W4) | m ≤ r ≤ d₁ − 1, and n + 1 − 2r ≤ s_H ≤ r + 1, s_H − 1 ≤ e₁(H) ≤ r | Theorem 1(7) |
| (W5) | unjoined points P, Q (no line of A through both): (m_P − 1) + (m_Q − 1) ≤ d₁, or = d₂ with no third multiple point on PQ | §11 |
| (W6) | the moduli of L move the points of every line | Theorem 1(8) |
| (W7) | n ≥ 18 once the 17-line run finishes (n ≤ 16 proved, n ≤ 14 Barakat–Kühne) | §10 |

For fixed d₁ the admissible exponent pairs are (d₁, d₁ + δ) with 0 ≤ δ ≤ d₁ − 4. Two ends of
the range behave differently and are treated separately below:

- **narrow classes**, δ = d₁ − 4, n = 3d₁ − 3: (6, 8), (7, 10), (8, 12), (9, 14), …; every line
  has exactly s_H = d₁ points (W2 is an equality);
- **balanced classes**, δ ∈ {0, 1}: (7, 7), (7, 8), (8, 8), (8, 9), …; lines may have anywhere
  from 4 or 5 up to d₁ points.

All the free, non-divisionally-free lattices found inside the window so far (five at n = 15,
none at n = 16) are balanced, have ε := d₁ − m = 2, and have at least five points on every line.

## 2. The narrow classes: what the equalities force

Throughout this section δ = d₁ − 4, n = 3d₁ − 3, d₂ = 2d₁ − 4, and (A₁, A₂) is a counterexample.

**Proposition 2.1 (numerical structure).**
1. Every line has exactly d₁ multiple points, m ≤ d₁ − 1, and r = mdr(A₂) = d₁ − 1.
2. Every Ziegler restriction of A₂ has exponents (d₁ − 1, d₂ + 1); every Ziegler restriction of
   A₁ has exponents (d₁, d₂). The Yoshinaga defect d₁d₂ − e₁e₂ equals d₁ − 3 on every line.
3. Σ_P m_P = 3d₁(d₁ − 1), Σ_P m_P² = 12(d₁ − 1)², N = d₁² − 2d₁ + 4.
4. τ(A₂) = τ(A₁) = (n − 1)(n − 1 − r) + r² − (d₁ − 3): in du Plessis–Wall terms, A₂ misses
   freeness by exactly d₁ − 3 (nearly free would be 1).
5. There is no real, hence no rational, counterexample in a narrow class, for any n.

*Proof.* (1) (W2) with δ + 4 = d₁ gives s_H = d₁. (W4) gives r ≥ s_H − 1 = d₁ − 1 and r ≤ d₁ − 1.
(2) (W4): s_H − 1 ≤ e₁ ≤ r, so e₁ = d₁ − 1 and e₂ = n − 1 − e₁ = d₂ + 1. For A₁ this is Yoshinaga's
criterion. The defect is d₁d₂ − (d₁ − 1)(d₂ + 1) = δ + 1.
(3) Σ m_P = Σ_H s_H = n d₁; Σ m_P(m_P − 1) = n(n − 1) (each pair of lines meets once); N = Σ m_P − b₂
with b₂ = d₁d₂ + n − 1. Substitute n = 3d₁ − 3, d₂ = 2d₁ − 4.
(4) τ is lattice-determined (Σ (m_P − 1)²), equal to (n − 1)² − d₁d₂ on the free member; expand
with r = d₁ − 1.
(5) Melchior: a real arrangement that is not a pencil has 2N − b₂ ≥ 3. Here 2N − b₂ =
2(d₁² − 2d₁ + 4) − (2d₁² − d₁ − 4) = 12 − 3d₁ < 3 for d₁ ≥ 4. ∎

So in a narrow class the two members differ on *every* line by the same amount, the whole
discrepancy being carried by the minimal derivation of A₂, which exists one degree below the
free exponent. The next statement says where that derivation sits geometrically.

**Proposition 2.2 (the minimal derivation of A₂ is pinned to the multiple points).** Let θ be a
logarithmic derivation of A₂ of degree r = d₁ − 1 not in S·E, and ω := θ × E, i.e.
ω = (θ₂z − θ₃y, θ₃x − θ₁z, θ₁y − θ₂x), a vector of three forms of degree d₁ with x ω₁ + y ω₂ + z ω₃ = 0,
equivalently a section of Ω¹_{P²}(d₁ + 1). Then for every line H of A₂ with normal vector n_H:

  ω|_H = f_H · n_H  with  f_H = c_H · ∏_{P ∈ H} α_P,  c_H ∈ ℂ*.

That is: along every line of A₂ the field ω vanishes exactly at the multiple points of that line,
each to first order, and nowhere else. In particular Z(ω) ⊇ all N multiple points, while
deg c₂(Ω¹(d₁ + 1)) = d₁² − d₁ + 1 = N + (d₁ − 3).

*Proof.* Put θ_H := θ − (θ(α_H)/α_H) E; it is logarithmic, of degree d₁ − 1, and θ_H(α_H) = 0, so
it is tangent to H and restricts to an element of the Ziegler module D(A₂^H, μ). The restriction
is nonzero (the kernel of restriction on {θ(α_H) = 0} is α_H·{θ(α_H) = 0}, which lives in degrees
≥ d₁; `TC_ATTEMPT.md` Thm 1(7), step 2). By 2.1(2) the Ziegler module has exponents
(d₁ − 1, d₂ + 1) with d₁ − 1 < d₂ + 1, so its degree-(d₁ − 1) part is one-dimensional, spanned by
the lowest basis element θ₁^H, and θ_H|_H = c θ₁^H. Now the dichotomy of Thm 1(7), step 3: if
θ₁^H were a multiple of E_H then e₂ ≤ s_H − 1 = d₁ − 1, against e₂ = d₂ + 1; so det(θ₁^H, E_H) is a
nonzero binary form of degree e₁ + 1 = d₁ = s_H vanishing at the s_H points (each α_P divides it
because α_P^{μ_P} | θ₁^H(α_P) and μ_P ≥ 1), hence equals c'·∏ α_P. Finally θ × E = θ_H × E, and on H
both θ_H and E are tangent to H, so θ_H × E is normal to H with coefficient det(θ_H|_H, E_H). The
Chern class: c(Ω¹_{P²}) = 1 − 3h + 3h², so c₂(Ω¹(k)) = k² − 3k + 3; with k = d₁ + 1 this is
d₁² − d₁ + 1, and N = d₁² − 2d₁ + 4 by 2.1(3). ∎

**Proposition 2.3 (the same object for the free member).** Let θ₁ be a basis element of degree d₁
of D₀(A₁) and ω₁ := θ₁ × E, a section of Ω¹(d₁ + 2). Then for every line H,
ω₁|_H = c_H · ∏_{P ∈ H} α_P · ℓ_H · n_H with ℓ_H a linear form: on each line the field vanishes at the
d₁ multiple points and at exactly one further point (or doubly at a multiple point). Moreover, if
(E, θ₁, θ₂) is a Saito basis, ω₁ · θ₂ = c·Q where Q is the defining polynomial of A₁.

*Proof.* As in 2.2, with exponents (d₁, d₂): the lowest Ziegler element is not a multiple of E_H
(else e₂ ≤ d₁ − 1 < d₂), so det(θ₁^H, E_H) has degree d₁ + 1 and vanishes at d₁ points. The last
claim is Saito's criterion: det(E, θ₁, θ₂) = c·Q, and (θ₁ × E)·θ₂ = det(θ₁, E, θ₂). ∎

**What this does and does not give.** The field ω of 2.2 is a global, degree-d₁ object attached
to the non-free member, and its restriction to every line is determined by the lattice up to the
constants c_H. The d₁ − 3 "extra" zeros of ω (counted in the length of Z(ω) when it is finite) lie
on no line of A₂ or are embedded at multiple points in directions transverse to all lines there.
A contradiction would follow if one could show that a section of Ω¹(d₁ + 1) cannot vanish at the
N points of such a configuration with only d₁ − 3 residual length — a statement about the
*positions* of the points, which is where every numerical argument stops (§5). I could not prove
it. For d₁ = 6 (n = 15, class (6, 8)) it is true by Theorem 10.1, since no counterexample exists;
for d₁ = 7 (n = 18, class (7, 10)) it is exactly what a counterexample would have to violate.

## 3. The other classes: what survives

For δ < d₁ − 4 the window (W2) is an interval, and the equalities of §2 are replaced by
inequalities. What remains true for every n:

**Proposition 3.1.** In a counterexample with exponents (d₁, d₂):
1. r = mdr(A₂) ≥ max(m, s_max − 1, (n + 1 − s_min)/2), where s_max, s_min are the largest and
   smallest numbers of points on a line; and r ≤ d₁ − 1.
2. On every line, e₁(H) ≤ r ≤ d₁ − 1; the Yoshinaga defect d₁d₂ − e₁(H)e₂(H) ≥ δ + 1, with
   equality exactly on the lines where e₁(H) = d₁ − 1.
3. The minimal derivation θ of A₂ defines ω = θ × E, a section of Ω¹(r + 2), and on every line H
   ω|_H = f_H n_H with f_H a binary form of degree r + 1 divisible by ∏_{P ∈ H} α_P; so ω vanishes
   at every multiple point of L and, on the line H, at r + 1 − s_H further points (with
   multiplicity). Length count: c₂(Ω¹(r + 2)) = r² + r + 1 (for r = d₁ − 1 this is the
   d₁² − d₁ + 1 of 2.2).
4. Deleting any line H from the pair gives A₁ ∖ H plus-one generated with mdr = d₁ and A₂ ∖ H with
   mdr ≤ d₁ − 1: every one-line deletion of a counterexample lattice is a lattice with n − 1 lines
   whose realisations have different minimal Jacobian degrees.

*Proof.* (1), (2) are (W4) and Yoshinaga's inequality e₁e₂ ≤ d₁d₂; the defect is a decreasing
function of e₁ on [0, (n − 1)/2], and e₁ ≤ d₁ − 1 gives the bound. (3) is the argument of 2.2
with degree r in place of d₁ − 1: the restriction of θ_H to H is a nonzero element of degree r of
the Ziegler module, and not a multiple of E_H: if θ_H|_H = g·E_H then α_P^{μ_P} | g·α_P forces
α_P^{μ_P − 1} | g for all P ∈ H, so r − 1 = deg g ≥ Σ_P (μ_P − 1) = n − 1 − s_H, i.e.
s_H ≥ n − r ≥ d₂ + 2 > d₁ ≥ s_H, a contradiction. Hence det(θ_H|_H, E_H) is a nonzero form of
degree r + 1 vanishing at the s_H points. (4) is `TC_ATTEMPT.md` §11, item 4 (Abe's
deletion theorem for the free member; derivations restrict for the other). ∎

**Proposition 3.2 (what the enumerations add).** For n = 15 and 16 every lattice in the window
was enumerated. The free ones have: ε = d₁ − m = 2, s_min ≥ 5, and exist only in the classes
(7, 7) (five lattices) — none in (6, 8) or (7, 8). The non-realisable ones are the overwhelming
majority (99.9998% of the candidates for (7, 7); 100% for (7, 8)). I take two things from this:
- the window is far from sharp: the real obstruction is *realisability* of dense linear spaces,
  not the numerical conditions;
- the ε = 1 case (m = d₁ − 1), where Dimca–Kühne–Pokora show A₂ would be plus-one generated,
  produced no free lattice at all at n = 15, 16, although it is the densest part of the window.
  A theorem "ε = 1 ⇒ TC" would be the natural next piece of the window; I did not find a proof.

## 4. The reduction: one statement about positions

Fix a lattice L in the window and an irreducible component C of its non-degenerate realisation
space over ℂ (realisations with exactly the incidences of L, modulo PGL₃). Suppose C contains a
free point, with exponents (d₁, d₂). For a line H of L let

  π_H : C → M_H

be the map sending a realisation to the configuration of its s_H multiple points on H, with their
Ziegler multiplicities μ, modulo PGL₂ (so M_H has dimension s_H − 3), and let

  Σ_H := { configurations in M_H whose multiarrangement has e₁ ≤ d₁ − 1 }.

**Theorem 4.1 (the jump loci coincide).**
1. Σ_H is closed in M_H (e₁ is lower semicontinuous: {e₁ ≥ d} is open).
2. For every line H, the non-free locus of C is Z(C) = π_H^{-1}(Σ_H). It is a proper closed subset.
3. Consequently π_H^{-1}(Σ_H) = π_{H'}^{-1}(Σ_{H'}) for all lines H, H' of L.

*Proof.* (1) The conditions α_P^{μ_P} | θ(α_P) in a fixed degree d are linear in the coefficients
of θ, with a matrix depending polynomially on the positions; D(A, μ)_d ≠ 0 is a rank condition,
closed; e₁ ≤ d iff D_d ≠ 0. (2) Yoshinaga: A is free iff the Ziegler restriction to H has
exponents (d₁, d₂); for a non-free A the restriction has e₁e₂ < d₁d₂ with e₁ + e₂ = n − 1, i.e.
e₁ ≤ d₁ − 1; so A ∈ Z(C) iff π_H(A) ∈ Σ_H. Z(C) ≠ C because C contains a free point, and it is
closed (openness of freeness, or (1)). (3) follows. ∎

The content of TC for L is thus that the closed set Z(C), which is cut out on C by *any single
line's* jump condition, is empty. Equivalently:

**Conjecture J (jumping forces an incidence).** Let A be an arrangement of n lines whose lattice
lies in the window of §1 and admits a free realisation with exponents (d₁, d₂). If for one line H
the Ziegler multirestriction (A^H, μ) has e₁ ≤ d₁ − 1, then A has an incidence not in L: two
points of A^H coincide, or some other line passes through a multiple point it should miss.

**Theorem 4.2.** For a lattice L with a free realisation: J for one line of L ⟺ J for every line
of L ⟺ TC holds for L. In particular TC in P² (all n) is equivalent to Conjecture J.

*Proof.* "J for H" says π_H^{-1}(Σ_H) ∩ (non-degenerate realisations) = ∅, i.e. Z(C) = ∅ for every
component C, by 4.1(2). ∎

The equivalence is trivial once 4.1 is in place; its value is in what it allows: to prove TC for
L it is enough to prove J for **one** line, and one may choose the line with the fewest points.
Then M_H has dimension s_H − 3 ≤ d₁ − 3, the multiplicity vector μ on H is known, and Σ_H is the
degeneracy locus of an explicit (n − 1) × 2d₁ matrix of binomial-coefficient entries in the point
coordinates (the Yoshinaga matrix of `work/tc15/exact_tc_lattice.py`). The whole question is
whether the preimage of this hypersurface-like locus under π_H can meet C away from the
degenerations.

**What is proved about J.**
- (J1) s_H ≤ 3: M_H is a point, Σ_H is empty or everything, so J holds (rigidity principle). In
  the window s_H ≥ 4, so this never applies directly, but (W6) is its consequence: the moduli of
  L must move the points on *every* line, otherwise π_H(C) is a point, Σ_H ∌ π_H(C) because C has a
  free point, and Z(C) = ∅.
- (J2) Disjointness criterion: if for two lines H, H' one shows π_H^{-1}(Σ_H) ∩ π_{H'}^{-1}(Σ_{H'})
  ⊆ Deg, then J holds on C (by 4.1(3) the two preimages are equal, so each is inside Deg).
- (J3) For the 15-line (7, 7) family with moduli, J is proved exactly: the non-free locus
  d⁷(d − 1)⁷(2d − 1)³(d² − 3d + 1)(d² + d − 1) lies in the degeneration locus
  d²³(d − 1)²³(2d − 1)⁷(d² − 3d + 1)⁴(d² + d − 1)⁴ (`TC_ATTEMPT.md` Theorem 2). Note the shape: every
  factor of the jump polynomial is a factor of the degeneration polynomial, with a *smaller*
  exponent; the jump happens at the parameter values where extra triple points appear, and
  nowhere else.
- (J4) For every lattice of the window with n ≤ 16, and for every lattice of n = 17 processed so
  far, J holds (no mixed lattice), vacuously in most cases because the lattices are rigid or not
  realisable.

## 5. Why the known tools cannot finish, and what the patterns suggest

**Every available tool is numerical.** Abe's bounds (points per line, exponent gaps, heavy
points), the division theorem, du Plessis–Wall and the Dimca–Sticlaru classification by τ and
mdr, Dimca–Kühne–Pokora's multiplicity bounds: each constrains integers attached to L and to
(A₁, A₂). Section 1 is the complete list of what they give, and §2–3 show the integers are
mutually consistent for every d₁ ≥ 6: the narrow class with d₁ = 7 (n = 18) is numerically as
consistent as the one with d₁ = 6 that Theorem 10.1 excluded by enumeration. No further
inequality can come from these sources without a new idea, because the free member A₁ already
realises every numerical value, and A₂ only differs by having an extra derivation of degree
r ≤ d₁ − 1 — a statement about the positions of the points (Theorem 4.1).

**The patterns, and what they would mean.**
1. *Non-realisability dominates.* Nearly all lattices in the window are abstract linear spaces
   with no complex realisation. A proof along these lines would be a realisability theorem for
   dense linear spaces (every line with ≥ δ + 4 points, every pair of lines meeting at a point
   of multiplicity ≤ d₁ − 1, b₂ = d₁d₂ + n − 1). Such theorems exist only for very symmetric
   objects (nets, Latin squares — the (6, 8) class at n = 15 was decided by Kneser's theorem).
2. *Free lattices in the window have ε = 2 and are "almost rigid".* All five free (7, 7)
   lattices at n = 15 have realisation spaces of dimension ≤ 1; the one with moduli is a
   configuration of 15 lines through 6 points (a dual point configuration with an incidence
   theorem). Corollary C of `TC_ATTEMPT.md` says free non-DF lattices with moduli need excess
   dimension, i.e. an incidence theorem. A proof of TC for lattices with *zero-dimensional*
   realisation space is already available and trivial: a single Galois orbit cannot contain both
   free and non-free members (freeness is Galois-invariant), so a mixed lattice needs ≥ 2
   components over ℚ or a positive-dimensional one (Prop. A in `TC_ATTEMPT.md` §11). The hard
   case is therefore lattices with moduli, and J3 is the only such case fully computed; it
   satisfied J with the jump locus *inside* the degeneration locus with smaller exponents.
3. *The shape of J3 suggests a mechanism.* On the 15-line family the Yoshinaga matrix drops rank
   exactly where an extra triple point is created, and the order of vanishing of the rank-drop
   polynomial (7, 7, 3, 1, 1) is smaller than that of the degeneration polynomial (23, 23, 7, 4, 4).
   If this is general — "the jump polynomial divides a power of the degeneration polynomial" —
   it would say that the Yoshinaga matrix of a line sees only the incidences of the whole
   arrangement, which is a form of Conjecture J. I have no proof and no counter-model: the
   9-line Dimca–Pokora pair shows that *mdr* can change without an incidence change, but there
   neither member is free, so it does not touch J.

**An induction that does not close.** By 3.1(4), a counterexample with n lines yields, by deleting
any line, an (n − 1)-line lattice realised with two different values of mdr (d₁ for the
plus-one-generated deletion of A₁, ≤ d₁ − 1 for the deletion of A₂). If "mdr is combinatorial
for lattices in the range mdr < (n − 2)/2" (Dimca's Conjecture 1.2) held, TC would follow for
all n by induction from the verified base. Dimca–Pokora's 9-line pair refutes the stronger form
of that conjecture, so this route is closed unless one restricts it to lattices that are
deletions of window lattices — a class nobody has characterised.

## 6. Conclusion

- TC in P² for all n is equivalent to Conjecture J (Theorem 4.2): a jump of the Ziegler
  exponents on one line forces an extra incidence. This is not a new idea in spirit (it is what
  Yoshinaga's criterion already says), but Theorem 4.1 makes it precise as a statement about the
  single map π_H : C → M_H for a well-chosen line, which is where a proof would have to work.
- In the narrow classes (δ = d₁ − 4) the counterexample is rigidly structured (Propositions
  2.1–2.3): the non-free member carries a section ω of Ω¹(d₁ + 1) vanishing at all N multiple
  points with only d₁ − 3 residual length, and simply at each multiple point along each line. I
  could not show this is impossible; for d₁ = 6 it is (Theorem 10.1), for d₁ = 7 it is open and is
  the (7, 10) class of n = 18.
- The real content of the verified cases (n ≤ 16, soon 17) is not that the numerical window is
  empty but that the dense linear spaces it allows are almost never realisable, and the few free
  ones that are realisable are rigid or satisfy J by an explicit computation. Turning this into
  a theorem needs either a realisability theorem for dense linear spaces or a proof of J, and I
  do not see the mechanism for either. This is where I stop.
