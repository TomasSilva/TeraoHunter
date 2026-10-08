# Terao's conjecture in P²: the foliation of the minimal derivation

**Status (8 Oct 2026).** This is not a proof of Terao's conjecture (TC). It is one new tool and
what it proves. Attached to the non-free member A₂ of a counterexample is the foliation F of P²
defined by its logarithmic derivation of least degree r = mdr(A₂); every line of A₂ is invariant.
Its inflection curve J (the Jacobian of the Gauss map of F, Pereira's first extactic curve) has
degree 3r, is divisible by the product Q of the n lines, and vanishes to order ≥ 2 at every
singular point of F. Comparing the zeros of J, the singular points of F and the lattice gives an
exact identity (Theorem 4.1):

  (d₁ − r)(d₂ − r) = Σ_{P ∉ A₂} μ_P − Σ_{P ∈ A₂} d_P,

whose right-hand side counts singularities of F that are not accounted for by the lines through
them. When deg(J/Q) = 3r − n ≤ 1 the right-hand side is forced to be 0 while the left-hand side
is ≥ 1. Hence (Theorem 5.1):

> **A counterexample to TC in P² has 3·mdr(A₂) ≥ n + 2; hence d₂ ≤ 2d₁ − 6, d₁ ≥ 7 and n ≥ 15.**

This removes, for every n, the two outermost exponent classes of the window of `TC_GENERAL.md` §1:
δ = d₁ − 4 (the narrow classes of `TC_GENERAL.md` §2: (6, 8), (7, 10), (8, 12), …) and δ = d₁ − 5
((6, 7), (7, 9), (8, 11), …), and inside the remaining classes it removes every value of mdr(A₂)
below (n + 2)/3. Concretely:

- 14 lines: the single class (6, 7) is excluded by theory. Together with Dimca–Ibadula–Măcinic
  (13 lines) this is a theoretical proof of TC for n ≤ 14; the Barakat–Kühne computation is no
  longer needed for the bound.
- 15 lines: (6, 8) is excluded by theory (Theorem 10.1 of `TC_ATTEMPT.md` did it by enumeration
  plus Kneser's theorem); (7, 7) survives, with mdr(A₂) = 6 forced.
- 17 lines: (7, 9) and the part r = 6 of (8, 8) — exactly the two classes the running
  enumeration covers — are excluded by theory. The only open 17-line class is part [B] of
  (8, 8): A₂ nearly free with exponents (7, 10).
- 18 lines: (7, 10) excluded; (8, 9) survives with mdr(A₂) = 7.

What is *not* proved: TC. For deg(J/Q) ≥ 2 the identity allows a nonzero right-hand side carried
by singular points of the curve J/Q, and §6 says exactly what those singularities would have to
be. Sections 1–5 are proved; §6 is analysis of the remaining cases; §7 lists the computer checks
of the structure on known free arrangements.

Notation as in `TC_ATTEMPT.md` and `TC_GENERAL.md`: (A₁, A₂) a counterexample with lattice L of
n lines; A₁ free with exponents (1, d₁, d₂), d₁ ≤ d₂, δ = d₂ − d₁, n = d₁ + d₂ + 1; A₂ not free
with defining polynomial Q = ∏ α_H; S = ℂ[x, y, z], E = x∂ₓ + y∂_y + z∂_z; D₀(A₂) the derivations
with θ(α_H) ∈ (α_H) for all H and θ(Q) = 0 (the Jacobian syzygies), r = mdr(A₂) its least degree;
m_P the multiplicity of a point, s_H = |A^H| the number of multiple points on a line,
b₂ = Σ_P (m_P − 1) = d₁d₂ + n − 1.

## 0. Inputs

Two facts from `TC_ATTEMPT.md` Theorem 1 are used, nothing else from the window:

- (I1) r = mdr(A₂) ≤ d₁ − 1. (du Plessis–Wall, as in Dimca–Kühne–Pokora: τ = Σ_P (m_P − 1)² is a
  lattice invariant; among arrangements with the τ of a free lattice the free ones are exactly
  those with mdr = d₁, and a non-free one has mdr < d₁.)
- (I2) s_H ≤ d₂ + 1 for every line H. (Ziegler: the multirestriction of the free member A₁ to H
  has exponents (d₁, d₂); a 2-multiarrangement with s points and exponents (e₁, e₂) has
  s ≤ max(e₁, e₂) + 1, since a basis element not proportional to the Euler derivation E_H has a
  determinant with E_H of degree e_i + 1 vanishing at all s points.)

Also r ≥ 1, because a derivation of degree 0 in D₀(A₂) would be a constant vector field tangent
to every line, i.e. all lines through one point.

## 1. The foliation of θ and its singular points

Fix θ ∈ D₀(A₂) of degree r, θ ≠ 0, and put ω := θ × E, i.e.

  ω = (θ₂z − θ₃y, θ₃x − θ₁z, θ₁y − θ₂x),

three forms of degree r + 1 with xω₁ + yω₂ + zω₃ = 0: a section of Ω¹_{P²}(r + 2). It defines
the foliation F of P² whose tangent line at a point [x] is the line through [x] and [θ(x)]
(normal vector ω(x)); F has degree r. Singular points of F: the points with θ(x) ∈ ℂx, i.e.
ω(x) = 0. In an affine chart (say z = 1) F is given by the planar vector field
X = (θ₁ − uθ₃, θ₂ − vθ₃)(u, v, 1); at a singular point P we write L_P for its linear part,
μ_P := dim 𝒪_P/(X₁, X₂) for its Milnor number, and we call P

- *nondegenerate* if det L_P ≠ 0 (equivalently μ_P = 1), and then *radial* if L_P is a nonzero
  scalar, *non-radial* otherwise;
- *degenerate* if det L_P = 0; of *order k* if the Taylor expansion of X at P starts in degree k.

For a singular point P let m'_P be the number of lines of A₂ through P: m'_P = m_P ≥ 2 for a
multiple point of A₂, m'_P = 1 for a *free* point (on one line, not a multiple point), m'_P = 0
for an *off-line* point.

**1.1 (invariant lines, `TC_GENERAL.md` 3.1(3)).** Every line H of A₂ is F-invariant: θ(α_H) ∈ (α_H)
says that θ maps the plane H ⊂ ℂ³ into itself, so ω = θ × E is normal to H along H:

  ω|_H = f_H · n_H,   f_H := det_H(θ|_H, E|_H) ∈ Sym^{r+1}(H*),

and f_H ≠ 0 (proof below). The zeros of f_H are exactly the singular points of F on H (for P ∈ H,
f_H(P) = 0 iff θ(P) ∈ ℂP); every multiple point of A₂ on H is one (two invariant lines through P
force θ(P) ∈ H ∩ H' = ℂP). Write ν_P(H) := ord_P f_H; then Σ_{P ∈ H} ν_P(H) = r + 1.

*Proof that f_H ≠ 0.* Put K_H := θ(α_H)/α_H and θ_H := θ − K_H E; it is tangent to H and
restricts to the Ziegler multiarrangement (A₂^H, μ), μ_P = m_P − 1, and f_H = det_H(θ_H|_H, E_H).
If f_H = 0 then θ_H|_H = g E_H with g ∈ ℂ[H]_{r−1}. If g = 0, α_H divides θ_H and θ_H/α_H is a
logarithmic derivation of A₂ of degree r − 1 not proportional to E, against the minimality of r
(`TC_ATTEMPT.md` Thm 1(7), step 2). If g ≠ 0, the Ziegler conditions α_P^{μ_P} | g α_P force
α_P^{μ_P − 1} | g, so r − 1 ≥ Σ_{P ∈ H} (μ_P − 1) = n − 1 − s_H, i.e. s_H ≥ n − r ≥ d₂ + 2 by (I1),
against (I2). ∎

**1.2 Lemma (saturation).** ω₁, ω₂, ω₃ have no common factor. Hence F has isolated singular points
and

  Σ_{P ∈ Sing F} μ_P = c₂(Ω¹_{P²}(r + 2)) = r² + r + 1.

*Proof.* Suppose f | ω with deg f = k ≥ 1 and write ω = f ω'. The kernel of (a, b, c) ↦ xa + yb + zc
on triples of forms is generated by the Koszul relations, so ω' = θ' × E for a vector field θ' of
degree r − k, and then θ − fθ' is killed by × E, i.e. θ = fθ' + gE with g ∈ S_{r−1}. No α_H
divides f: otherwise ω would vanish identically on H, i.e. f_H = 0, against 1.1. Hence for every
H, θ(α_H) = fθ'(α_H) + gα_H ∈ (α_H) gives θ'(α_H) ∈ (α_H): θ' is logarithmic for A₂, of degree
r − k < r, and not a multiple of E (else ω' = 0). Subtracting its Euler component (θ'(Q) = hQ with
h ∈ S_{r−k−1}; θ' − (h/n)E ∈ D₀(A₂)) gives a nonzero element of D₀(A₂) of degree < r, a
contradiction. The count is c₂ of the rank-2 bundle Ω¹(r + 2): c(Ω¹_{P²}) = 1 − 3h + 3h², so
c₂(Ω¹(k)) = k² − 3k + 3, which is r² + r + 1 for k = r + 2; μ_P is the length of the zero scheme
of ω at P. ∎

**1.3 Local facts.** Let P be a singular point.
(a) Every line of A₂ through P is an invariant line of the planar field X at P, hence an
eigendirection of L_P. So m'_P ≥ 3 forces L_P to be scalar (radial, or L_P = 0), and m'_P = 2
forces L_P to be diagonal in the basis given by the two lines.
(b) For a line H ∋ P: ν_P(H) = 1 iff the eigenvalue ε_H(P) of L_P along H is nonzero.
(Differentiate f_H = det_H(θ|_H, E|_H) at P: with θ(P) = λP, d(f_H)_P(w) = det_H((Dθ_P − λ)w, P)
for w ∈ H, and (Dθ_P − λ)w ≡ ε_H(P)·w mod ℂP.)
(c) For a line H ∋ P with local equation v, invariance gives X₂ = v·q, and
μ_P = I_P(X₁, v) + I_P(X₁, q) = ν_P(H) + I_P(X₁, q) ≥ ν_P(H).
(d) If L_P = diag(0, e) with e ≠ 0 in a basis whose first vector is the direction of a line H of
A₂ (a *saddle-node along H*), then μ_P = ν_P(H) by (c), because q(P) = e ≠ 0.

## 2. The inflection curve of F

Let M := (∂ω_i/∂x_j) be the Jacobian matrix of ω: ℂ³ → ℂ³ and

  J := det M ∈ S_{3r}.

Projectively, ω is the Gauss map G: P² ⇢ P²*, [x] ↦ [ω(x)] = T_{[x]}F, given by the net of
degree-(r + 1) forms spanned by ω₁, ω₂, ω₃, and J is the Jacobian of this net.

**Proposition 2.1.**
1. In the chart z = 1, with X the planar field of §1: J(u, v, 1) = (r + 1)·det(X, DX·X). So J is
   the locus of inflection points of the leaves of F (where the curvature of the integral curves
   of X vanishes): the first extactic curve of F in the sense of Pereira.
2. J ≠ 0.
3. α_H | J for every line H of A₂; hence Q | J and n ≤ 3r.
4. J vanishes to order ≥ 2 at every singular point of F.

Write J = Q·J', with J' ∈ S_j, j := 3r − n ≥ 0.

*Proof.* (1) Put x = λ(u, v, 1); the change of coordinates (u, v, λ) ↦ x has Jacobian λ², and
∂ω/∂λ = (r + 1)λ^r ω̄, ∂ω/∂u = λ^{r+1}∂_u ω̄, ∂ω/∂v = λ^{r+1}∂_v ω̄ with ω̄ := ω(u, v, 1), so
J·λ² = (r + 1)λ^{3r+2} det(∂_u ω̄, ∂_v ω̄, ω̄), i.e. J(u, v, 1) = (r + 1) det(∂_u ω̄, ∂_v ω̄, ω̄). Now
ω̄ = (X₂, −X₁, vX₁ − uX₂), and adding v times the second row and u times the first row to the third
row leaves a third row (−X₂, X₁, 0), so the determinant is
X₁²∂_uX₂ − X₂²∂_vX₁ + X₁X₂(∂_vX₂ − ∂_uX₁) = det(X, DX·X), the numerator of the curvature of the
integral curves of X.
(2) If J ≡ 0 every integral curve of X has zero curvature, so every leaf of F is contained in a
line. A foliation whose leaves are lines is a pencil: its leaves form a curve Γ ⊂ P²*, and the
number of leaves through a generic point, which is 1, equals deg Γ. The pencil through [v] is
defined by the form E × v, whose coefficients are linear, so a saturated ω defining it has r = 0,
against r ≥ 1 and Lemma 1.2.
(3) On an invariant line the leaf is the line itself, of curvature zero, so J|_H = 0 by (1).
(Directly: ω maps the plane H into ℂn_H, so dω_x has a kernel on H for every x ∈ H.) The lines
are distinct, so their product divides J, and deg Q = n ≤ 3r = deg J.
(4) Euler's relation for the forms ω_i of degree r + 1 reads M(x)·x = (r + 1)ω(x). Multiply by the
adjugate: J(x)·x = (r + 1)·adj(M(x))·ω(x). At a singular point P expand in t = x − P: ω(P + t) =
M(P)t + O(|t|²) and adj(M(P + t)) = adj(M(P)) + O(|t|), so the right-hand side is
(r + 1)adj(M(P))M(P)t + O(|t|²) = (r + 1)J(P)t + O(|t|²) = O(|t|²), because J(P) = 0 (M(P) kills
P ≠ 0). Since P + t ≠ 0, J(P + t) = O(|t|²). ∎

The bound n ≤ 3·deg F for the number of invariant lines of a foliation, and the curve J, are due
to Pereira (extactic curves); the proof above is self-contained. Applied to A₂: a counterexample
has **mdr(A₂) ≥ n/3**, which is at least as strong as every lower bound on r in
`TC_GENERAL.md` (W4), with equality exactly in the narrow classes.

## 3. The multiplicity of J at the singular points

Everything in this section is local at a singular point P, computed from
J = (r + 1)·det(X, DX·X) in a chart centred at P. Write X = X_k + X_{k+1} + ⋯ with X_i
homogeneous of degree i in t = (u, v) and X_k ≠ 0 (k = 1 iff L_P ≠ 0), E₀ := u∂_u + v∂_v the
local radial field, and for a homogeneous field Y let ⟨t, Y⟩ := det(t, Y) = uY₂ − vY₁, the binary
form whose roots are the lines through P invariant under Y. A line through P invariant under X is
invariant under every X_i, so it is among the roots of the first ⟨t, X_i⟩ that is not identically
zero; such an i exists, because X = φ·E₀ with φ a unit would give J ≡ 0 near P.

**Proposition 3.1.** Let P be a singular point of F and ℓ_P the number of lines through P
invariant under F (so ℓ_P ≥ m'_P).
(a) If P is nondegenerate and non-radial: mult_P(J) = 2, and ℓ_P ≤ 2.
(b) If P is nondegenerate and radial, let i_P ≥ 2 be the least i with X_i not a multiple of E₀.
    Then mult_P(J) = i_P + 1 ≥ 3 and ℓ_P ≤ i_P + 1. In particular mult_P(J) ≥ max(3, ℓ_P).
(c) If L_P ≠ 0 is singular: mult_P(J) ≥ 3 and ℓ_P ≤ 2.
(d) If L_P = 0, of order k ≥ 2: if X_k is not a multiple of E₀ then ℓ_P ≤ k + 1 and
    mult_P(J) ≥ 3k − 1; if X_k = g·E₀ (dicritical leading term) and i_P > k is the least i with X_i
    not a multiple of E₀, then ℓ_P ≤ i_P + 1 and mult_P(J) = 2k + i_P − 1. In both cases
    mult_P(J) ≥ ℓ_P + 2.

*Proof.* Expand det(X, DX·X): DX_a·X_b has degree a + b − 1 and det(X_c, DX_a·X_b) degree
a + b + c − 1.
(a), (c): the degree-2 part is det(Lt, L²t). For L invertible with distinct eigenvalues λ₁, λ₂ and
t = αe₁ + βe₂ it is αβλ₁λ₂(λ₂ − λ₁)det(e₁, e₂) ≠ 0; for L = λ(I + N), N ≠ 0 nilpotent, λ ≠ 0, it is
λ²det(t, Nt) ≠ 0; for L scalar, for L of rank 1 (then L² ∝ L) and for L nilpotent (L² = 0) it
vanishes. A nonzero non-scalar L has at most two eigendirections, so ℓ_P ≤ 2.
(b): write X = φ·t + X_i + ⋯ with φ := c + g₁ + ⋯ + g_{i−2} the sum of the multiples of E₀
(c ≠ 0, i = i_P). Then DX·X = (φ·dφ(t) + φ² + dφ(X_i))·t + (i + 1)φX_i + (higher), using Euler
DX_i·t = iX_i, so det(X, DX·X) = ((i + 1)φ² − φ·dφ(t) − φ² − dφ(X_i))·⟨t, X_i⟩ + ⋯, whose lowest
term is i·c²·⟨t, X_i⟩, of degree i + 1, nonzero.
(d): if X_k is not a multiple of E₀ the lowest possible term is det(X_k, DX_k·X_k), of degree
3k − 1, and 3k − 1 ≥ (k + 1) + 2 ≥ ℓ_P + 2 for k ≥ 2. If X_k = gE₀ with deg g = k − 1 and i = i_P:
every X_a with k ≤ a < i is a multiple of E₀, and D(gE₀)·(hE₀) is a multiple of t, so
det(X_c, DX_a·X_b) vanishes unless one of a, b, c equals i; the lowest surviving degree is
2k + i − 1, from (c, a, b) ∈ {(i, k, k), (k, i, k), (k, k, i)}, with sum
(−k + i + 1)·g²·⟨t, X_i⟩ ≠ 0 (Euler: dg(t) = (k − 1)g, DX_i·t = iX_i). And
ℓ_P ≤ deg⟨t, X_i⟩ = i + 1 ≤ 2k + i − 3 for k ≥ 2. ∎

**Corollary 3.2 (multiplicity of J' = J/Q).** Since mult_P(Q) = m'_P ≤ ℓ_P:
- off-line points (m'_P = 0): mult_P(J') ≥ 2;
- free points (m'_P = 1): mult_P(J') ≥ 1, with equality only if P is nondegenerate non-radial;
- double points of A₂: mult_P(J') = 0 if P is nondegenerate non-radial, ≥ 1 if P is radial or
  L_P is singular and nonzero, ≥ 3 if L_P = 0;
- multiple points with m_P ≥ 3: P is radial or L_P = 0 (1.3(a)); mult_P(J') ≥ 0 in the first
  case and ≥ 2 in the second.

## 4. The excess identity

For a singular point P of F put

  d_P := Σ_{H ∈ A₂, H ∋ P} ν_P(H) − (m'_P − 1) − μ_P.

**Theorem 4.1.** (d₁ − r)(d₂ − r) = Σ_{P off-line} μ_P − Σ_{P on a line of A₂} d_P.

*Proof.* Count the zeros of all the f_H in two ways. By 1.1, Σ_H Σ_{P ∈ H} ν_P(H) = n(r + 1), and
the double sum equals Σ_P Σ_{H ∋ P} ν_P(H) = Σ_P (d_P + m'_P − 1 + μ_P) over all singular points.
Here Σ_P μ_P = r² + r + 1 (Lemma 1.2) and Σ_P (m'_P − 1) = b₂ − #{off-line points}, because the
multiple points of A₂ are exactly the singular points with m'_P ≥ 2 and contribute Σ(m_P − 1) = b₂,
free points contribute 0 and off-line points −1. For an off-line point d_P = 1 − μ_P. Hence

  n(r + 1) = Σ_{on} d_P − Σ_{off} μ_P + b₂ + r² + r + 1,

and with b₂ = d₁d₂ + n − 1 and d₁ + d₂ = n − 1 this is
d₁d₂ − r(n − 1 − r) = Σ_{off} μ_P − Σ_{on} d_P, whose left side is (d₁ − r)(d₂ − r). ∎

**Lemma 4.2 (tame points).** d_P = 0 if P is nondegenerate (whatever m'_P), and if P is a
saddle-node along a line H of A₂ (L_P = diag(0, e), e ≠ 0, zero eigendirection along H) with
m'_P ∈ {1, 2}.

*Proof.* Nondegenerate: μ_P = 1, and ν_P(H) = 1 for every H ∋ P by 1.3(b) (both eigenvalues are
nonzero), so d_P = m'_P − (m'_P − 1) − 1 = 0. Saddle-node along H: μ_P = ν_P(H) by 1.3(d); if
m'_P = 2 the other line H' is the e-eigendirection, so ν_P(H') = 1 and d_P = ν_P(H) + 1 − 1 − μ_P = 0;
if m'_P = 1, d_P = ν_P(H) − 0 − μ_P = 0. ∎

So for a counterexample the right-hand side of 4.1 is at least δ + 1 ≥ 1 (d₁ − r ≥ 1 by (I1) and
d₂ − r ≥ δ + 1), and by Lemma 4.2 it is carried entirely by off-line singular points, by free
points that are degenerate with their zero eigendirection transverse to the line, and by
degenerate points of A₂ with m'_P ≥ 3 or with L_P = 0. By Corollary 3.2 every one of these is a
point of multiplicity ≥ 2 of the curve J' = 0, of degree j = 3r − n. (Radial multiple points and
saddle-nodes along a line also lie on J' when m'_P ≤ 2, but they contribute nothing.)

For a free arrangement (d₁ = r) the identity says the right-hand side is 0; §7 checks this on
examples.

## 5. The theorem

**Theorem 5.1.** A counterexample (A₁, A₂) to Terao's conjecture in P² satisfies

  3·mdr(A₂) ≥ n + 2,  i.e.  deg(J/Q) ≥ 2.

Consequently (with (I1)) n ≤ 3d₁ − 5, i.e. d₂ ≤ 2d₁ − 6; and d₁ ≥ 7, n ≥ 15.

*Proof.* Let j = 3r − n. If j = 0, J' is a nonzero constant and mult_P(J') = 0 at every singular
point; by Corollary 3.2 there are no off-line and no free points, every double point of A₂ is
nondegenerate non-radial, and every point with m_P ≥ 3 is nondegenerate (radial). If j = 1, J' is a
linear form and mult_P(J') ≤ 1 everywhere: no off-line points; free points are nondegenerate
non-radial; double points are nondegenerate, or have L_P = diag(0, e) with e ≠ 0 in the basis of
the two lines (1.3(a)), i.e. are saddle-nodes along a line; points with m_P ≥ 3 are nondegenerate.
In both cases every singular point is tame (Lemma 4.2) and there are no off-line points, so the
right-hand side of Theorem 4.1 is 0, while (d₁ − r)(d₂ − r) ≥ δ + 1 ≥ 1. Hence j ≥ 2.
For the consequences: n + 2 ≤ 3r ≤ 3(d₁ − 1) gives n ≤ 3d₁ − 5 and δ = n − 2d₁ − 1 ≤ d₁ − 6; d₁ ≥ 6
by Abe (`TC_GENERAL.md` (W1)) and d₁ = 6 forces δ = 0, n = 13, excluded by
Dimca–Ibadula–Măcinic; so d₁ ≥ 7 and n = 2d₁ + δ + 1 ≥ 15. ∎

**Remark 5.2 (the extremal cases as a model).** The proof shows what j ∈ {0, 1} would force: a
saturated foliation of degree r with n = 3r (resp. 3r − 1) invariant lines whose singular points
are all nondegenerate multiple points of the arrangement (resp. plus nondegenerate free points on
the line J' = 0), and then Theorem 4.1 computes the lattice: (d₁, d₂) = (r, 2r − 1) with
N = r² + r + 1 (resp. (r, 2r − 2)). These are precisely the free arrangements with 3r or 3r − 1
lines whose minimal derivation is saturated: the full monomial arrangements
xyz(x^{r−1} − y^{r−1})(y^{r−1} − z^{r−1})(z^{r−1} − x^{r−1}) for j = 0, and e.g. their one-line
deletions for j = 1 (§7). So j = 0 and j = 1 do occur, but only when d₁ = r, i.e. for the free
member itself; the inequality of Theorem 5.1 is sharp in that sense.

**Corollary 5.3 (the window after Theorem 5.1).** A counterexample has d₁ ≥ 7, 0 ≤ δ ≤ d₁ − 6,
n ∈ [2d₁ + 1, 3d₁ − 5], and r = mdr(A₂) ∈ [⌈(n + 2)/3⌉, d₁ − 1]. For the classes with n ≤ 20:

| class | n | admissible r | j = 3r − n | (d₁ − r)(d₂ − r) |
|---|---|---|---|---|
| (6, 6) | 13 | – | – | excluded (Dimca–Ibadula–Măcinic) |
| (6, 7) | 14 | none | – | **excluded (5.1)**; was Barakat–Kühne (computer) |
| (6, 8) | 15 | none | – | **excluded (5.1)**; was Thm 10.1 (computer + Kneser) |
| (7, 7) | 15 | 6 | 3 | 1 |
| (7, 8) | 16 | 6 | 2 | 2 |
| (7, 9) | 17 | none | – | **excluded (5.1)** |
| (8, 8) | 17 | 7 | 4 | 1 |
| (7, 10) | 18 | none | – | **excluded (5.1)** |
| (8, 9) | 18 | 7 | 3 | 2 |
| (8, 10) | 19 | 7 | 2 | 3 |
| (9, 9) | 19 | 7, 8 | 2, 5 | 4, 1 |
| (8, 11) | 20 | none | – | **excluded (5.1)** |
| (9, 10) | 20 | 8 | 4 | 2 |

For (7, 7) and (7, 8) TC is known by enumeration (`TC_ATTEMPT.md` Theorems 10.2, 10.4), so the
first open case overall remains 17 lines, class (8, 8) with r = 7 — which is part [B] of
`TC_ATTEMPT.md` §11 (A₂ nearly free with exponents (7, 10)), the part *not* covered by the running
enumeration. The running enumeration ((7, 9), and (8, 8) with r = 6) has become an independent
computational check of Theorem 5.1 rather than the decisive step for n = 17.

## 6. What survives: j ≥ 2

Theorem 4.1 and Corollary 3.2 turn TC (for the lattices in the window) into the statement

  (★) the singular points of the curve J' = 0, of degree j = 3r − n ≤ d₁ − δ − 4, cannot carry
      the excess (d₁ − r)(d₂ − r) ≥ δ + 1.

This is where the argument stops. What can be said about the first open values of j:

**j = 2 (classes (7, 8), (8, 10), and (9, 9) with r = 7).** J' is a conic. A smooth conic has no
singular point, so the right-hand side of 4.1 would be 0: impossible. If J' is a pair of distinct
lines, the excess is carried by their intersection point P₀ alone, with mult_{P₀}(J) = m'_{P₀} + 2
exactly; Proposition 3.1 leaves three possibilities: P₀ off-line and nondegenerate non-radial
(then μ_{P₀} = 1, the excess is 1, so δ = 0, d₁ = d₂ = r + 1, n = 2r + 3 = 3r − 2, r = 5: the
13-line class, already excluded); P₀ a free point which is a saddle-node of order δ + 2 with its
zero eigendirection transverse to its line and mult_{P₀}(J) = 3; or P₀ a degenerate point of A₂
with m_{P₀} ≥ 3 of a specific kind (non-dicritical of order 2 with m_{P₀} = 3, or dicritical of
order 2 with m_{P₀} = i_{P₀} + 1) whose blown-up singularities carry the excess. If J' is a double
line ℓ², every point of ℓ has mult_P(J') = 2 and the excess can be spread along ℓ: if ℓ is not
F-invariant, ω|_ℓ is a section of 𝒪_ℓ(r), so ℓ carries at most r singular points of F; if ℓ is
F-invariant, A₂ ∪ ℓ is a saturated arrangement of 3r − 1 invariant lines to which the j = 1
analysis applies, and the excess equals the number of singular points on ℓ not on A₂. I have not
been able to exclude the last two situations.

**j ≥ 3.** J' is a cubic or worse and (★) asks which singular curves of degree j can occur as J/Q.
The data of §3–4 give at each excess point a precise local model (a node of J' at an off-line
nondegenerate point contributes exactly 1; a saddle-node free point of order p with transverse
centre contributes p − 1 and has mult_P(J) ≥ 3; a degenerate point of A₂ with m_P ≥ 3 contributes
−d_P with mult_P(J) ≥ m_P + 2), so for the smallest open cases the excess is 1 or 2 and must sit
at one or two prescribed singularities of J'. For (7, 7) at 15 lines J' is a cubic with exactly one
unit of excess; TC holds there by Theorem 10.2, so some further global constraint kills this
configuration, but I do not know which.

**Further structure, not used above.** The cofactors K_H := θ(α_H)/α_H ∈ S_{r−1} satisfy
Σ_H K_H = 0, and at a singular point P with θ(P) = λ_P P: K_H(P) = λ_P for H ∌ P, while
K_H(P) − λ_P is the eigenvalue of L_P transverse to H for H ∋ P. Summing over H:
Σ_{H ∋ P} (transverse eigenvalue of H at P) = −n λ_P at every singular point; so θ(P) = 0 at
every off-line point, and at a radial point with m_P lines the common eigenvalue is −nλ_P/m_P. On
a line H the Lagrange relation for K_H|_H at the r + 1 zeros of f_H, combined with the
Camacho–Sad formula Σ_{P ∈ H} CS(F, H, P) = 1, reads Σ_{P ∈ H} λ_P/ε_H(P) = −1, where
ε_H(P) = −c_H ∏_{P' ≠ P} det_H(P', P) is the eigenvalue along H (f_H = c_H ∏_{P'} det_H(P', ·)).
At a radial point all ε_H(P), H ∋ P, coincide; around a cycle of lines joined at radial points
this is a Carnot-type multiplicative condition on the positions of the other singular points on
those lines. These identities constrain positions; I did not get a contradiction from them.

## 7. Computer checks of the structure

`work/foliation/inflect.m2` (Macaulay2) computes, for an arrangement over ℚ: the least-degree
Jacobian syzygy θ, ω = θ × E and its gcd, J = det(∂ω/∂x), the factorisation of J/Q, the degree of
the zero scheme of ω, the singular points off the lines and those that are not multiple points,
the profile, b₂ and Δ := b₂ − n + 1 − r(n − 1 − r), which equals (d₁ − r)(d₂ − r) when the lattice
has a free member. Results (r = syzygy degree minus (n − 1)):

| arrangement | n | r | j | J/Q | deg Z(ω) | off-line | non-multiple | Δ |
|---|---|---|---|---|---|---|---|---|
| xyz(x²−y²)(y²−z²)(z²−x²), exp. (3, 5) | 9 | 3 | 0 | const | 13 = r²+r+1 | 0 | 0 | 0 |
| xyz(x⁴−y⁴)(y⁴−z⁴)(z⁴−x⁴), exp. (5, 9) | 15 | 5 | 0 | const | 31 | 0 | 0 | 0 |
| xyz(x−y)(x−z), exp. (2, 2) | 5 | 2 | 1 | a line | 7 | 0 | 1 | 0 |
| xy(x²−y²)(z²−x²)(z²−y²), exp. (3, 4) | 8 | 3 | 1 | z | 13 | 0 | 2 | 0 |
| B₃ minus x = 0, exp. (3, 4) | 8 | 3 | 1 | x (the deleted line) | 13 | 0 | 2 | 0 |
| xyz(x−y)(x−z)(y−z)(x+y−z)(x−y+z)(x−y−z) | 9 | 3 | 0 | const | 13 | 0 | 0 | 0 |

The last arrangement was entered as a non-free guess; Theorem 4.1 with j = 0 predicts that any
saturated extremal arrangement is free with exponents (r, 2r − 1) and N = r² + r + 1, and indeed
its profile came out as 3 quadruple, 4 triple, 6 double points, b₂ = 23, χ₀(0) = 15 = 3·5.
Arrangements with a line carrying ≥ n − r points (generic arrangements, arrangements with a point
of multiplicity > n/2, one-line perturbations of the examples above) have an unsaturated minimal
derivation (gcd(ω) ≠ 1) and are outside the scope of §1–5, as expected: saturation was proved from
(I1) and (I2), which hold in a counterexample but not in general.

## 8. Conclusion

- New theorem: a counterexample to TC in P² has 3·mdr(A₂) ≥ n + 2, hence d₂ ≤ 2d₁ − 6, d₁ ≥ 7,
  n ≥ 15. The narrow classes of `TC_GENERAL.md` §2 and the next classes δ = d₁ − 5 are gone for
  every n; the 14-line case has a theoretical proof; at 17 lines only the nearly-free part [B] of
  (8, 8) is open.
- The tool is geometric, not numerical: the inflection curve J of the foliation of the minimal
  derivation, J = Q·J', and the identity (d₁ − r)(d₂ − r) = excess of J'. It is the first
  constraint in this project that mixes the lattice numbers with the *positions* (through the
  singular points of J'), which is what `TC_GENERAL.md` §5 said a proof needs.
- TC itself is not proved. Within this framework it is equivalent, for the window lattices, to
  (★): the singular points of a curve of degree 3r − n ≤ d₁ − δ − 4 cannot carry excess
  (d₁ − r)(d₂ − r). The first open instances are j = 2 (a conic J' that is a line pair or a double
  line) and j = 3.

## Sources

- J. V. Pereira, Vector fields, invariant varieties and linear systems, Ann. Inst. Fourier 51
  (2001) 1385–1405 (extactic curves; the bound 3·deg F on invariant lines).
- C. Camacho, P. Sad, Invariant varieties through singularities of holomorphic vector fields,
  Ann. of Math. 115 (1982) (the index formula quoted in §6).
- M. Brunella, Birational geometry of foliations, IMPA (Milnor numbers under blow-up, used in the
  j = 2 discussion of §6).
- The sources of `TC_ATTEMPT.md` for (I1), (I2) and the window.
