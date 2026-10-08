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
exhaustively (C(57, 17) ≈ 10¹⁴); a random sample of 17-, 18- and 19-line subsets is being
analysed (`pg27sample.c`).
