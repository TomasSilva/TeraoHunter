# Over a finite field, free does not imply divisionally free

**Answer.** The claim "every free arrangement over F_q is divisionally free (DF)" is **false for
every finite field**. The smallest ranks of counterexamples are:

| field | smallest rank of a free, non-DF arrangement | example (§4) |
|-------|---------------------------------------------|--------------|
| F_q, q ≥ 4 | 3 | A_q: 3(q − 1) lines, exponents (1, q, 2q − 4) |
| F_3 | 4 | 22 planes, exponents (1, 3, 9, 9) |
| F_2 | 6 | 50 hyperplanes, exponents (1, 2, 4, 8, 16, 19) |

Below these ranks the implication holds:
- in rank ≤ 2 by definition;
- for F_3 in rank 3, and for F_2 in ranks 3, 4, 5, by exhaustive computation (§6);
- for F_2 in rank 3, also by a hand proof (§6).

Every higher rank also has counterexamples, obtained by adding coordinate hyperplanes (§5).

All the counterexamples come from one mechanism: **blocking sets**. The deletion theorem makes
"covering" arrangements free (§2). If the covering is minimal, no restriction can have a
characteristic polynomial that divides χ(A) (§3).

## 1. Notation

* V = F_q^ℓ and S = F_q[x_1, …, x_ℓ]. B is the set of all hyperplanes of V, so |B| = [ℓ]_q, where
  [m]_q = 1 + q + … + q^{m−1}.
* For A ⊆ B, A* ⊆ PG(V*) is the set of normals (points α_H). A^H is the restriction of A to H.
* **Divisionally free** (Abe, Invent. Math. 204 (2016), Thm–Def. 4.3):
  - arrangements of rank ≤ 2 are DF;
  - for ℓ ≥ 3, A is DF iff there is H ∈ A with χ(A^H; t) | χ(A; t) and A^H DF.
  In particular, **if no H ∈ A has χ(A^H) | χ(A), then A is not DF.**
* A set A* ⊆ PG(V*) is a **blocking set (with respect to lines)** if it meets every line. A line
  meeting A* in exactly one point h is **tangent at h**. The blocking set is **minimal** if every
  point of A* lies on a tangent line.
* **Dictionary.** A codimension-2 subspace X ⊂ V corresponds to the line X^⊥ ⊂ PG(V*), and
  X ⊂ H ⟺ α_H ∈ X^⊥. Hence:
  - A* is blocking ⟺ every codimension-2 subspace of V lies in some H ∈ A.
  - A* has a tangent at α_H ⟺ some hyperplane X of H lies in no other member of A ⟺ A^H is not
    the full arrangement of all [ℓ−1]_q hyperplanes of H.

## 2. The covering lemma

**Lemma 1.** Let A ⊆ B, and suppose A* is a blocking set. Then A is free with exponents

  (1, q, q², …, q^{ℓ−2}, d), d = |A| − [ℓ−1]_q.

In particular every A′ with A ⊆ A′ ⊆ B is free too, since A′* is also blocking.

*Attribution.* In rank 3, "A* is blocking" says that every point of PG(2, q) lies on a line of A,
i.e. χ(A, q) = 0. In that case Lemma 1 is Yoshinaga's theorem (Free arrangements over finite
fields, Proc. Japan Acad. 82 (2006)); Abe gives another proof (J. Singularities 8 (2014),
Cor. 6.2), where the Frobenius argument also appears (Prop. 6.1). I have not found the rank ≥ 4
version, which covers codimension-2 subspaces rather than points, in the literature.

*Proof.*
1. **B is free.** The derivations θ_j = Σ_i x_i^{q^j} ∂_i (0 ≤ j ≤ ℓ − 1) are logarithmic for
   every F_q-hyperplane, because θ_j(α) = α^{q^j}. Their determinant is the Moore determinant,
   which is a nonzero constant times ∏_{H∈B} α_H. So by Saito's criterion, which holds over any
   field, B is free with exponents (1, q, …, q^{ℓ−1}).
2. **Delete the missing hyperplanes.** List B \ A = {H_1, …, H_m} in any order, and set
   A_k = B \ {H_1, …, H_k}. Assume A_{k−1} is free with exponents (1, q, …, q^{ℓ−2}, q^{ℓ−1} − k + 1).
3. **The restriction is full.** A hyperplane X of H_k is a codimension-2 subspace of V, so it lies
   in some H ∈ A ⊆ A_{k−1}, and H ≠ H_k because H_k ∉ A. Hence A_{k−1}^{H_k} is the arrangement of
   all hyperplanes of H_k ≅ F_q^{ℓ−1}. By step 1 it is free with exponents (1, q, …, q^{ℓ−2}).
4. **Apply the deletion theorem.** Terao's deletion theorem (any field) now says that A_k is free
   with exponents (1, q, …, q^{ℓ−2}, q^{ℓ−1} − k).
5. **Finish.** Take k = m = |B| − |A|. Since |B| − q^{ℓ−1} = [ℓ−1]_q, the last exponent is d. ∎

*Remark (an algebraic certificate that avoids the deletion theorem).* Let C = B \ A and
Q_C = ∏_{H∈C} α_H. Suppose θ_{ℓ−1} ≡ Σ_{j<ℓ−1} f_j θ_j modulo Q_C·Der(S) for some f_j ∈ S. Put
θ′ = (θ_{ℓ−1} − Σ f_j θ_j)/Q_C.
* θ′ is logarithmic along A, because α_H ∤ Q_C for H ∈ A.
* det(θ_0, …, θ_{ℓ−2}, θ′) = det(θ_0, …, θ_{ℓ−1})/Q_C = Q_A.

So A is free with the exponents of Lemma 1, by Saito's criterion. This is a module-membership test
over S/(Q_C), and §7 uses it as an independent computer check.

*Remark.* Being blocking is visible in χ. Write b₂ = Σ_{X ∈ L₂(A)} (|A_X| − 1) for the
t^{ℓ−2}-coefficient of χ. The number of codimension-2 subspaces covered by A is
|A|·[ℓ−1]_q − b₂(A). So A* is blocking ⟺ |A|·[ℓ−1]_q − b₂(A) equals the Gaussian binomial
[ℓ choose 2]_q.

## 3. The theorem

**Theorem 2.** Let A* ⊆ PG(V*) be a minimal blocking set, and suppose d = |A| − [ℓ−1]_q ≥ q^{ℓ−2}.
Then A is free with exponents (1, q, …, q^{ℓ−2}, d), and A is not divisionally free.

*Proof.* Freeness is Lemma 1. Now suppose χ(A^H) | χ(A) for some H ∈ A.
1. Since χ(A) = (t − 1)(t − q)⋯(t − q^{ℓ−2})(t − d), we get χ(A^H) = χ(A)/(t − e) for some root e.
2. For a central arrangement, the sum of the roots of χ equals the number of hyperplanes. So
   |A^H| = |A| − e ≥ |A| − max(q^{ℓ−2}, d) = |A| − d = [ℓ−1]_q.
3. But α_H lies on a tangent line, so some hyperplane of H is missing from A^H. Therefore
   |A^H| ≤ [ℓ−1]_q − 1.

This is a contradiction, so no H ∈ A has χ(A^H) | χ(A), and A ∉ DF. ∎

## 4. Examples

### 4.1 q ≥ 4, rank 3: the arrangement A_q

Let A_q = {x_i = a x_j : 1 ≤ i < j ≤ 3, a ∈ F_q^*}, with defining polynomial

  Q_q = (x^{q−1} − y^{q−1})(y^{q−1} − z^{q−1})(z^{q−1} − x^{q−1}).

**The normals form a minimal blocking set.** The normals e_i − a e_j fill the three sides of the
coordinate triangle in PG(2, q)*, minus the three vertices.
* *Blocking.* A side contains q − 1 normals. A line through exactly one vertex meets the opposite
  side at a non-vertex point. A line through no vertex meets every side at a non-vertex point.
* *Minimal.* For h on side e₁e₂, the line h e₃ meets the three sides only in h and in the vertex
  e₃. So it is tangent at h.

**Conclusion.** Here d = (3q − 3) − (q + 1) = 2q − 4, and d ≥ q ⟺ q ≥ 4. By Theorem 2, A_q is
free with exponents (1, q, 2q − 4) and not DF for every q ≥ 4. (For q = 3, A_3 is the braid
arrangement, which is DF.) For q = 4, A_4 is the dual Hesse arrangement.

*Direct proof of freeness (independent of Lemma 1).*
* The three derivations below lie in D(A_q):
  - θ₀ = E (Euler);
  - θ₁ = Σ x_i^q ∂_i, since θ₁(α) = α^q;
  - θ₂ = (yz)^{q−2} ∂_x + (xz)^{q−2} ∂_y + (xy)^{q−2} ∂_z. On x = a y we have
    θ₂(x − a y) = z^{q−2} y^{q−2}(1 − a^{q−1}) = 0, and the other lines follow by symmetry.
* Put u = x^{q−1}, v = y^{q−1}, w = z^{q−1}. Then
  det(θ_k(x_i)) = uv² − uw² − u²v + u²w + vw² − v²w = (u − v)(v − w)(w − u) = Q_q.
* The degrees add up: 1 + q + (2q − 4) = |A_q|. By Saito's criterion these three derivations form
  a basis.

**Not DF, by counting.** Every line of A_q contains exactly q points: one coordinate point, of
multiplicity q − 1, and q − 1 triple points (ab : b : 1). So |A_q^H| − 1 = q − 1 ∉ {q, 2q − 4}.

*Second family.* Let Q = r². Over F_Q, consider the lines ax + by + cz = 0 with
a^{r+1} + b^{r+1} + c^{r+1} = 0. Their normals form the Hermitian unital, which is a minimal
blocking set of size r³ + 1. So by Theorem 2 these lines are free with exponents (1, r², r³ − r²)
and not DF. Checked for F₄, F₉ and F₁₆; the F₄ case is again the dual Hesse arrangement.

### 4.2 q = 3, rank 4: an elliptic quadric

Work in the dual space PG(3, 3) with coordinates y. Set:
* Q(y) = y₁² + y₂² + y₃y₄. Since y₁² + y₂² is anisotropic over F₃, Q is an elliptic quadric with
  10 points.
* p = (0:0:1:0) ∈ Q, with tangent plane T_p = {y₄ = 0}. Then T_p ∩ Q = {p}.
* L = {y₂ = y₄ = 0}, a tangent line at p.

Let **A* = L ∪ {y : y₄ ≠ 0, Q(y) ≠ 0}**, i.e. the 22 planes

  x₁, x₃, x₁ + x₃, x₁ − x₃, and a x₁ + b x₂ + c x₃ + x₄ for a, b, c ∈ F₃ with c ≠ −(a² + b²).

*Blocking.* Suppose a line M misses A*.
* If M ⊄ T_p, its three points off T_p are not in A*, so they lie on Q. But a line meeting a
  quadric in 3 points lies on it, and Q contains no lines.
* If M ⊂ T_p, it meets L ⊆ A*.

*Minimal.* Note that A* ∩ T_p = L and A* ∩ Q = {p}.
* *If h ∈ L:* every line through h in T_p, other than L, is tangent at h.
* *If h ∉ T_p ∪ Q:*
  1. The polar plane of h meets Q in a conic of 4 points. The remaining 6 points of Q lie on
     exactly 3 secants through h.
  2. A secant through h that meets L lies in the plane π = ⟨h, L⟩ ≠ T_p, and π ∩ Q is a conic.
     Since q is odd, h lies on 0 or 2 tangents of that conic, so at most 2 lines through h in π
     are secants.
  3. So some secant M through h misses L. Its other three points are two points of Q ∖ {p} and
     the point M ∩ T_p. That last point is not on L, and it is not on Q, because T_p ∩ Q = {p}
     and p ∉ M. None of these three points is in A*, so M is tangent at h.

*Conclusion.* d = 22 − 13 = 9 ≥ 9. By Theorem 2, A is free with exponents (1, 3, 9, 9), so
χ = (t − 1)(t − 3)(t − 9)², and A is not DF. Its restrictions have 10, 11 or 12 lines (never 13).
Its stabiliser in PGL(4, 3) has order 36.

This construction does not generalise to q ≥ 4. The set is still blocking, but it is not
minimal: a tangent at a point h off T_p ∪ Q would need its q − 1 ≥ 3 points off T_p (other than
h) to lie on Q. Checked for q = 5, 7.

### 4.3 q = 2, rank 6: a complete cap

Over F₂ a line has 3 points. So A* is blocking ⟺ its complement K is a **cap** (no three
collinear points). And A* is minimal ⟺ K is a **complete** cap (every point off K is a + b with
a, b ∈ K).

Take the 13-point complete cap of PG(5, 2)

  K = {e₁, …, e₆, e₃+e₄+e₅, e₂+e₄+e₆, e₁+e₄+e₅, e₁+e₃+e₅, e₁+e₃+e₄+e₆, e₁+e₂+e₃+e₆, e₁+e₂+e₃+e₄}

and let A be the 50 hyperplanes ker α with α ∉ K ∪ {0}. Equivalently, Q_A is the product of all
63 nonzero linear forms divided by the 13 forms in K.

Here d = 50 − 31 = 19 ≥ 16. By Theorem 2, A is free with exponents (1, 2, 4, 8, 16, 19) and not
DF.

More generally, every complete cap K ⊂ PG(ℓ − 1, 2) with |K| ≤ 2^{ℓ−2} works.

## 5. All higher ranks

**Lemma 3.** For any A, A × Φ₁ (add one coordinate hyperplane in a new variable) is free iff A is
free, and it is DF iff A is DF.

*Proof.*
* *Freeness.* A product of arrangements is free iff both factors are, and its exponents are the
  union of theirs. So A × Φ₁ is free iff A is, with exponents exp(A) ∪ {1}.
* *DF ⇐.* χ(A × Φ₁) = χ(A)(t − 1), and restricting to the new hyperplane gives A. So if A is DF,
  so is A × Φ₁.
* *DF ⇒.* Suppose A × Φ₁ is DF, witnessed by H. If rank A ≤ 2, A is DF anyway.
  - If H is the new hyperplane, the witness condition says that (A × Φ₁)^H = A is DF.
  - If H ∈ A, then (A × Φ₁)^H = A^H × Φ₁ is DF, and χ(A^H)(t − 1) | χ(A)(t − 1). By induction on
    rank, A^H is DF, and χ(A^H) | χ(A), so A is DF via H. ∎

So A_q × Φ_k (q ≥ 4), the §4.2 example × Φ_k (q = 3) and the §4.3 example × Φ_k (q = 2) are free
and not DF in every rank above the thresholds. (G(q−1, q−1, 4) over F₅ and F₇ also gives
irreducible rank-4 examples, by computation.)

## 6. Below the thresholds

**Exhaustive computation.** The search covered all hyperplane subsets of PG(ℓ − 1, q) up to
PGL(ℓ, q), using orderly generation with nauty.
* Over F₂ and F₃ the lattice determines the arrangement up to PGL, because binary and ternary
  matroids are uniquely representable. So the counts below are counts of lattices.
* DF was decided by Abe's recursive definition, and freeness by an exact computation over F_q.
* For F₃ in rank 4, exact local freeness was tested first. All 7 locally free candidates turned
  out to be free.

| field | rank | essential | split χ | DF | non-DF, split χ | non-DF and free |
|-------|------|----------:|--------:|---:|----------------:|----------------:|
| F₂ | 3 | 6 | 5 | 5 | 0 | 0 |
| F₃ | 3 | 25 | 19 | 18 | 1 (AG(2,3), not free) | 0 |
| F₂ | 4 | 36 | 19 | 19 | 0 | 0 |
| F₂ | 5 | 1,326 | 164 | 151 | 13 (none free; Macaulay2 agrees) | 0 |
| F₃ | 4 | 111,802 | 6,267 | 3,390 | 2,877 | **7** |

**The 7 free non-DF classes in PG(3, 3).**

| planes | exponents | classes |
|-------:|-----------|--------:|
| 22 | (1, 3, 9, 9) | 2, one of which is §4.2 |
| 23 | (1, 3, 9, 10) | 2 |
| 25 | (1, 3, 9, 12) | 3 |

For each of them A* is a minimal blocking set with d ≥ 9, so Theorem 2 accounts for every free
non-DF arrangement in PG(3, 3).

In fact the match is exact (`work/covering_minimal_pg33.py`).
* Up to PGL, the minimal blocking sets of PG(3, 3) have sizes 13 (a plane), 19, 20 (×2),
  21 (×9), 22 (×2), 23 (×2) and 25 (×3).
* The seven of size ≥ 22 — exactly the threshold d ≥ 9 — are the seven free non-DF classes.
* All 3,098 covering classes of PG(3, 3) have χ = (t − 1)(t − 3)(t − 9)(t − (|A| − 13)), as
  Lemma 1 predicts (`work/covering_chi_exhaustive.py`).

**Why the construction gives nothing smaller.**
* In rank 3, Theorem 2 needs a minimal blocking set of PG(2, q) with at least 2q + 1 points.
  Bruen–Thas give the bound q√q + 1 for minimal blocking sets, which is below 2q + 1 for q = 2, 3
  and equal to it for q = 4 (the unital).
* Over F₂ in ranks 4 and 5, Theorem 2 would need a complete cap of size at most 4, respectively 8.
  Random search found complete caps only of sizes 5, 8 in PG(3, 2) and 9, 10, 16 in PG(4, 2).

These are heuristics; the exhaustive table is the proof.

**Rank 3, q = 2, by hand.** Every essential free arrangement in F₂³ is DF.

*Ingredients* (all valid over any field):
* *Frobenius bound.* θ_F = Σ x_i^q ∂_i ∈ D(A) and θ_F ∉ S·E. So if A is free with exponents
  (1, d₁, d₂), then d₁ ≤ q. (Compare Abe 2014, Prop. 6.1.)
* *Ziegler.* If A is free with exponents (1, d₁, d₂), then the Ziegler multirestriction
  (A^H, m^H) is free with exponents (d₁, d₂).
* *DF in rank 3.* A is DF iff some H has |A^H| − 1 ∈ {d₁, d₂}.

*Lemma (d₁ = 1 ⟹ DF).* The multirestriction to any H has a nonzero derivation θ(u) = Mu of
degree 1.
* If θ is a multiple of E, every m_P = 1, so |A^H| = |A| − 1 = d₂ + 1.
* Otherwise M is not scalar, and every α_P is an eigenvector of Mᵀ or is killed by θ. A non-scalar
  2×2 matrix has at most two eigenlines, so |A^H| ≤ 2. Essentiality gives |A^H| ≥ 2, so
  |A^H| = 2 = d₁ + 1. ∎

*Theorem (q = 2).* By the Frobenius bound d₁ ≤ 2. If d₁ = 1, use the lemma. So let d₁ = 2 and
suppose A is not DF.
1. Lines of PG(2, 2) have 3 points, and 3 = d₁ + 1 is excluded, so |A^H| = 2 for every H.
2. Then the Ziegler restriction is Boolean with multiplicities {m₁, m₂} = {2, d₂}.
3. A point of PG(2, 2) lies on at most 3 lines, so m ≤ 2. Hence d₂ = 2 and |A| = 5.
4. So every line contains exactly two points, each on exactly 3 lines. Counting incidences between
   lines and triple points gives 3·t₃ = 2·5 = 10, which is impossible. ∎

## 7. Verification

* **Lemma 1, exhaustively** (`work/covering_m2_all.py`, Macaulay2). Every essential covering
  class is free with exactly the predicted exponents:
  - all 3,097 in PG(3, 3);
  - all 38 in PG(4, 2).

  Random tests with the exact engine also agree: 40/40 in PG(2, 5) and 51/51 in PG(3, 2)
  (`work/covering_lemma_test*.py`).
* **The certificate of §2** returns true for the free 22-plane example and false for a non-free
  16-plane control (`work/membership_validate.m2`).
* **Cross-check of the PG(3, 3) search.** The unfiltered run (`work/free_df_parallel.py`, exact
  freeness for all 2,877 non-DF split-χ classes) found the same 7 classes. Among those 2,877
  classes, the covering ones are exactly the 7 free ones. In PG(4, 2), none of the 13 non-free
  non-DF classes is covering.
* **A_q.**
  - Determinant identity checked symbolically for q = 4, 5, 7.
  - Freeness checked with the exact engine for q = 4, 5, 7, 8, 9, 11, 13, 16.
  - Freeness checked with Macaulay2 for q = 4, 5, 7, 8, 9, 11.
  - Script: `work/Aq_verify.py`.
* **§4.2.** Free with generator degrees (1, 3, 9, 9) by the exact engine and by Macaulay2 over
  ZZ/3, with restriction sizes {10, 11, 12} (`work/f3_hit_analyze.py`, `work/f3_quadric_model.py`).
  The model is projectively equivalent to the class found by the search.
* **§4.3.** Combinatorial checks (`work/f2_rank6_cap.py`, `work/f2_cap_basis.py`):
  - K is a complete cap;
  - χ(A) = (t−1)(t−2)(t−4)(t−8)(t−16)(t−19);
  - all 50 restrictions have 28 or 30 hyperplanes, and none has χ dividing χ(A).
  - Freeness, independently of the deletion theorem: the certificate of §2 holds. Macaulay2 finds
    θ_{F⁵} in the span of E, θ_F, θ_{F²}, θ_{F³}, θ_{F⁴} modulo Q_K, in 1,951 s
    (`work/f2_rank6_membership.m2`).
* **Exhaustive tables.** `work/exhaustive_free_df.py` (F₂ ranks 3–5, F₃ rank 3) and
  `work/find_hit_3_4.py`, `work/find_all_3_4.py` (F₃ rank 4). The 7 classes are stored in
  `work/f3_rank4_free_nondf.json`.

## 8. Context

Over F₂ and F₃, Terao's conjecture holds trivially, because the lattice determines the
arrangement. Over any field it holds for DF arrangements, because Abe's DF is combinatorial.
Theorem 2 shows that the DF route cannot reach every free arrangement over any finite field.

Even so, the freeness in Theorem 2 is still combinatorial within F_q: both hypotheses can be read
off the lattice (and q), as the remark after Lemma 1 shows.

**No Terao pairs from the F₃ examples.** I checked the restrictions of all 7 PG(3, 3) classes
(`work/f3_rigidity_f9.py`):
* each has exactly one frame-normalised realisation over GF(3), GF(9) and GF(27);
* none has a realisation over GF(4), GF(5), GF(7), GF(8), GF(13), GF(16) or GF(25).

So these lattices are rigid and appear to be characteristic-3 only. They cannot give a Terao
counterexample.
