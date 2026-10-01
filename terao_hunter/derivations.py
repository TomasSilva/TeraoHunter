"""
derivations.py -- the gauge-reduced logarithmic derivation system.

Mathematics
-----------
Let A have lines alpha_0..alpha_{n-1}, Q = prod alpha_i, and fix one line
H = alpha_h.  Put

    D_H(A) = { theta in D(A) : theta(alpha_h) = 0 }.

Then D(A) = S.E (+) D_H(A)   (theta -> theta - (theta(alpha_h)/alpha_h) E),
and D_H(A) ~ D_0(A) = { theta : theta(Q) = 0 } = Jacobian syzygies, as graded
modules, via

    theta  |->  sigma = theta - (theta(Q) / (n Q)) E,
    theta(Q)/Q = sum_i theta(alpha_i)/alpha_i   (a polynomial).

sigma is nonzero whenever theta is: theta = gE with theta(alpha_h) = 0 forces
g = 0.  In particular an Euler multiple can never appear as a syzygy.

Consequences used by the search:
  * mdr(Q) = min { r : D_H(A)_r != 0 };
  * D_H(A) is torsion free, so D_H(A)_r = 0 implies D_H(A)_s = 0 for s < r;
  * A free with exponents (1, d1, d2)  <=>  D_H(A) is free on generators of
    degrees d1, d2.

Coordinates: if alpha_h has nonzero coefficient in the variable x_e, then
theta_e = -(l_h[u] theta_u + l_h[v] theta_v) / l_h[e], and the unknowns are
the two polynomials (theta_u, theta_v) of degree r: (r+1)(r+2) columns.
For every other line the condition "alpha_i divides theta(alpha_i)" is
imposed by vanishing at r+1 distinct points of the line, giving
(n-1)(r+1) rows.  The points P_m = P1 + m P2 are defined by an EXACT rule,
so the F_p matrix is the reduction of the characteristic-zero matrix.
"""

from __future__ import annotations

from typing import List, Sequence, Tuple

import numpy as np

from .arrangement import Arrangement
from .linalg import BadPrime, reduce_scalar
from .polys import Poly, monomials, p_add, p_linear, p_scale


def gauge_indices(arr: Arrangement, h: int) -> Tuple[int, int, int]:
    """(e, u, v): eliminated component e (l_h[e] != 0) and free components."""
    lh = arr.lines[h]
    e = max(i for i in range(3) if lh[i] != 0)   # prefer z, then y, then x
    u, v = [i for i in range(3) if i != e]
    return e, u, v


def _line_points_rule(l: Sequence) -> Tuple[int, int, int]:
    """Exact rule choosing two independent points of the line l:
    k = first index with l[k] != 0; P1 = l x e_a, P2 = l x e_b, {a,b} = rest."""
    k = next(i for i in range(3) if l[i] != 0)
    a, b = [i for i in range(3) if i != k]
    return k, a, b


def _cross_unit(l, a):
    """l x e_a."""
    if a == 0:
        return (0 * l[0], l[2], -l[1])
    if a == 1:
        return (-l[2], 0 * l[0], l[0])
    return (l[1], -l[0], 0 * l[0])


def reduced_coefficients(arr: Arrangement, h: int):
    """For every line i != h the pair (A_i, B_i) with
    l_h[e] * theta(alpha_i) = A_i theta_u + B_i theta_v  on D_H."""
    e, u, v = gauge_indices(arr, h)
    lh = arr.lines[h]
    out = []
    for i, l in enumerate(arr.lines):
        if i == h:
            continue
        out.append((i, lh[e] * l[u] - l[e] * lh[u],
                    lh[e] * l[v] - l[e] * lh[v]))
    return out


def exact_gauge_matrix(arr: Arrangement, h: int, r: int) -> List[list]:
    """Characteristic-zero matrix (reference implementation, small r)."""
    mons = monomials(r)
    rows = []
    for (i, Ai, Bi) in reduced_coefficients(arr, h):
        l = arr.lines[i]
        _k, a, b = _line_points_rule(l)
        P1, P2 = _cross_unit(l, a), _cross_unit(l, b)
        for m in range(r + 1):
            P = tuple(P1[t] + m * P2[t] for t in range(3))
            vals = [(P[0] ** mi) * (P[1] ** mj) * (P[2] ** mk)
                    for (mi, mj, mk) in mons]
            rows.append([Ai * val for val in vals] + [Bi * val for val in vals])
    return rows


def gauge_matrix_mod_p(arr: Arrangement, h: int, r: int, p: int, root: int
                       ) -> np.ndarray:
    """The same matrix reduced mod (p, s - root), built directly in F_p."""
    if r + 1 >= p:
        raise BadPrime(p)
    L = np.array([[reduce_scalar(c, p, root) for c in l] for l in arr.lines],
                 dtype=np.int64)
    e, u, v = gauge_indices(arr, h)
    lh = L[h]
    others = [i for i in range(arr.n) if i != h]
    Lo = L[others]
    Acoef = (lh[e] * Lo[:, u] % p - Lo[:, e] * lh[u] % p) % p
    Bcoef = (lh[e] * Lo[:, v] % p - Lo[:, e] * lh[v] % p) % p
    P1 = np.zeros_like(Lo)
    P2 = np.zeros_like(Lo)
    for t, i in enumerate(others):
        _k, a, b = _line_points_rule(arr.lines[i])
        P1[t] = [int(x) % p for x in _cross_unit([int(c) for c in L[i]], a)]
        P2[t] = [int(x) % p for x in _cross_unit([int(c) for c in L[i]], b)]
    ms = np.arange(r + 1, dtype=np.int64)
    # points: shape (n-1, r+1, 3)
    pts = (P1[:, None, :] + ms[None, :, None] * P2[:, None, :]) % p
    pts = pts.reshape(-1, 3)
    npts = pts.shape[0]
    powers = np.ones((npts, 3, r + 1), dtype=np.int64)
    for k in range(1, r + 1):
        powers[:, :, k] = powers[:, :, k - 1] * pts % p
    mons = np.array(monomials(r), dtype=np.int64)
    vals = powers[:, 0, mons[:, 0]] * powers[:, 1, mons[:, 1]] % p
    vals = vals * powers[:, 2, mons[:, 2]] % p
    Arep = np.repeat(Acoef, r + 1)[:, None]
    Brep = np.repeat(Bcoef, r + 1)[:, None]
    return np.concatenate([Arep * vals % p, Brep * vals % p], axis=1)


def gauge_vector_to_theta(arr: Arrangement, h: int, r: int, vec: Sequence
                          ) -> Tuple[Poly, Poly, Poly]:
    """(theta_x, theta_y, theta_z) in the ORIGINAL coordinates."""
    mons = monomials(r)
    N = len(mons)
    e, u, v = gauge_indices(arr, h)
    f = {m: c for m, c in zip(mons, vec[:N]) if c != 0}
    g = {m: c for m, c in zip(mons, vec[N:2 * N]) if c != 0}
    lh = arr.lines[h]
    inv = -1 / lh[e]
    third = p_add(p_scale(f, lh[u] * inv), p_scale(g, lh[v] * inv))
    comps = [None, None, None]
    comps[u], comps[v], comps[e] = f, g, third
    return tuple(comps)


def apply_theta_to_line(theta: Sequence[Poly], l: Sequence) -> Poly:
    out: Poly = {}
    for comp, c in zip(theta, l):
        if c != 0:
            out = p_add(out, p_scale(comp, c))
    return out


def euler() -> Tuple[Poly, Poly, Poly]:
    return (p_linear(1, 0, 0), p_linear(0, 1, 0), p_linear(0, 0, 1))


def choose_gauge_line(arr: Arrangement) -> int:
    """Any line works; take one with the most intersection points (its
    restriction conditions are the ones we drop)."""
    counts = arr.line_point_counts()
    return max(range(arr.n), key=lambda i: (counts[i], -i))
