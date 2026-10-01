"""
arrangement.py -- exact projective line arrangements in CP^2.

Conventions
-----------
* A line is the coefficient triple (a, b, c) of a x + b y + c z over an exact
  number field, normalised so that the first nonzero coordinate is 1.
  Equality of normalised triples is equality of projective lines.
* A rank-two flat is the set of ALL lines through one intersection point,
  double points included.  `flats()` returns the full labelled list.
* b2 = sum_p (m_p - 1) is the coefficient of t in
        chi(A, t) = t^3 - n t^2 + b2 t - (b2 - n + 1)
                  = (t - 1) (t^2 - (n - 1) t + (b2 - (n - 1))).
  The constant coefficient of the REDUCED polynomial is  c0 = b2 - (n - 1);
  for target exponents (1, d1, d2):  c0* = d1 d2  and  b2* = (n - 1) + d1 d2.
  The two numbers are kept apart everywhere (`b2`, `reduced_constant`).

Admissibility (the search rejects everything else): at least 3 lines, no
repeated line, essential (no common point), and not a near-pencil.
"""

from __future__ import annotations

from fractions import Fraction
from math import isqrt
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .fields import NFElem, NumberField, QQ, field_of, is_exact_scalar

Line = Tuple[object, object, object]


def normalise(v: Sequence) -> Tuple:
    """Projective normal form: first nonzero coordinate equal to 1."""
    for c in v:
        if c != 0:
            if c == 1:
                return tuple(v)
            inv = 1 / c if not isinstance(c, int) else Fraction(1, c)
            return tuple(x * inv for x in v)
    raise ValueError("zero vector is not a projective point/line")


def cross(u: Sequence, v: Sequence) -> Tuple:
    return (u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0])


def det3(u, v, w):
    return (u[0] * (v[1] * w[2] - v[2] * w[1])
            - u[1] * (v[0] * w[2] - v[2] * w[0])
            + u[2] * (v[0] * w[1] - v[1] * w[0]))


def dot(u, v):
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


class InvalidArrangement(ValueError):
    pass


class Arrangement:
    """Immutable exact line arrangement.  Line order is meaningful (labels)."""

    def __init__(self, lines: Iterable[Sequence], field: Optional[NumberField]
                 = None):
        raw = [tuple(l) for l in lines]
        for l in raw:
            if len(l) != 3 or not all(is_exact_scalar(c) for c in l):
                raise InvalidArrangement("inexact or malformed line data")
        flat = [c for l in raw for c in l]
        K = field_of(flat)
        if field is not None and not field.is_rational:
            if not K.is_rational and K != field:
                raise InvalidArrangement("coordinates not in declared field")
            K = field
        self.field: NumberField = K
        self.lines: Tuple[Line, ...] = tuple(
            normalise(tuple(K.coerce(c) for c in l)) for l in raw)
        self._points: Optional[Dict[Tuple, Tuple[int, ...]]] = None

    # -- basic structure --------------------------------------------------
    def __len__(self):
        return len(self.lines)

    @property
    def n(self) -> int:
        return len(self.lines)

    def has_repeated_line(self) -> bool:
        return len(set(self.lines)) != len(self.lines)

    def points(self) -> Dict[Tuple, Tuple[int, ...]]:
        """{normalised point: sorted tuple of line labels through it}."""
        if self._points is None:
            if self.has_repeated_line():
                raise InvalidArrangement("repeated line")
            acc: Dict[Tuple, set] = {}
            L = self.lines
            for i in range(len(L)):
                for j in range(i + 1, len(L)):
                    p = normalise(cross(L[i], L[j]))
                    acc.setdefault(p, set()).update((i, j))
            self._points = {p: tuple(sorted(s)) for p, s in acc.items()}
        return self._points

    def flats(self) -> List[Tuple[int, ...]]:
        """Full labelled list of rank-two flats (double points included)."""
        return sorted(self.points().values())

    def multiplicities(self) -> Dict[int, int]:
        t: Dict[int, int] = {}
        for f in self.points().values():
            t[len(f)] = t.get(len(f), 0) + 1
        return dict(sorted(t.items()))

    def m_max(self) -> int:
        return max(len(f) for f in self.points().values())

    def b2(self) -> int:
        return sum(len(f) - 1 for f in self.points().values())

    def reduced_constant(self) -> int:
        """Constant coefficient of chi(t)/(t-1); equals d1*d2 when free."""
        return self.b2() - (self.n - 1)

    def line_point_counts(self) -> List[int]:
        """s_H for every line: number of DISTINCT intersection points on H."""
        cnt = [0] * self.n
        for f in self.points().values():
            for i in f:
                cnt[i] += 1
        return cnt

    # -- admissibility ----------------------------------------------------
    def admissibility(self) -> Tuple[bool, str]:
        if self.n < 3:
            return False, "fewer_than_three_lines"
        if self.has_repeated_line():
            return False, "repeated_line"
        m = self.m_max()
        if m == self.n:
            return False, "pencil_nonessential"
        if m == self.n - 1:
            return False, "near_pencil"
        return True, "ok"

    def is_admissible(self) -> bool:
        return self.admissibility()[0]

    # -- exponents --------------------------------------------------------
    def candidate_exponents(self) -> Optional[Tuple[int, int]]:
        """(d1, d2), d1 <= d2, with d1 + d2 = n - 1 and d1 d2 = b2 - (n-1);
        None when the reduced characteristic polynomial has no such roots
        (then A is not free, by Terao's factorisation theorem)."""
        return exponents_from_b2(self.n, self.b2())

    # -- transformations --------------------------------------------------
    def permuted(self, perm: Sequence[int]) -> "Arrangement":
        """New arrangement whose line k is the old line perm[k]."""
        return Arrangement([self.lines[i] for i in perm], self.field)

    def transformed(self, T: Sequence[Sequence]) -> "Arrangement":
        """Image under the projective map x -> T x (T invertible 3x3).
        Lines are covectors: l -> l T^{-1}; we use adj(T) (same up to
        scale)."""
        adj = _adjugate(T)
        if det3(*T) == 0:
            raise ValueError("singular projective transformation")
        new = [tuple(sum(l[i] * adj[i][j] for i in range(3))
                     for j in range(3)) for l in self.lines]
        return Arrangement(new, self.field)

    def replace_lines(self, remove: Sequence[int], add: Sequence[Sequence]
                      ) -> "Arrangement":
        keep = [l for i, l in enumerate(self.lines) if i not in set(remove)]
        return Arrangement(keep + [tuple(a) for a in add], self.field)

    def delete_line(self, i: int) -> "Arrangement":
        return self.replace_lines([i], [])

    def conjugated(self) -> "Arrangement":
        """Galois conjugate (quadratic fields); same abstract field."""
        K = self.field
        return Arrangement([tuple(K.conjugate(c) for c in l)
                            for l in self.lines], K)

    # -- numerics (diagnostics only) ---------------------------------------
    def embedded(self):
        import numpy as np
        K = self.field
        arr = np.array([[K.embed(c) for c in l] for l in self.lines],
                       dtype=np.complex128)
        if K.is_real_embedded:
            return arr.real.copy()
        return arr

    def height(self) -> int:
        h = 0
        for l in self.lines:
            for c in l:
                for q in (c.coeffs if isinstance(c, NFElem) else (c,)):
                    h = max(h, abs(q.numerator), q.denominator)
        return h

    # -- serialisation ----------------------------------------------------
    def to_json(self) -> dict:
        K = self.field
        return {"field": K.to_json() if not K.is_rational else "QQ",
                "lines": [[K.elem_to_json(c) for c in l] for l in self.lines]}

    @staticmethod
    def from_json(d: dict) -> "Arrangement":
        K = NumberField.from_json(d.get("field", "QQ"))
        return Arrangement([[K.elem_from_json(c) for c in l]
                            for l in d["lines"]], K)

    def key(self) -> str:
        """Coordinate-level identity of the unordered line set."""
        K = self.field
        body = sorted(str([K.elem_to_json(c) for c in l]) for l in self.lines)
        return str((K.minpoly, body))


def _adjugate(T):
    (a, b, c), (d, e, f), (g, h, i) = T
    return ((e * i - f * h, c * h - b * i, b * f - c * e),
            (f * g - d * i, a * i - c * g, c * d - a * f),
            (d * h - e * g, b * g - a * h, a * e - b * d))


def exponents_from_b2(n: int, b2: int) -> Optional[Tuple[int, int]]:
    c0 = b2 - (n - 1)
    disc = (n - 1) ** 2 - 4 * c0
    if disc < 0 or c0 < 0:
        return None
    r = isqrt(disc)
    if r * r != disc or (n - 1 - r) % 2:
        return None
    d1 = (n - 1 - r) // 2
    return (d1, n - 1 - d1)


def b2_target(n: int, d1: int, d2: int) -> int:
    """b2* = (n - 1) + d1 d2  (NOT the reduced constant coefficient d1 d2)."""
    if not (1 <= d1 <= d2 and d1 + d2 == n - 1):
        raise ValueError("need 1 <= d1 <= d2 and d1 + d2 = n - 1")
    return (n - 1) + d1 * d2


def line_through(p: Sequence, q: Sequence) -> Tuple:
    return normalise(cross(p, q))
