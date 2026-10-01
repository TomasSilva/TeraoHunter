"""
fields.py -- exact characteristic-zero number fields embedded in C.

A field is K = Q[s]/(f) with f monic, irreducible, integer coefficients,
together with ONE fixed complex embedding s -> root (a numerical label that
selects the embedding; it is never used in exact computations).

Scalars
-------
* Elements of Q are `fractions.Fraction`.
* Elements of a field of degree m >= 2 are `NFElem` (tuple of m Fractions in
  the power basis 1, s, ..., s^(m-1)).  An NFElem whose higher coefficients
  vanish compares and hashes equal to the corresponding Fraction, so mixed
  containers are consistent and `v == 0` is an exact zero test.

The six project fields are presets (`PRESETS`), but nothing below is
specific to them: any monic irreducible integer polynomial defines a field,
and `NumberField.from_minpoly` is what the realisation module uses when an
incidence system has solutions in an extension.  Restricting a search to the
presets is a search limitation, not a mathematical exclusion.

Exact numerical data (embeddings) are kept strictly apart: `embed` returns a
Python complex and is used only by numerical diagnostics.
"""

from __future__ import annotations

import cmath
import math
import random
from fractions import Fraction
from typing import Iterable, List, Optional, Sequence, Tuple

Scalar = "Fraction | NFElem"


# --------------------------------------------------------------------------
# small integer-polynomial helpers (coefficients low -> high)
# --------------------------------------------------------------------------

def _trim(p: list) -> list:
    while p and p[-1] == 0:
        p.pop()
    return p


def _poly_roots_complex(coeffs: Sequence[int]) -> List[complex]:
    """All complex roots (numpy companion matrix); numerical labels only."""
    import numpy as np
    return [complex(r) for r in np.roots(list(reversed(coeffs)))]


def _polymulmod_p(a, b, f, p):
    """a*b mod (f, p); f monic, lists low->high."""
    m = len(f) - 1
    res = [0] * (len(a) + len(b) - 1) if a and b else []
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                res[i + j] = (res[i + j] + ai * bj) % p
    for k in range(len(res) - 1, m - 1, -1):
        c = res[k]
        if c:
            for j in range(m + 1):
                res[k - m + j] = (res[k - m + j] - c * f[j]) % p
    return _trim(res[:m])


def _polypowmod_p(base, e, f, p):
    result = [1]
    base = list(base)
    while e:
        if e & 1:
            result = _polymulmod_p(result, base, f, p)
        base = _polymulmod_p(base, base, f, p)
        e >>= 1
    return result


def _polydivmod_p(a, b, p):
    a = [x % p for x in a]
    _trim(a)
    b = _trim([x % p for x in b])
    inv = pow(b[-1], p - 2, p)
    q = [0] * max(0, len(a) - len(b) + 1)
    while len(a) >= len(b) and a:
        c = a[-1] * inv % p
        k = len(a) - len(b)
        q[k] = c
        for j, bj in enumerate(b):
            a[k + j] = (a[k + j] - c * bj) % p
        _trim(a)
    return q, a


def _polygcd_p(a, b, p):
    a, b = _trim([x % p for x in a]), _trim([x % p for x in b])
    while b:
        a, b = b, _polydivmod_p(a, b, p)[1]
    if a:
        inv = pow(a[-1], p - 2, p)
        a = [x * inv % p for x in a]
    return a


def roots_mod_p(f: Sequence[int], p: int, rng: Optional[random.Random] = None
                ) -> List[int]:
    """Distinct roots of the integer polynomial f in F_p (p an odd prime).

    gcd(f, x^p - x) isolates the split part; Cantor--Zassenhaus equal-degree
    splitting finds the roots.  Deterministic given `rng` seed.
    """
    rng = rng or random.Random(p)
    fp = _trim([c % p for c in f])
    if len(fp) <= 1:
        return []
    inv = pow(fp[-1], p - 2, p)
    fp = [c * inv % p for c in fp]
    if len(fp) == 2:
        return [(-fp[0]) % p]
    xp = _polypowmod_p([0, 1], p, fp, p)
    diff = list(xp) + [0] * max(0, 2 - len(xp))
    diff[1] = (diff[1] - 1) % p
    g = _polygcd_p(fp, diff, p)
    out: List[int] = []

    def split(h):
        if len(h) <= 1:
            return
        if len(h) == 2:
            out.append((-h[0]) % p)
            return
        while True:
            a = rng.randrange(p)
            t = _polypowmod_p([a, 1], (p - 1) // 2, h, p)
            t = list(t) + [0] * max(0, 1 - len(t))
            t[0] = (t[0] - 1) % p
            d = _polygcd_p(h, t, p)
            if 1 < len(d) < len(h):
                split(d)
                split(_polydivmod_p(h, d, p)[0])
                return

    split(g)
    return sorted(set(out))


# --------------------------------------------------------------------------
# number fields
# --------------------------------------------------------------------------

class NumberField:
    """Q[s]/(f) with a fixed complex embedding.  Degree 1 is Q itself."""

    def __init__(self, minpoly: Sequence[int], embedding: Optional[complex]
                 = None, name: Optional[str] = None, check: bool = True):
        f = [int(c) for c in minpoly]
        if len(f) < 2 or f[-1] != 1:
            raise ValueError("minpoly must be monic of degree >= 1 "
                             "(integer coefficients, low -> high)")
        self.minpoly: Tuple[int, ...] = tuple(f)
        self.degree = len(f) - 1
        if self.degree == 1:
            if f[0] != 0:
                raise ValueError("use minpoly [0, 1] for Q")
            self.embedding = 0j
        else:
            if check and not _is_irreducible_over_Q(f):
                raise ValueError(f"minpoly {f} is reducible over Q")
            roots = _poly_roots_complex(f)
            if embedding is None:
                # principal choice: largest imaginary part, then largest real
                embedding = max(roots, key=lambda r: (round(r.imag, 9),
                                                      round(r.real, 9)))
            best = min(roots, key=lambda r: abs(r - embedding))
            if abs(best - embedding) > 1e-6 * (1 + abs(best)):
                raise ValueError("embedding is not close to a root of minpoly")
            self.embedding = _polish_root(f, best)
        self.name = name or ("QQ" if self.degree == 1 else
                             "QQ[s]/(" + _poly_str(f) + ")")
        # reduction table: s^k for k in [m, 2m-2] in the power basis
        m = self.degree
        self._red: List[Tuple[Fraction, ...]] = []
        if m >= 2:
            cur = tuple(Fraction(-c) for c in f[:m])          # s^m
            self._red.append(cur)
            for _ in range(m - 2):
                top = cur[-1]
                shifted = (Fraction(0),) + cur[:-1]
                cur = tuple(shifted[i] + top * self._red[0][i]
                            for i in range(m))
                self._red.append(cur)

    # -- identity ---------------------------------------------------------
    def __eq__(self, other):
        return (isinstance(other, NumberField)
                and self.minpoly == other.minpoly
                and abs(self.embedding - other.embedding) < 1e-9)

    def __hash__(self):
        return hash(("NumberField", self.minpoly,
                     round(self.embedding.real, 6),
                     round(self.embedding.imag, 6)))

    def __repr__(self):
        return self.name

    @property
    def is_rational(self) -> bool:
        return self.degree == 1

    @property
    def is_real_embedded(self) -> bool:
        return abs(self.embedding.imag) < 1e-12

    # -- elements ---------------------------------------------------------
    def __call__(self, *coeffs) -> "Scalar":
        """K(a0, a1, ...) = a0 + a1 s + ...  (rationals / ints / strings)."""
        if len(coeffs) == 1 and isinstance(coeffs[0], (list, tuple)):
            coeffs = tuple(coeffs[0])
        cs = [_to_fraction(c) for c in coeffs]
        return self._make(cs)

    def _make(self, cs: Sequence[Fraction]):
        m = self.degree
        cs = list(cs)
        if len(cs) > m:
            if any(cs[m:]):
                raise ValueError("too many coefficients for this field")
            cs = cs[:m]
        if all(c == 0 for c in cs[1:]):
            return cs[0] if cs else Fraction(0)
        cs += [Fraction(0)] * (m - len(cs))
        return NFElem(self, tuple(cs))

    @property
    def gen(self):
        if self.degree == 1:
            raise ValueError("Q has no generator")
        return self._make([0, 1])

    def coerce(self, v) -> "Scalar":
        """Exact coercion; floats are rejected (no silent rationalisation)."""
        if isinstance(v, NFElem):
            if v.field != self:
                raise ValueError(f"element of {v.field} used in {self}")
            return v
        return _to_fraction(v)

    def contains(self, v) -> bool:
        if isinstance(v, NFElem):
            return v.field == self
        return isinstance(v, (int, Fraction))

    def embed(self, v) -> complex:
        if isinstance(v, NFElem):
            acc, pw = 0j, 1 + 0j
            for c in v.coeffs:
                acc += (c.numerator / c.denominator) * pw
                pw *= self.embedding
            return acc
        return complex(float(v))

    def all_embeddings(self) -> List[complex]:
        """All complex roots of the minimal polynomial (numerical labels)."""
        if self.degree == 1:
            return [0j]
        return [_polish_root(self.minpoly, r)
                for r in _poly_roots_complex(self.minpoly)]

    def with_embedding(self, root: complex) -> "NumberField":
        """Same abstract field, different complex embedding.

        Applying a different embedding to the same exact coordinates is how
        Galois-conjugate complex arrangements are represented: exact data
        (lattice, freeness, certificates) are embedding-independent.
        """
        return NumberField(self.minpoly, root, name=self.name, check=False)

    def conjugate(self, v):
        """Nontrivial Galois automorphism; quadratic fields only."""
        if not isinstance(v, NFElem):
            return v
        if self.degree != 2:
            raise NotImplementedError("explicit automorphisms are only "
                                      "implemented for quadratic fields")
        # s -> -b - s  for  s^2 + b s + c
        b = Fraction(self.minpoly[1])
        a0, a1 = v.coeffs
        return self._make([a0 - a1 * b, -a1])

    # -- modular reduction ------------------------------------------------
    def roots_mod(self, p: int) -> List[int]:
        if self.degree == 1:
            return [0]
        return roots_mod_p(self.minpoly, p)

    # -- serialisation ----------------------------------------------------
    def to_json(self) -> dict:
        return {"type": "number_field", "name": self.name,
                "minpoly": list(self.minpoly), "generator": "s",
                "embedding": [self.embedding.real, self.embedding.imag]}

    @staticmethod
    def from_json(tag) -> "NumberField":
        if tag in (None, "QQ"):
            return QQ
        if isinstance(tag, dict) and tag.get("type") == "number_field":
            emb = tag.get("embedding")
            emb = complex(emb[0], emb[1]) if emb else None
            K = NumberField(tag["minpoly"], emb, name=tag.get("name"))
            return QQ if K.degree == 1 else K
        if isinstance(tag, dict) and tag.get("type") == "quadratic":
            return quadratic_field(int(tag["d"]))       # legacy database tag
        raise ValueError(f"unrecognised field tag {tag!r}")

    @staticmethod
    def from_minpoly(minpoly: Sequence[int], embedding=None, name=None):
        return NumberField(minpoly, embedding, name)

    def elem_to_json(self, v):
        if isinstance(v, NFElem):
            return [str(c) for c in v.coeffs]
        return str(_to_fraction(v))

    def elem_from_json(self, d):
        if isinstance(d, list):
            return self._make([Fraction(c) for c in d])
        if isinstance(d, str):
            return Fraction(d)
        if isinstance(d, int):
            return Fraction(d)
        raise ValueError(f"inexact or malformed scalar {d!r}")


class NFElem:
    """Element of a NumberField of degree >= 2 (immutable, exact, hashable)."""

    __slots__ = ("field", "coeffs")

    def __init__(self, field: NumberField, coeffs: Tuple[Fraction, ...]):
        self.field = field
        self.coeffs = coeffs

    # -- helpers ----------------------------------------------------------
    def _other(self, o) -> Optional[Tuple[Fraction, ...]]:
        if isinstance(o, NFElem):
            if o.field != self.field:
                raise ValueError("mixing elements of different number fields")
            return o.coeffs
        if isinstance(o, (int, Fraction)):
            return (Fraction(o),) + (Fraction(0),) * (self.field.degree - 1)
        return None

    @property
    def is_rational(self):
        return all(c == 0 for c in self.coeffs[1:])

    # -- arithmetic -------------------------------------------------------
    def __add__(self, o):
        oc = self._other(o)
        if oc is None:
            return NotImplemented
        return self.field._make([a + b for a, b in zip(self.coeffs, oc)])

    __radd__ = __add__

    def __neg__(self):
        return NFElem(self.field, tuple(-a for a in self.coeffs))

    def __sub__(self, o):
        oc = self._other(o)
        if oc is None:
            return NotImplemented
        return self.field._make([a - b for a, b in zip(self.coeffs, oc)])

    def __rsub__(self, o):
        oc = self._other(o)
        if oc is None:
            return NotImplemented
        return self.field._make([b - a for a, b in zip(self.coeffs, oc)])

    def __mul__(self, o):
        if isinstance(o, (int, Fraction)):
            return self.field._make([a * o for a in self.coeffs])
        oc = self._other(o)
        if oc is None:
            return NotImplemented
        m = self.field.degree
        prod = [Fraction(0)] * (2 * m - 1)
        for i, a in enumerate(self.coeffs):
            if a:
                for j, b in enumerate(oc):
                    if b:
                        prod[i + j] += a * b
        res = prod[:m]
        for k in range(m, 2 * m - 1):
            c = prod[k]
            if c:
                red = self.field._red[k - m]
                for i in range(m):
                    res[i] += c * red[i]
        return self.field._make(res)

    __rmul__ = __mul__

    def inverse(self):
        """Exact inverse by solving the multiplication-matrix system."""
        m = self.field.degree
        K = self.field
        # columns: self * s^j
        cols = []
        cur = self
        for j in range(m):
            cc = cur.coeffs if isinstance(cur, NFElem) else \
                (Fraction(cur),) + (Fraction(0),) * (m - 1)
            cols.append(cc)
            cur = cur * K.gen
        A = [[cols[j][i] for j in range(m)] + [Fraction(int(i == 0))]
             for i in range(m)]
        for c in range(m):
            piv = next((r for r in range(c, m) if A[r][c] != 0), None)
            if piv is None:
                raise ZeroDivisionError("zero divisor: minpoly reducible?")
            A[c], A[piv] = A[piv], A[c]
            pv = A[c][c]
            A[c] = [v / pv for v in A[c]]
            for r in range(m):
                if r != c and A[r][c] != 0:
                    fac = A[r][c]
                    A[r] = [vr - fac * vc for vr, vc in zip(A[r], A[c])]
        return K._make([A[i][m] for i in range(m)])

    def __truediv__(self, o):
        if isinstance(o, (int, Fraction)):
            if o == 0:
                raise ZeroDivisionError
            return self * (Fraction(1) / Fraction(o))
        if isinstance(o, NFElem):
            return self * o.inverse()
        return NotImplemented

    def __rtruediv__(self, o):
        if isinstance(o, (int, Fraction)):
            return self.inverse() * o
        return NotImplemented

    def __pow__(self, n):
        if not isinstance(n, int):
            return NotImplemented
        if n < 0:
            return self.inverse() ** (-n)
        result, base = Fraction(1), self
        while n:
            if n & 1:
                result = base * result
            base = base * base
            n >>= 1
        return result

    # -- identity ---------------------------------------------------------
    def __eq__(self, o):
        if isinstance(o, NFElem):
            return o.field == self.field and o.coeffs == self.coeffs
        if isinstance(o, (int, Fraction)):
            return self.is_rational and self.coeffs[0] == o
        return NotImplemented

    def __hash__(self):
        if self.is_rational:
            return hash(self.coeffs[0])
        return hash((self.field.minpoly, self.coeffs))

    def __bool__(self):
        return any(c != 0 for c in self.coeffs)

    def __repr__(self):
        return "<" + ", ".join(str(c) for c in self.coeffs) + ">"

    def __complex__(self):
        return self.field.embed(self)

    def __float__(self):
        raise TypeError("refusing implicit float conversion of an exact "
                        "number-field element; use field.embed()")


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def _to_fraction(v) -> Fraction:
    if isinstance(v, Fraction):
        return v
    if isinstance(v, bool):
        raise TypeError("bool is not an exact scalar")
    if isinstance(v, int):
        return Fraction(v)
    if isinstance(v, str):
        return Fraction(v.strip())
    if isinstance(v, NFElem):
        if v.is_rational:
            return v.coeffs[0]
        raise ValueError("non-rational number-field element")
    if isinstance(v, float) or isinstance(v, complex):
        raise TypeError("floats are not exact scalars")
    num = getattr(v, "numerator", None)      # numpy ints, sympy Rationals (p/q)
    if hasattr(v, "p") and hasattr(v, "q"):
        return Fraction(int(v.p), int(v.q))
    if num is not None and not isinstance(v, float):
        return Fraction(int(v.numerator), int(v.denominator))
    raise TypeError(f"cannot coerce {type(v).__name__} to an exact scalar")


def _poly_str(f):
    terms = []
    for k in range(len(f) - 1, -1, -1):
        c = f[k]
        if c == 0:
            continue
        mon = "" if k == 0 else ("s" if k == 1 else f"s^{k}")
        coef = ("" if (abs(c) == 1 and k) else str(abs(c)))
        terms.append(("-" if c < 0 else "+") + coef + mon)
    s = "".join(terms)
    return s[1:] if s.startswith("+") else s


def _polish_root(f, r: complex, iters: int = 30) -> complex:
    for _ in range(iters):
        val = sum(c * r ** k for k, c in enumerate(f))
        der = sum(k * c * r ** (k - 1) for k, c in enumerate(f) if k)
        if der == 0:
            break
        step = val / der
        r = r - step
        if abs(step) < 1e-16 * (1 + abs(r)):
            break
    return complex(r)


def _is_irreducible_over_Q(f) -> bool:
    """Exact irreducibility.  Degree <= 3: rational-root test (complete).
    Higher degree: sympy's exact factorisation."""
    deg = len(f) - 1
    if deg == 1:
        return True
    if deg <= 3:
        c0 = f[0]
        if c0 == 0:
            return False
        for d in range(1, abs(c0) + 1):
            if c0 % d == 0:
                for r in (d, -d):
                    if sum(c * r ** k for k, c in enumerate(f)) == 0:
                        return False
        return True
    import sympy as sp
    x = sp.symbols("x")
    poly = sp.Poly(sum(int(c) * x ** k for k, c in enumerate(f)), x)
    return poly.is_irreducible


def quadratic_field(d: int) -> NumberField:
    """Q(sqrt d) with the project's principal embedding
    (sqrt d > 0 for d > 0, sqrt d = +i sqrt|d| for d < 0)."""
    d = int(d)
    emb = complex(math.sqrt(d), 0) if d > 0 else complex(0, math.sqrt(-d))
    return NumberField([-d, 0, 1], emb, name=f"QQ(sqrt{d})")


QQ = NumberField([0, 1], name="QQ")

PRESETS = {
    "QQ": QQ,
    "QQ(sqrt2)": quadratic_field(2),
    "QQ(sqrt3)": quadratic_field(3),
    "QQ(sqrt5)": quadratic_field(5),
    "QQ(i)": quadratic_field(-1),
    "QQ(sqrt-3)": quadratic_field(-3),
}
PRESETS["QQ(sqrt-1)"] = PRESETS["QQ(i)"]


def field_of(values: Iterable) -> NumberField:
    """The common field of a collection of scalars (QQ if all rational)."""
    K = QQ
    for v in values:
        if isinstance(v, NFElem):
            if K.is_rational:
                K = v.field
            elif K != v.field:
                raise ValueError("mixed number fields in one object")
    return K


def is_exact_scalar(v) -> bool:
    return isinstance(v, (Fraction, NFElem)) or (
        isinstance(v, int) and not isinstance(v, bool))
