"""
polys.py -- sparse exact polynomials in (x, y, z) over a number field.

A polynomial is a dict {(i, j, k): coefficient} with exact scalar
coefficients (Fraction / NFElem); zero coefficients are never stored, so
`not p` is an exact zero test.  Everything here is deliberately elementary:
the independent verifier does NOT import this module (it uses sympy), so the
two paths share no polynomial code.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Dict, List, Sequence, Tuple

Poly = Dict[Tuple[int, int, int], object]

ONE = Fraction(1)
ZERO = Fraction(0)


def monomials(d: int) -> List[Tuple[int, int, int]]:
    """Monomials of degree d, fixed order: i descending, then j descending."""
    if d < 0:
        return []
    return [(i, j, d - i - j) for i in range(d, -1, -1)
            for j in range(d - i, -1, -1)]


def dim_S(d: int) -> int:
    return (d + 1) * (d + 2) // 2 if d >= 0 else 0


def p_add(a: Poly, b: Poly) -> Poly:
    out = dict(a)
    for m, c in b.items():
        v = out.get(m, ZERO) + c
        if v == 0:
            out.pop(m, None)
        else:
            out[m] = v
    return out


def p_scale(a: Poly, s) -> Poly:
    if s == 0:
        return {}
    return {m: c * s for m, c in a.items()}


def p_sub(a: Poly, b: Poly) -> Poly:
    return p_add(a, p_scale(b, -ONE))


def p_mul(a: Poly, b: Poly) -> Poly:
    out: Poly = {}
    if len(a) > len(b):
        a, b = b, a
    for (i1, j1, k1), c1 in a.items():
        for (i2, j2, k2), c2 in b.items():
            m = (i1 + i2, j1 + j2, k1 + k2)
            out[m] = out.get(m, ZERO) + c1 * c2
    return {m: c for m, c in out.items() if c != 0}


def p_linear(a, b, c) -> Poly:
    out = {}
    for m, v in (((1, 0, 0), a), ((0, 1, 0), b), ((0, 0, 1), c)):
        if v != 0:
            out[m] = v
    return out


def p_degree(a: Poly):
    """Degree if homogeneous, None for the zero polynomial; raises if the
    polynomial is not homogeneous."""
    degs = {sum(m) for m in a}
    if not degs:
        return None
    if len(degs) != 1:
        raise ValueError("polynomial is not homogeneous")
    return degs.pop()


def p_diff(a: Poly, var: int) -> Poly:
    out = {}
    for m, c in a.items():
        e = m[var]
        if e:
            mm = list(m)
            mm[var] -= 1
            out[tuple(mm)] = c * e
    return out


def p_eval(a: Poly, pt: Sequence) -> object:
    x, y, z = pt
    acc = ZERO
    for (i, j, k), c in a.items():
        acc = acc + c * (x ** i) * (y ** j) * (z ** k)
    return acc


def p_product_of_lines(lines: Sequence[Sequence]) -> Poly:
    q: Poly = {(0, 0, 0): ONE}
    for (a, b, c) in lines:
        q = p_mul(q, p_linear(a, b, c))
    return q


def p_substitute_linear(a: Poly, forms: Sequence[Poly]) -> Poly:
    """a(F0, F1, F2) for three linear forms (as Polys).  Powers are cached."""
    d = p_degree(a)
    if d is None:
        return {}
    pows = []
    for F in forms:
        lst = [{(0, 0, 0): ONE}]
        for _ in range(d):
            lst.append(p_mul(lst[-1], F))
        pows.append(lst)
    cache: Dict[Tuple[int, int], Poly] = {}
    out: Poly = {}
    for (i, j, k), c in a.items():
        key = (i, j)
        if key not in cache:
            cache[key] = p_mul(pows[0][i], pows[1][j])
        out = p_add(out, p_scale(p_mul(cache[key], pows[2][k]), c))
    return out


def p_divide_by_linear(a: Poly, line: Sequence):
    """Exact division of a homogeneous polynomial by the linear form
    l = a0 x + a1 y + a2 z.  Returns (quotient, remainder_is_zero).

    Lex-order division with respect to a variable in which l has a nonzero
    coefficient; exact arithmetic throughout.
    """
    var = next(i for i in range(3) if line[i] != 0)
    lead = line[var]
    rest = [(i, line[i]) for i in range(3) if i != var and line[i] != 0]
    work = dict(a)
    quo: Poly = {}
    # process monomials by decreasing exponent of `var`
    while work:
        m = max(work, key=lambda t: t[var])
        if m[var] == 0:
            return quo, False
        c = work.pop(m)
        qm = list(m)
        qm[var] -= 1
        qm = tuple(qm)
        qc = c / lead
        quo[qm] = quo.get(qm, ZERO) + qc
        for (i, li) in rest:
            mm = list(qm)
            mm[i] += 1
            mm = tuple(mm)
            v = work.get(mm, ZERO) - qc * li
            if v == 0:
                work.pop(mm, None)
            else:
                work[mm] = v
    return {m: c for m, c in quo.items() if c != 0}, True


def p_restrict_to_line(a: Poly, line: Sequence) -> Dict[Tuple[int, int], object]:
    """Restriction of a to the line l = 0 as a binary form in the two
    remaining variables (the variable with nonzero coefficient is solved
    for).  The restriction is zero iff l divides a."""
    var = next(i for i in range(3) if line[i] != 0)
    others = [i for i in range(3) if i != var]
    lead = line[var]
    # var = -(l_o0 * o0 + l_o1 * o1) / lead
    c0 = -line[others[0]] / lead
    c1 = -line[others[1]] / lead
    d = p_degree(a)
    if d is None:
        return {}
    # binomial expansion of (c0 u + c1 v)^e, cached
    from math import comb
    out: Dict[Tuple[int, int], object] = {}
    pw0 = [ONE]
    pw1 = [ONE]
    for _ in range(d):
        pw0.append(pw0[-1] * c0)
        pw1.append(pw1[-1] * c1)
    for m, c in a.items():
        e = m[var]
        u, v = m[others[0]], m[others[1]]
        for t in range(e + 1):
            coef = c * comb(e, t) * pw0[t] * pw1[e - t]
            if coef == 0:
                continue
            key = (u + t, v + e - t)
            val = out.get(key, ZERO) + coef
            if val == 0:
                out.pop(key, None)
            else:
                out[key] = val
    return out


def p_to_json(a: Poly, K) -> list:
    return [[list(m), K.elem_to_json(c)] for m, c in sorted(a.items())]


def p_from_json(data, K) -> Poly:
    out = {}
    for m, c in data:
        v = K.elem_from_json(c)
        if v != 0:
            out[tuple(int(e) for e in m)] = v
    return out


def det3_polys(cols: Sequence[Sequence[Poly]]) -> Poly:
    """Determinant of the 3x3 matrix whose COLUMNS are the given triples."""
    (a, d, g), (b, e, h), (c, f, i) = cols  # column-major naming
    # matrix rows: [a b c; d e f; g h i]
    t1 = p_mul(a, p_sub(p_mul(e, i), p_mul(f, h)))
    t2 = p_mul(b, p_sub(p_mul(d, i), p_mul(f, g)))
    t3 = p_mul(c, p_sub(p_mul(d, h), p_mul(e, g)))
    return p_add(p_sub(t1, t2), t3)
