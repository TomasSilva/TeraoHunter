"""
linalg.py -- exact linear algebra with a modular fast path.

SOUNDNESS OF THE MODULAR DIRECTION (tested in tests/test_linalg.py)
-------------------------------------------------------------------
Let M have entries in a number field K, p a prime, and phi: R -> F_p a ring
homomorphism defined on a subring R of K containing every entry of M
(R = Z_(p)[s], phi(s) = a root of the minimal polynomial mod p; this needs
every denominator to be prime to p -- otherwise `BadPrime`).  Then

    rank_{F_p} phi(M)  <=  rank_K M,

because a minor that is nonzero after phi is nonzero before.  Therefore:

  * FULL COLUMN RANK mod p  ==>  full column rank in characteristic zero
    (a proof; used for "no syzygy of degree r exists").
  * RANK DEFICIENCY mod p proves NOTHING about characteristic zero
    (p may be a bad prime).  It is only a hint telling us which exact kernel
    vector to reconstruct.
  * A reconstructed characteristic-zero kernel vector is accepted only after
    exact verification by the caller; together with nullity_p == (number of
    verified independent vectors) for one prime, the exact nullity is known.

Timeouts raise `Unresolved`; they are never converted into a mathematical
statement.
"""

from __future__ import annotations

import time
from fractions import Fraction
from math import isqrt
from typing import Callable, Dict, Iterator, List, Optional, Sequence, Tuple

import numpy as np

from .fields import NFElem, NumberField


class BadPrime(Exception):
    """The reduction homomorphism is not defined at this prime."""


class Unresolved(Exception):
    """A bounded computation did not finish; no mathematical conclusion."""


# --------------------------------------------------------------------------
# primes
# --------------------------------------------------------------------------

def is_prime(n: int) -> bool:
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for q in small:
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in small:                       # deterministic below 3.3e24
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def split_primes(K: NumberField, start: int = (1 << 31) - 1,
                 full: bool = True) -> Iterator[Tuple[int, List[int]]]:
    """Primes below `start` (descending) with the roots of minpoly mod p.

    full=True : only primes where the minimal polynomial has deg(K) distinct
                roots (needed to reconstruct elements of K).
    full=False: any prime with at least one root (enough for rank bounds).
    All primes are < 2^31 so products fit in int64.
    """
    p = start if start % 2 else start - 1
    while p > 1 << 20:
        if is_prime(p):
            roots = K.roots_mod(p)
            if (len(roots) == K.degree) if full else roots:
                yield p, roots
        p -= 2
    raise Unresolved("prime supply exhausted")


def reduce_scalar(v, p: int, root: int) -> int:
    """phi(v) in F_p; raises BadPrime if a denominator vanishes mod p."""
    if isinstance(v, NFElem):
        acc, pw = 0, 1
        for c in v.coeffs:
            if c:
                den = c.denominator % p
                if den == 0:
                    raise BadPrime(p)
                acc = (acc + c.numerator % p * pow(den, p - 2, p) % p * pw) % p
            pw = pw * root % p
        return acc
    c = v if isinstance(v, Fraction) else Fraction(v)
    den = c.denominator % p
    if den == 0:
        raise BadPrime(p)
    return c.numerator % p * pow(den, p - 2, p) % p


# --------------------------------------------------------------------------
# F_p elimination (numpy int64, p < 2^31)
# --------------------------------------------------------------------------

def rref_mod_p(A: np.ndarray, p: int) -> Tuple[np.ndarray, List[int]]:
    """Reduced row echelon form over F_p.  Returns (R, pivot_columns)."""
    A = np.array(A, dtype=np.int64) % p
    rows, cols = A.shape
    pivots: List[int] = []
    r = 0
    for c in range(cols):
        if r == rows:
            break
        nz = np.nonzero(A[r:, c])[0]
        if nz.size == 0:
            continue
        piv = r + int(nz[0])
        if piv != r:
            A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = A[r] * inv % p
        col = A[:, c].copy()
        col[r] = 0
        idx = np.nonzero(col)[0]
        if idx.size:
            A[idx] = (A[idx] - (col[idx, None] * A[r][None, :]) % p) % p
        pivots.append(c)
        r += 1
    return A[:r], pivots


def rank_mod_p(A: np.ndarray, p: int) -> int:
    return len(rref_mod_p(A, p)[1])


def kernel_mod_p(A: np.ndarray, p: int) -> Tuple[np.ndarray, List[int], List[int]]:
    """Canonical (RREF) kernel basis over F_p.

    Returns (basis, pivots, free_columns); basis[t] has 1 at free column t,
    0 at the other free columns.
    """
    R, pivots = rref_mod_p(A, p)
    cols = A.shape[1]
    pivset = set(pivots)
    free = [c for c in range(cols) if c not in pivset]
    basis = np.zeros((len(free), cols), dtype=np.int64)
    for t, fc in enumerate(free):
        basis[t, fc] = 1
        if pivots:
            basis[t, pivots] = (-R[:, fc]) % p
    return basis, pivots, free


# --------------------------------------------------------------------------
# rational reconstruction / CRT
# --------------------------------------------------------------------------

def rational_reconstruct(a: int, M: int) -> Optional[Fraction]:
    """n/d with n = a d (mod M), |n|, d <= sqrt(M/2); None if impossible."""
    a %= M
    bound = isqrt(M // 2)
    r0, r1 = M, a
    s0, s1 = 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    if s1 < 0:
        r1, s1 = -r1, -s1
    from math import gcd
    if gcd(r1, s1) != 1:
        return None
    return Fraction(r1, s1)


def _vandermonde_inverse(roots: Sequence[int], p: int) -> List[List[int]]:
    m = len(roots)
    A = [[pow(r, k, p) for k in range(m)] + [int(i == j) for j in range(m)]
         for i, r in enumerate(roots)]
    for c in range(m):
        piv = next(r for r in range(c, m) if A[r][c] % p)
        A[c], A[piv] = A[piv], A[c]
        inv = pow(A[c][c], p - 2, p)
        A[c] = [v * inv % p for v in A[c]]
        for r in range(m):
            if r != c and A[r][c]:
                f = A[r][c]
                A[r] = [(vr - f * vc) % p for vr, vc in zip(A[r], A[c])]
    return [row[m:] for row in A]


class KernelResult:
    def __init__(self, vectors, nullity_mod_p, free_columns, primes, selected):
        self.vectors = vectors              # list of lists of exact scalars
        self.nullity_mod_p = nullity_mod_p  # UPPER bound for the exact nullity
        self.free_columns = free_columns
        self.primes = primes
        self.selected = selected            # indices into the RREF basis


def modular_kernel(K: NumberField,
                   build: Callable[[int, int], np.ndarray],
                   select: Optional[Callable[[np.ndarray, int, int], Sequence[int]]]
                   = None,
                   verify: Optional[Callable[[List[list]], bool]] = None,
                   max_primes: int = 400, deadline: Optional[float] = None,
                   check_every: int = 2) -> KernelResult:
    """Reconstruct canonical kernel vectors over K by CRT over split primes.

    build(p, root)   -> the matrix reduced mod (p, s - root); may raise
                        BadPrime.
    select(basis, p, root) -> indices of the RREF-basis vectors wanted (chosen
                        on the first good prime using the first root; default
                        all).  The choice is part of the *search*, the result
                        is validated by `verify`.
    verify(vectors)  -> exact characteristic-zero check (REQUIRED for any
                        mathematical use; if None the caller must verify).

    Returns a KernelResult; raises Unresolved on budget exhaustion.
    An empty kernel mod p returns immediately with no vectors: that IS a
    proof of trivial kernel in characteristic zero (full column rank).
    """
    m = K.degree
    state: Dict[Tuple[int, ...], dict] = {}
    best_rank = -1
    used: List[int] = []
    for p, roots in split_primes(K, full=True):
        if deadline is not None and time.time() > deadline:
            raise Unresolved("modular kernel: deadline")
        if len(used) >= max_primes:
            raise Unresolved("modular kernel: prime budget exhausted")
        try:
            kernels = [kernel_mod_p(build(p, r), p) for r in roots]
        except BadPrime:
            continue
        piv_sets = {tuple(kr[1]) for kr in kernels}
        if len(piv_sets) != 1:
            continue                       # inconsistent across roots: bad p
        pivots = kernels[0][1]
        free = kernels[0][2]
        if len(pivots) < best_rank:
            continue                       # rank dropped: certainly bad
        if len(pivots) > best_rank:
            best_rank = len(pivots)
            state.clear()                  # all earlier primes were bad
        if not free:
            return KernelResult([], 0, [], [p], [])
        key = tuple(pivots)
        st = state.get(key)
        if st is None:
            sel = list(select(kernels[0][0], p, roots[0])) if select else \
                list(range(len(free)))
            st = state[key] = {"sel": sel, "M": 1, "x": None, "primes": [],
                               "free": free}
        sel = st["sel"]
        vals = [kr[0][sel] for kr in kernels]            # each (k, ncols)
        if m == 1:
            coeff = [vals[0]]
        else:
            W = _vandermonde_inverse(roots, p)
            coeff = []
            for k in range(m):
                acc = np.zeros_like(vals[0])
                for j in range(m):
                    acc = (acc + vals[j] * W[k][j] % p) % p
                coeff.append(acc)
        flat = [int(v) for arr in coeff for v in arr.ravel()]
        if st["x"] is None:
            st["x"], st["M"] = flat, p
        else:
            M = st["M"]
            Minv = pow(M % p, p - 2, p)
            st["x"] = [x + M * (((r - x) * Minv) % p)
                       for x, r in zip(st["x"], flat)]
            st["M"] = M * p
        st["primes"].append(p)
        used.append(p)
        if len(st["primes"]) % check_every and len(st["primes"]) > 1:
            continue
        recon = [rational_reconstruct(x, st["M"]) for x in st["x"]]
        if any(r is None for r in recon):
            continue
        k, ncols = len(sel), vals[0].shape[1]
        vectors = []
        for t in range(k):
            vec = []
            for c in range(ncols):
                cs = [recon[kk * k * ncols + t * ncols + c] for kk in range(m)]
                vec.append(K._make(cs))
            vectors.append(vec)
        if verify is None or verify(vectors):
            return KernelResult(vectors, len(free), free, st["primes"], sel)
    raise Unresolved("modular kernel: no primes")


# --------------------------------------------------------------------------
# small exact (characteristic-zero) elimination -- reference implementation
# --------------------------------------------------------------------------

def exact_rref(rows: Sequence[Sequence], ncols: int):
    A = [list(r) for r in rows]
    pivots = []
    r = 0
    for c in range(ncols):
        piv = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        pv = A[r][c]
        A[r] = [v / pv for v in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [vi - f * vr for vi, vr in zip(A[i], A[r])]
        pivots.append(c)
        r += 1
        if r == len(A):
            break
    return A[:r], pivots


def exact_rank(rows, ncols) -> int:
    return len(exact_rref(rows, ncols)[1])


def exact_nullspace(rows, ncols) -> List[list]:
    R, pivots = exact_rref(rows, ncols)
    pivset = set(pivots)
    out = []
    for fc in (c for c in range(ncols) if c not in pivset):
        v = [Fraction(0)] * ncols
        v[fc] = Fraction(1)
        for i, pc in enumerate(pivots):
            v[pc] = -R[i][fc]
        out.append(v)
    return out
