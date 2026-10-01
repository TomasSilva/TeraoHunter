"""
certificates.py -- exact freeness / non-freeness certificates (search side).

Statuses are explicit; nothing here converts a failure into a mathematical
claim:

  FREE_CERTIFIED       exact Saito certificate found AND re-verified
  NONFREE_CERTIFIED    exact non-freeness certificate found AND re-verified
  UNRESOLVED           budget / timeout / bad primes / unsupported input
  ERROR                unexpected exception (recorded verbatim)

Certificates
------------
saito:
    theta1, theta2 homogeneous of degrees d1 <= d2, d1 + d2 = n - 1,
    alpha_i | theta_j(alpha_i) for every line, det(E, theta1, theta2) = c Q
    with c != 0.  (Saito's criterion.)

jacobian_syzygy:
    (f, g, h) != 0 homogeneous of degree r with f Q_x + g Q_y + h Q_z = 0.
    It proves mdr(Q) <= r -- existence, NOT minimality.  Together with a
    lattice whose reduced characteristic polynomial is (t-d1)(t-d2), d1<=d2,
    and r < d1 it proves non-freeness: a free arrangement with this lattice
    has exponents (1, d1, d2) (Terao factorisation), hence mdr = d1.
    A Jacobian syzygy is never an Euler multiple (E(Q) = nQ != 0), so Euler
    multiples cannot produce false positives.

terao_factorisation:
    the reduced characteristic polynomial t^2 - (n-1) t + (b2 - n + 1) has no
    admissible integer roots.  Lattice-only; such arrangements are useful
    non-free anchors but can never be the second member of a counterexample
    (a free arrangement never shares their lattice).

The independent re-verification lives in `terao_hunter.verifier` and shares no
polynomial or linear-algebra code with this module.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Optional, Tuple

from .arrangement import Arrangement
from .derivations import (apply_theta_to_line, choose_gauge_line, euler,
                          gauge_matrix_mod_p, gauge_vector_to_theta)
from .linalg import (BadPrime, Unresolved, modular_kernel, rank_mod_p,
                     reduce_scalar, split_primes)
from .polys import (dim_S, det3_polys, monomials, p_add, p_degree, p_diff,
                    p_divide_by_linear, p_from_json, p_mul,
                    p_product_of_lines, p_restrict_to_line, p_scale, p_sub,
                    p_to_json)

FREE_CERTIFIED = "free_certified"
NONFREE_CERTIFIED = "nonfree_certified"
UNRESOLVED = "unresolved"
ERROR = "error"

CERT_VERSION = 1


@dataclass
class Classification:
    status: str
    exponents: Optional[Tuple[int, int]] = None      # forced by the lattice
    certificate: Optional[dict] = None
    reason: str = ""
    facts: dict = field(default_factory=dict)        # exact side facts
    hints: dict = field(default_factory=dict)        # NON-proof diagnostics
    seconds: float = 0.0


# --------------------------------------------------------------------------
# modular rank facts
# --------------------------------------------------------------------------

def syzygy_nullity_mod_p(arr: Arrangement, r: int, h: Optional[int] = None,
                         n_primes: int = 1, skip: int = 0):
    """[(p, root, nullity_p)] for D_H(A)_r.  nullity_p is an UPPER bound for
    the characteristic-zero nullity; nullity_p == 0 is a proof that
    D_H(A)_r = 0, i.e. there is no nonzero Jacobian syzygy of degree <= r."""
    if r < 0:
        return []
    h = choose_gauge_line(arr) if h is None else h
    ncols = (r + 1) * (r + 2)
    out = []
    seen = 0
    for p, roots in split_primes(arr.field, full=False):
        seen += 1
        if seen <= skip:
            continue
        try:
            M = gauge_matrix_mod_p(arr, h, r, p, roots[0])
        except BadPrime:
            continue
        out.append((p, roots[0], ncols - rank_mod_p(M, p)))
        if len(out) >= n_primes:
            break
    return out


def prove_no_syzygy_below(arr: Arrangement, d: int, h: Optional[int] = None,
                          max_primes: int = 3) -> dict:
    """Try to PROVE mdr >= d (no nonzero Jacobian syzygy of degree < d).

    proved=True  : full column rank mod a valid prime at degree d-1 (sound).
    proved=False : every tried prime was rank deficient -- a HINT that a
                   syzygy exists, never a proof.
    """
    if d <= 0:
        return {"proved": True, "degree": d - 1, "reason": "vacuous"}
    facts = syzygy_nullity_mod_p(arr, d - 1, h, n_primes=max_primes)
    for (p, root, nul) in facts:
        if nul == 0:
            return {"proved": True, "degree": d - 1, "prime": p, "root": root,
                    "reason": "full_column_rank_mod_p"}
    return {"proved": False, "degree": d - 1,
            "nullities_mod_p": [[p, nul] for (p, _r, nul) in facts],
            "reason": "rank_deficient_mod_all_tried_primes(hint_only)"}


def graded_nullity_hint(arr: Arrangement, degrees, h: Optional[int] = None):
    """{r: nullity_p of D_H(A)_r} for one prime.  Upper bounds only; exact
    only where they are 0."""
    out = {}
    for r in degrees:
        facts = syzygy_nullity_mod_p(arr, r, h, n_primes=1)
        out[int(r)] = facts[0][2] if facts else None
    return out


# --------------------------------------------------------------------------
# in-package exact checks (the independent verifier re-does these in sympy)
# --------------------------------------------------------------------------

def _is_log_along_all(arr: Arrangement, theta) -> bool:
    for l in arr.lines:
        if p_restrict_to_line(apply_theta_to_line(theta, l), l):
            return False
    return True


def check_saito(arr: Arrangement, cert: dict) -> Tuple[bool, str]:
    K = arr.field
    d1, d2 = cert["d1"], cert["d2"]
    if not (isinstance(d1, int) and isinstance(d2, int) and 0 <= d1 <= d2
            and d1 + d2 == arr.n - 1):
        return False, "degree_bookkeeping"
    thetas = []
    for key, d in (("theta1", d1), ("theta2", d2)):
        th = tuple(p_from_json(c, K) for c in cert[key])
        if len(th) != 3 or not any(th):
            return False, f"{key}_zero_or_malformed"
        for comp in th:
            if comp and p_degree(comp) != d:
                return False, f"{key}_wrong_degree"
        if not _is_log_along_all(arr, th):
            return False, f"{key}_not_logarithmic"
        thetas.append(th)
    c = K.elem_from_json(cert["c"])
    if c == 0:
        return False, "c_is_zero"
    det = det3_polys([euler(), thetas[0], thetas[1]])
    Q = p_product_of_lines(arr.lines)
    if p_sub(det, p_scale(Q, c)):
        return False, "determinant_identity_fails"
    return True, "ok"


def check_syzygy(arr: Arrangement, cert: dict) -> Tuple[bool, str]:
    K = arr.field
    comps = tuple(p_from_json(c, K) for c in cert["components"])
    r = cert["degree"]
    if len(comps) != 3 or not any(comps):
        return False, "zero_or_malformed"
    for comp in comps:
        if comp and p_degree(comp) != r:
            return False, "wrong_degree"
    Q = p_product_of_lines(arr.lines)
    total = {}
    for var, comp in enumerate(comps):
        total = p_add(total, p_mul(comp, p_diff(Q, var)))
    if total:
        return False, "syzygy_identity_fails"
    return True, "ok"


# --------------------------------------------------------------------------
# finding certificates
# --------------------------------------------------------------------------

def find_syzygy_certificate(arr: Arrangement, r: int, deadline=None,
                            h: Optional[int] = None, max_primes: int = 400
                            ) -> Optional[dict]:
    """Exact nonzero Jacobian syzygy of degree r, or None when D_H(A)_r = 0
    is PROVED.  Raises Unresolved when neither is established."""
    h = choose_gauge_line(arr) if h is None else h
    K = arr.field

    def build(p, root):
        return gauge_matrix_mod_p(arr, h, r, p, root)

    def make(vec):
        theta = gauge_vector_to_theta(arr, h, r, vec)
        hpoly = {}
        for l in arr.lines:
            val = apply_theta_to_line(theta, l)
            if not val:
                continue
            quo, ok = p_divide_by_linear(val, l)
            if not ok:
                return None                       # not logarithmic: reject
            hpoly = p_add(hpoly, quo)
        from fractions import Fraction
        E = euler()
        scale = Fraction(-1, arr.n)
        return tuple(p_add(theta[i], p_scale(p_mul(hpoly, E[i]), scale))
                     for i in range(3))

    holder = {}

    def verify(vectors):
        sig = make(vectors[0])
        if sig is None or not any(sig):
            return False
        cert = {"type": "jacobian_syzygy", "version": CERT_VERSION,
                "degree": int(r),
                "components": [p_to_json(c, K) for c in sig],
                "arrangement": arr.to_json(), "gauge_line": int(h)}
        ok, _ = check_syzygy(arr, cert)
        if ok:
            holder["cert"] = cert
        return ok

    res = modular_kernel(K, build, select=lambda basis, p, root: [0],
                         verify=verify, max_primes=max_primes,
                         deadline=deadline)
    if not res.vectors:
        return None
    holder["cert"]["primes_used"] = len(res.primes)
    return holder["cert"]


def find_saito_certificate(arr: Arrangement, d1: int, d2: int, deadline=None,
                           h: Optional[int] = None, max_primes: int = 400):
    """(certificate | None, reason).  None + reason starting with 'proved_'
    is an exact negative; anything else without a certificate is a failure
    of the search (raises Unresolved) -- never a proof."""
    if not (1 <= d1 <= d2 and d1 + d2 == arr.n - 1):
        raise ValueError("need 1 <= d1 <= d2, d1 + d2 = n - 1")
    h = choose_gauge_line(arr) if h is None else h
    K = arr.field
    need1 = 2 if d1 == d2 else 1

    def build1(p, root):
        return gauge_matrix_mod_p(arr, h, d1, p, root)

    def verify_log(r):
        def _v(vectors):
            return all(_is_log_along_all(
                arr, gauge_vector_to_theta(arr, h, r, v)) for v in vectors)
        return _v

    res1 = modular_kernel(K, build1, select=lambda b, p, r: list(range(
        min(need1, b.shape[0]))), verify=verify_log(d1),
        max_primes=max_primes, deadline=deadline)
    if len(res1.vectors) < need1:
        # nullity_p < need1 for a valid prime => exact nullity < need1
        return None, "proved_degree_d1_kernel_too_small"
    theta1 = gauge_vector_to_theta(arr, h, d1, res1.vectors[0])
    if d1 == d2:
        theta2 = gauge_vector_to_theta(arr, h, d2, res1.vectors[1])
    else:
        mons2 = monomials(d2)
        N2 = len(mons2)

        def select2(basis, p, root):
            # choose a basis vector with det(E, theta1, v)(pt) != 0 mod p
            import random as _r
            rng = _r.Random(p)
            pt = [rng.randrange(1, p) for _ in range(3)]
            t1 = [_eval_mod(c, pt, p, root) for c in theta1]
            from .derivations import gauge_indices
            e, u, v = gauge_indices(arr, h)
            lh = [reduce_scalar(c, p, root) for c in arr.lines[h]]
            inv = pow(lh[e], p - 2, p)
            mv = [pow(pt[0], a, p) * pow(pt[1], b, p) % p * pow(pt[2], c, p)
                  % p for (a, b, c) in mons2]
            for j in range(basis.shape[0]):
                fu = sum(int(basis[j, k]) * mv[k] for k in range(N2)) % p
                fv = sum(int(basis[j, N2 + k]) * mv[k] for k in range(N2)) % p
                comp = [0, 0, 0]
                comp[u], comp[v] = fu, fv
                comp[e] = (-(lh[u] * fu + lh[v] * fv) * inv) % p
                d = (pt[0] * (t1[1] * comp[2] - t1[2] * comp[1])
                     - pt[1] * (t1[0] * comp[2] - t1[2] * comp[0])
                     + pt[2] * (t1[0] * comp[1] - t1[1] * comp[0])) % p
                if d:
                    return [j]
            return [0]

        res2 = modular_kernel(K, lambda p, root: gauge_matrix_mod_p(
            arr, h, d2, p, root), select=select2, verify=verify_log(d2),
            max_primes=max_primes, deadline=deadline)
        if not res2.vectors:
            return None, "proved_degree_d2_kernel_trivial"
        theta2 = gauge_vector_to_theta(arr, h, d2, res2.vectors[0])
    det = det3_polys([euler(), theta1, theta2])
    Q = p_product_of_lines(arr.lines)
    if not det:
        raise Unresolved("saito: chosen pair has zero determinant")
    m0 = next(iter(Q))
    c = det.get(m0, 0) / Q[m0]
    cert = {"type": "saito", "version": CERT_VERSION, "d1": int(d1),
            "d2": int(d2), "c": K.elem_to_json(c),
            "theta1": [p_to_json(t, K) for t in theta1],
            "theta2": [p_to_json(t, K) for t in theta2],
            "arrangement": arr.to_json(), "gauge_line": int(h)}
    ok, why = check_saito(arr, cert)
    if not ok:
        raise Unresolved(f"saito: candidate pair rejected ({why})")
    return cert, "certified"


def _eval_mod(poly, pt, p, root):
    acc = 0
    for (a, b, c), coef in poly.items():
        acc += reduce_scalar(coef, p, root) * pow(pt[0], a, p) % p * \
            pow(pt[1], b, p) % p * pow(pt[2], c, p)
    return acc % p


# --------------------------------------------------------------------------
# classification
# --------------------------------------------------------------------------

def classify(arr: Arrangement, time_budget: float = 60.0,
             max_primes: int = 400, want_certificate: bool = True
             ) -> Classification:
    """Exact classification with explicit failure modes."""
    t0 = time.time()
    deadline = t0 + time_budget
    try:
        ok, why = arr.admissibility()
        if not ok:
            return Classification(ERROR, reason=f"inadmissible:{why}")
        exps = arr.candidate_exponents()
        if exps is None:
            cert = {"type": "terao_factorisation", "version": CERT_VERSION,
                    "n": arr.n, "b2": arr.b2(),
                    "reduced_constant": arr.reduced_constant(),
                    "arrangement": arr.to_json()}
            return Classification(NONFREE_CERTIFIED, None, cert,
                                  "reduced_char_poly_has_no_integer_roots",
                                  seconds=time.time() - t0)
        d1, d2 = exps
        if d1 < 1:
            return Classification(UNRESOLVED, exps,
                                  reason="degenerate_exponent_zero")
        proof = prove_no_syzygy_below(arr, d1)
        facts = {"no_syzygy_below_d1": proof}
        if proof["proved"]:
            if not want_certificate:
                return Classification(UNRESOLVED, exps, None,
                                      "mdr_ge_d1_proved_certificate_skipped",
                                      facts, seconds=time.time() - t0)
            cert, reason = find_saito_certificate(arr, d1, d2, deadline,
                                                  max_primes=max_primes)
            if cert is not None:
                return Classification(FREE_CERTIFIED, exps, cert, reason,
                                      facts, seconds=time.time() - t0)
            return Classification(UNRESOLVED, exps, None, reason, facts,
                                  seconds=time.time() - t0)
        cert = find_syzygy_certificate(arr, d1 - 1, deadline,
                                       max_primes=max_primes)
        if cert is not None:
            cert["forced_exponents"] = [d1, d2]
            return Classification(NONFREE_CERTIFIED, exps, cert,
                                  "syzygy_below_forced_d1", facts,
                                  seconds=time.time() - t0)
        facts["no_syzygy_below_d1"] = {"proved": True, "degree": d1 - 1,
                                       "reason": "kernel_trivial_mod_later_prime"}
        cert, reason = find_saito_certificate(arr, d1, d2, deadline,
                                              max_primes=max_primes)
        if cert is not None:
            return Classification(FREE_CERTIFIED, exps, cert, reason, facts,
                                  seconds=time.time() - t0)
        return Classification(UNRESOLVED, exps, None, reason, facts,
                              seconds=time.time() - t0)
    except Unresolved as e:
        return Classification(UNRESOLVED, reason=f"unresolved:{e}",
                              seconds=time.time() - t0)
    except Exception as e:                       # noqa: BLE001
        return Classification(ERROR, reason=f"exception:{type(e).__name__}:{e}",
                              seconds=time.time() - t0)


def exact_mdr(arr: Arrangement, time_budget: float = 60.0, r_max=None):
    """Exact minimal degree of a Jacobian syzygy with its two-part proof:
    (i) full column rank mod a valid prime in degree r-1, (ii) an exact
    verified syzygy in degree r.  Raises Unresolved otherwise."""
    deadline = time.time() + time_budget
    r_max = arr.n - 1 if r_max is None else r_max
    for r in range(1, r_max + 1):
        facts = syzygy_nullity_mod_p(arr, r, n_primes=1)
        if facts and facts[0][2] == 0:
            continue
        cert = find_syzygy_certificate(arr, r, deadline)
        if cert is None:
            continue
        lower = prove_no_syzygy_below(arr, r)
        if not lower["proved"]:
            raise Unresolved("mdr: lower bound not proved")
        return r, {"lower_bound": lower, "syzygy": cert}
    raise Unresolved("mdr: no syzygy found up to r_max")
