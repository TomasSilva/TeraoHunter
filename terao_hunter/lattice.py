"""
lattice.py -- full intersection lattices and their identification.

A rank-3 lattice is recorded as the FULL labelled list of rank-two flats
(double points included) on line labels 0..n-1.

Identification policy (never weakened anywhere in the package):

  * cheap invariants and the colour-refinement (WL) key are used ONLY for
    bucketing;
  * two lattices are declared EQUAL only with an explicit line permutation
    that maps the full flat list of one onto the full flat list of the other
    (re-verified by `verify_isomorphism`, which shares no search state);
  * they are declared DIFFERENT only by an isomorphism-invariant that
    differs, or by an exhaustive individualisation-refinement search;
  * a budget/timeout gives UNRESOLVED -- neither "equal" nor "different".
"""

from __future__ import annotations

import hashlib
import time
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

EQUAL = "equal"
DIFFERENT = "different"
UNRESOLVED = "unresolved"

Flat = Tuple[int, ...]


class Lattice:
    def __init__(self, n: int, flats: Sequence[Sequence[int]], check: bool = True):
        self.n = int(n)
        self.flats: Tuple[Flat, ...] = tuple(sorted(tuple(sorted(f))
                                                    for f in flats))
        if check:
            self._check()
        self._wl: Optional[Tuple[str, List[str]]] = None

    def _check(self):
        """Every pair of lines lies in exactly one flat (rank-3 geometry)."""
        seen = set()
        for f in self.flats:
            if len(f) < 2 or len(set(f)) != len(f):
                raise ValueError(f"malformed flat {f}")
            for a in range(len(f)):
                for b in range(a + 1, len(f)):
                    pr = (f[a], f[b])
                    if pr in seen:
                        raise ValueError(f"pair {pr} lies in two flats")
                    seen.add(pr)
        if len(seen) != self.n * (self.n - 1) // 2:
            raise ValueError("flat list does not cover every pair of lines")

    @staticmethod
    def from_arrangement(arr) -> "Lattice":
        return Lattice(arr.n, arr.flats(), check=False)

    # -- invariants (bucketing only) ---------------------------------------
    def multiple_flats(self) -> List[Flat]:
        return [f for f in self.flats if len(f) >= 3]

    def t_vector(self) -> Dict[int, int]:
        t: Dict[int, int] = {}
        for f in self.flats:
            t[len(f)] = t.get(len(f), 0) + 1
        return dict(sorted(t.items()))

    def b2(self) -> int:
        return sum(len(f) - 1 for f in self.flats)

    def m_max(self) -> int:
        return max(len(f) for f in self.flats)

    def line_point_counts(self) -> List[int]:
        c = [0] * self.n
        for f in self.flats:
            for i in f:
                c[i] += 1
        return c

    def line_profiles(self) -> List[Tuple[int, ...]]:
        prof: List[List[int]] = [[] for _ in range(self.n)]
        for f in self.flats:
            for i in f:
                prof[i].append(len(f))
        return [tuple(sorted(p)) for p in prof]

    def coarse_invariant(self) -> Tuple:
        return (self.n, tuple(self.t_vector().items()),
                tuple(sorted(self.line_profiles())))

    # -- colour refinement --------------------------------------------------
    def _graph(self):
        """Bipartite Levi graph: lines 0..n-1, multiple points n.. ."""
        mf = self.multiple_flats()
        adj: List[List[int]] = [[] for _ in range(self.n + len(mf))]
        for k, f in enumerate(mf):
            for i in f:
                adj[i].append(self.n + k)
                adj[self.n + k].append(i)
        return adj

    def wl_colours(self, rounds: Optional[int] = None) -> List[str]:
        return self._wl_data()[1]

    def wl_key(self) -> str:
        """Isomorphism-INVARIANT bucket key (collisions possible)."""
        return self._wl_data()[0]

    def _wl_data(self):
        if self._wl is None:
            adj = self._graph()
            col = _refine(adj, _initial_colours(self, adj))
            hist = sorted(col)
            key = hashlib.sha256(
                (str(self.n) + "|" + "|".join(hist)).encode()).hexdigest()[:32]
            self._wl = (key, col)
        return self._wl

    # -- serialisation -------------------------------------------------------
    def to_json(self) -> dict:
        return {"n": self.n, "flats": [list(f) for f in self.flats]}

    @staticmethod
    def from_json(d: dict) -> "Lattice":
        return Lattice(d["n"], d["flats"])

    def __eq__(self, other):
        """LABELLED equality (same flats on the same labels)."""
        return (isinstance(other, Lattice) and self.n == other.n
                and self.flats == other.flats)

    def __hash__(self):
        return hash((self.n, self.flats))

    def delete(self, i: int) -> "Lattice":
        relabel = {j: (j if j < i else j - 1) for j in range(self.n) if j != i}
        flats = []
        for f in self.flats:
            g = tuple(relabel[j] for j in f if j != i)
            if len(g) >= 2:
                flats.append(g)
        return Lattice(self.n - 1, flats)


def _h(s: str) -> str:
    return hashlib.sha1(s.encode()).hexdigest()[:16]


def _initial_colours(lat: Lattice, adj) -> List[str]:
    col = []
    for v in range(len(adj)):
        kind = "L" if v < lat.n else "P"
        col.append(_h(f"{kind}{len(adj[v])}"))
    return col


def _refine(adj, col: List[str], max_rounds: int = 64) -> List[str]:
    """Colour refinement with isomorphism-invariant colour NAMES (hashes of
    the refinement history), so colours are comparable across graphs."""
    ncls = len(set(col))
    for _ in range(max_rounds):
        new = [_h(col[v] + ":" + ",".join(sorted(col[w] for w in adj[v])))
               for v in range(len(adj))]
        k = len(set(new))
        col = new
        if k == ncls:
            break
        ncls = k
    return col


# --------------------------------------------------------------------------
# explicit isomorphism
# --------------------------------------------------------------------------

@dataclass
class IsoResult:
    status: str                       # EQUAL / DIFFERENT / UNRESOLVED
    perm: Optional[List[int]] = None  # perm[i] = label in B of line i of A
    reason: str = ""
    nodes: int = 0


def verify_isomorphism(A: Lattice, B: Lattice, perm: Sequence[int]) -> bool:
    """Independent check: perm is a bijection of labels carrying the FULL
    flat list of A onto the FULL flat list of B."""
    if A.n != B.n or sorted(perm) != list(range(A.n)):
        return False
    image = sorted(tuple(sorted(perm[i] for i in f)) for f in A.flats)
    return image == list(B.flats)


def find_isomorphism(A: Lattice, B: Lattice, node_budget: int = 20000,
                     time_budget: Optional[float] = 20.0) -> IsoResult:
    if A.coarse_invariant() != B.coarse_invariant():
        return IsoResult(DIFFERENT, reason="coarse_invariant")
    if A.wl_key() != B.wl_key():
        return IsoResult(DIFFERENT, reason="colour_refinement")
    if A.flats == B.flats:
        return IsoResult(EQUAL, list(range(A.n)), "identical_labelling")
    adjA, adjB = A._graph(), B._graph()
    deadline = None if time_budget is None else time.time() + time_budget
    counter = [0]

    class _Budget(Exception):
        pass

    def search(colA: List[str], colB: List[str]) -> Optional[List[int]]:
        counter[0] += 1
        if counter[0] > node_budget or (deadline and time.time() > deadline):
            raise _Budget
        cellsA: Dict[str, List[int]] = {}
        for v, c in enumerate(colA):
            cellsA.setdefault(c, []).append(v)
        cellsB: Dict[str, List[int]] = {}
        for v, c in enumerate(colB):
            cellsB.setdefault(c, []).append(v)
        if {c: len(v) for c, v in cellsA.items()} != \
                {c: len(v) for c, v in cellsB.items()}:
            return None
        target = None
        for c, vs in cellsA.items():
            if len(vs) > 1 and vs[0] < A.n:           # individualise LINES
                if target is None or len(vs) < len(cellsA[target]):
                    target = c
        if target is None:
            perm = [cellsB[colA[i]][0] for i in range(A.n)]
            return perm if verify_isomorphism(A, B, perm) else None
        v = cellsA[target][0]
        for w in cellsB[target]:
            tag = _h(f"ind{counter[0]}:{target}")
            nA = list(colA)
            nB = list(colB)
            nA[v] = tag
            nB[w] = tag
            res = search(_refine(adjA, nA), _refine(adjB, nB))
            if res is not None:
                return res
        return None

    try:
        perm = search(A.wl_colours(), B.wl_colours())
    except _Budget:
        return IsoResult(UNRESOLVED, reason="budget_exhausted",
                         nodes=counter[0])
    if perm is None:
        return IsoResult(DIFFERENT, reason="exhaustive_search",
                         nodes=counter[0])
    return IsoResult(EQUAL, perm, "explicit_permutation", counter[0])


class LatticeRegistry:
    """Buckets by WL key; classes inside a bucket are separated by explicit
    isomorphism only.  `classify` never merges on a hash."""

    def __init__(self, node_budget: int = 20000, time_budget: float = 20.0):
        self.buckets: Dict[str, List[dict]] = {}
        self.node_budget = node_budget
        self.time_budget = time_budget

    def classify(self, lat: Lattice, create: bool = True):
        """Returns (class_id | None, IsoResult | None, is_new).

        class_id is "<wlkey>.<index>".  If some comparison inside the bucket
        is UNRESOLVED and none is EQUAL, a new class flagged
        `identity_unresolved_against` is created (never silently merged,
        never silently declared new).
        """
        key = lat.wl_key()
        bucket = self.buckets.setdefault(key, [])
        unresolved = []
        for rec in bucket:
            res = find_isomorphism(lat, rec["lattice"], self.node_budget,
                                   self.time_budget)
            if res.status == EQUAL:
                return rec["id"], res, False
            if res.status == UNRESOLVED:
                unresolved.append(rec["id"])
        if not create:
            return None, None, False
        cid = f"{key}.{len(bucket)}"
        bucket.append({"id": cid, "lattice": lat,
                       "identity_unresolved_against": unresolved})
        return cid, IsoResult(EQUAL, list(range(lat.n)), "representative"), True

    def get(self, cid: str) -> Optional[dict]:
        key = cid.split(".")[0]
        for rec in self.buckets.get(key, []):
            if rec["id"] == cid:
                return rec
        return None

    def __len__(self):
        return sum(len(b) for b in self.buckets.values())

    def to_json(self):
        return [{"id": r["id"], "lattice": r["lattice"].to_json(),
                 "identity_unresolved_against": r["identity_unresolved_against"]}
                for b in self.buckets.values() for r in b]

    @staticmethod
    def from_json(data, **kw) -> "LatticeRegistry":
        reg = LatticeRegistry(**kw)
        for r in data:
            key = r["id"].split(".")[0]
            reg.buckets.setdefault(key, []).append(
                {"id": r["id"], "lattice": Lattice.from_json(r["lattice"]),
                 "identity_unresolved_against":
                     r.get("identity_unresolved_against", [])})
        return reg
