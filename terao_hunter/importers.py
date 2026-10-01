"""
importers.py -- read-only adapters for existing project data.

* FreeLineArrangements `certified_free_arrangements.jsonl.gz` records
  (string coordinates; quadratic elements as bracket tokens "[a+bs]" with
  s = sqrt(d)).
* The stored Saito certificates of that database are NOT trusted: anchors
  imported from it are re-certified by this package before use.

Nothing here writes to the source files.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
from fractions import Fraction
from pathlib import Path
from typing import Iterator, Optional

from .arrangement import Arrangement
from .fields import NumberField, QQ, quadratic_field

DEFAULT_DB = Path("/Users/tomassilva/Desktop/GitHub/FreeLineArrangements/"
                  "certified_free_arrangements.jsonl.gz")
DB_SHA256 = "d666d08f3a48fa4278a734c191d977a5d6bc47c09120d2e616936e6134ed2a79"

_TOKEN = re.compile(r"^(?P<a>[+-]?\d+(?:/\d+)?(?=[+-]))?(?P<b>[+-]?(?:\d+(?:/\d+)?)?)s$")


def parse_legacy_scalar(tok, K: NumberField):
    """'3/2', '[1/2-1/2s]', '-[...]', '[-1s]', '[s]' -> exact scalar."""
    if isinstance(tok, int):
        return Fraction(tok)
    t = tok.strip()
    neg = False
    if t.startswith("-["):
        neg, t = True, t[1:]
    if t.startswith("[") and t.endswith("]"):
        if K.is_rational:
            raise ValueError(f"quadratic token {tok!r} in a rational record")
        m = _TOKEN.match(t[1:-1].replace(" ", ""))
        if not m:
            raise ValueError(f"cannot parse quadratic token {tok!r}")
        a = Fraction(m.group("a")) if m.group("a") else Fraction(0)
        b = m.group("b")
        b = Fraction(1) if b in ("", "+") else Fraction(-1) if b == "-" \
            else Fraction(b)
        v = K(a, b)
        return -v if neg else v
    return Fraction(t)


def legacy_field(tag) -> NumberField:
    if tag in (None, "QQ"):
        return QQ
    if isinstance(tag, dict) and tag.get("type") == "quadratic":
        return quadratic_field(int(tag["d"]))
    return NumberField.from_json(tag)


def arrangement_from_legacy(rec: dict) -> Arrangement:
    K = legacy_field(rec.get("field", rec.get("coefficient_field", "QQ")))
    return Arrangement([[parse_legacy_scalar(c, K) for c in l]
                        for l in rec["lines"]], K)


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def iter_database(path: Optional[Path] = None, balanced_only: bool = False,
                  n_min: int = 0, n_max: int = 10 ** 9) -> Iterator[dict]:
    """Yield raw records (the `_meta` line is skipped).  Read-only."""
    path = Path(path or DEFAULT_DB)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            rec = json.loads(line)
            if "_meta" in rec:
                continue
            if not (n_min <= rec["n"] <= n_max):
                continue
            if balanced_only and rec["exponents"][1] != rec["exponents"][2]:
                continue
            yield rec
