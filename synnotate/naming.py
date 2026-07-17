# VENDORED VERBATIM from src/cohyp/annotation/naming.py -- do NOT edit here.
# synnotate_pkg must stay importable without the research tree, but the naming rule must
# never fork from the one that BUILT the corpus, or the CLI would label new genomes by a
# different rule than the model was trained under. Re-copy on any change upstream.
"""Single source of truth for RefSeq product-name handling.

Used by BOTH vocabulary construction and gene labelling so the two can never drift apart.
Deterministic, dependency-free, unit-tested (run: python -m doctest src/cohyp/annotation/naming.py -v).

METHOD (paper-ready)
--------------------
A CDS receives exactly one label, by a fixed precedence:

  1. eggNOG orthology  : KO -> Pfam cascade gives a `group_key`; if that key is in the frozen
                         vocabulary map, the gene takes its label.
  2. RefSeq product    : else, if the assembly's RefSeq/PGAP product string `is_informative`,
                         match it onto the vocabulary (exact, then `canonical`); unmatched -> OTHER.
  3. UNKNOWN (dark)    : else.

  => INVARIANT: a gene is dark iff NEITHER eggNOG NOR RefSeq names a function for it.

`is_informative` asks one question: *does this string name a FUNCTION?* A string is uninformative when
it names nothing, or names only a FAMILY/DOMAIN rather than a function. Only markers that EXPLICITLY declare unknown function are
uninformative ("hypothetical protein", "DUF####", "domain-containing protein"). We deliberately do NOT
attempt to strip locus-tag-derived family names ("Rv3235 family protein") -- see the NOTE below: it is
not decidable from the string, and every attempt also destroys real functions ("recombinase family
protein"). Frequency (>=min_count) is the only further gate.

`canonical` exists ONLY to match naming variants of the same function onto one label
("FtsL-like putative cell division protein" -> "cell division protein FtsL"). It is never itself a
label. It is deliberately a token SET, not an edit distance: FtsL/FtsB/FtsQ differ by one character but
are different proteins, so no fuzzy/similarity matching is used anywhere.

There are NO per-family special cases (no mobile-element collapsing, no leader-peptide rules) and ONE
frequency threshold (`--label-min-count`, default 50), applied identically to both label sources.
"""
from __future__ import annotations
import re

VERSION = "1.0"

# (a) strings that explicitly declare unknown function
_EXPLICIT_UNKNOWN = (
    "hypothetical protein", "uncharacterized protein", "uncharacterised protein",
    "protein of unknown function", "domain of unknown function", "putative uncharacterized",
)
# (b) names a domain/motif rather than a function
_DOMAIN_ONLY = re.compile(r"domain-containing|motif-containing", re.I)
# (c) unnamed-family accessions -- the name explicitly declares unknown function
_DUF = re.compile(r"\b(?:DUF|UPF)\s*\d+", re.I)

# NOTE (rejected rule): we do NOT try to strip "<identifier> family protein" names such as
# "Rv3235 family protein" (a locus-tag-derived family). It is tempting, because those name a family
# rather than a function -- but it is not decidable from the string. "IS110" (a curated ISfinder
# family) and "Rv3235" (an arbitrary locus tag) are structurally identical (2 letters + 3-4 digits),
# as are "VOC"/"SRPBCC" (curated superfamilies) vs "YncE"/"HGxxPAAW" (placeholders). A regex for this
# flagged 4,708 strings / 1.75M genes including "recombinase family protein", "transposase family
# protein", "glycosyltransferase family protein" and "porin family protein" -- all real functions.
# Any working version would be a hand-curated prefix list, which is neither complete nor reproducible.
# We therefore treat ONLY explicit unknown-function markers as uninformative (rules a-c), which is the
# rule RefSeq itself makes explicit, and let the >=min_count frequency gate handle the rest.

# tokens dropped when building the canonical form (they never distinguish a function)
_NOISE = frozenset({"protein", "family", "putative", "probable", "predicted", "domain",
                    "containing", "like", "conserved", "uncharacterized", "putative"})
_MULTISPECIES = re.compile(r"^MULTISPECIES:\s*", re.I)
_NONALNUM = re.compile(r"[^a-z0-9]+")


def is_informative(product: str) -> bool:
    """True iff the RefSeq product string names a FUNCTION.

    >>> is_informative("cell division protein FtsL")
    True
    >>> is_informative("FtsL-like putative cell division protein")
    True
    >>> is_informative("IS110 family transposase")          # identifier + real function -> keep
    True
    >>> is_informative("NepR family anti-sigma factor")
    True
    >>> is_informative("hypothetical protein")
    False
    >>> is_informative("DUF2167 domain-containing protein")
    False
    >>> is_informative("recombinase family protein")        # the family name IS the function
    True
    >>> is_informative("VOC family protein")                # curated superfamily
    True
    >>> is_informative("Rv3235 family protein")             # locus-tag family: NOT decidable, kept
    True
    >>> is_informative("DUF6879 family protein")            # explicit unknown-function marker
    False
    >>> is_informative("SCO7613 C-terminal domain-containing membrane protein")
    False
    >>> is_informative("")
    False
    """
    if not product:
        return False
    p = _MULTISPECIES.sub("", str(product)).strip()
    if not p or p.lower() in ("na", "n/a", "none", "-"):
        return False
    pl = p.lower()
    if any(h in pl for h in _EXPLICIT_UNKNOWN):
        return False
    if _DOMAIN_ONLY.search(pl) or _DUF.search(pl):
        return False
    return True


def canonical(product: str) -> str:
    """Canonical form for matching naming variants onto one vocabulary label.

    Lowercase -> drop MULTISPECIES: -> punctuation to space -> drop noise tokens -> sorted unique set.

    >>> canonical("FtsL-like putative cell division protein") == canonical("cell division protein FtsL")
    True
    >>> canonical("cell division protein FtsL") == canonical("cell division protein FtsB")
    False
    >>> canonical("MULTISPECIES: cell division protein FtsL")
    'cell division ftsl'
    """
    s = _MULTISPECIES.sub("", str(product)).lower()
    s = _NONALNUM.sub(" ", s)
    return " ".join(sorted({t for t in s.split() if t and t not in _NOISE}))


def build_canonical_index(labels) -> dict:
    """canonical(label) -> label, dropping canonical forms claimed by >1 distinct label.

    Ambiguous forms are dropped rather than resolved arbitrarily, so matching is deterministic.

    >>> idx = build_canonical_index(["cell division protein FtsL", "cell division protein FtsB"])
    >>> idx[canonical("FtsL-like putative cell division protein")]
    'cell division protein FtsL'
    """
    by = {}
    for L in labels:
        by.setdefault(canonical(L), set()).add(L)
    return {c: next(iter(s)) for c, s in by.items() if len(s) == 1}


def match_label(product: str, labels: set, canon_index: dict):
    """Map a product string onto a vocabulary label. Returns (label, how) or (None, reason).

    >>> labs = {"cell division protein FtsL"}
    >>> idx = build_canonical_index(labs)
    >>> match_label("cell division protein FtsL", labs, idx)
    ('cell division protein FtsL', 'exact')
    >>> match_label("FtsL-like putative cell division protein", labs, idx)
    ('cell division protein FtsL', 'canonical')
    >>> match_label("hypothetical protein", labs, idx)
    (None, 'uninformative')
    >>> match_label("glyoxalase", labs, idx)
    (None, 'unmatched')
    """
    if not is_informative(product):
        return None, "uninformative"
    p = _MULTISPECIES.sub("", str(product)).strip()
    if p in labels:
        return p, "exact"
    c = canonical(p)
    if c in canon_index:
        return canon_index[c], "canonical"
    return None, "unmatched"
