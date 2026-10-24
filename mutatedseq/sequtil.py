"""Small, dependency-free sequence helpers.

The standard genetic code is used.  Ambiguous bases (N, R, Y, ...) translate to
X; this is important because Ensembl pads 5'-partial CDS models with Ns.
"""

from __future__ import annotations

import re
from typing import Iterable, Optional, Tuple

_COMPLEMENT = str.maketrans(
    "ACGTRYMKBDHVNacgtrymkbdhvn",
    "TGCAYRKMVHDBNtgcayrkmvhdbn",
)

# Standard genetic code, DNA codons.  The order of bases in each position is
# T, C, A, G, which matches the classic NCBI table layout.
_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {}
for _i, _b1 in enumerate(_BASES):
    for _j, _b2 in enumerate(_BASES):
        for _k, _b3 in enumerate(_BASES):
            CODON_TABLE[_b1 + _b2 + _b3] = _AAS[_i * 16 + _j * 4 + _k]

STOP_CODONS = frozenset(("TAA", "TAG", "TGA"))
START_CODONS = frozenset(("ATG", "CTG", "TTG", "GTG"))  # alternate starts

_AA_3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "X": "Xaa", "*": "Ter", "U": "Sec", "O": "Pyl",
}


def reverse_complement(seq: str) -> str:
    """Reverse complement DNA, preserving case."""
    return seq.translate(_COMPLEMENT)[::-1]


def translate(seq: str, stop_at_stop: bool = True) -> Tuple[str, bool]:
    """Translate a nucleotide sequence.

    Returns ``(protein, stop_found)``.  Only complete codons are translated;
    a trailing 1-2 bases are represented by ``X`` so that partial models do not
    silently disappear.  ``stop_found`` is True when an in-frame stop codon was
    encountered and translation was cut there.
    """
    seq = seq.upper().replace("U", "T")
    out = []
    stop_found = False
    stopped_early = False
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        aa = CODON_TABLE.get(codon, "X")
        if aa == "*":
            stop_found = True
            if stop_at_stop:
                stopped_early = True
                break
            out.append("*")
        else:
            out.append(aa)
    # Only pad a trailing partial codon when translation actually reached the
    # end.  A CDS followed by a 3' UTR has a non-multiple-of-3 total length;
    # after cutting at the stop, adding X would create a spurious amino acid.
    if len(seq) % 3 and not stopped_early:
        out.append("X")
    return "".join(out), stop_found


def first_stop_codon(seq: str) -> Optional[Tuple[int, str]]:
    """Return ``(nt_offset, codon)`` for the first in-frame stop, else None."""
    seq = seq.upper().replace("U", "T")
    for i in range(0, len(seq) - 2, 3):
        if seq[i : i + 3] in STOP_CODONS:
            return i, seq[i : i + 3]
    return None


def aa3(aa: str) -> str:
    return _AA_3.get(aa, aa)


def _common_prefix_suffix(ref: str, alt: str) -> Tuple[int, int]:
    p = 0
    while p < min(len(ref), len(alt)) and ref[p] == alt[p]:
        p += 1
    s = 0
    while (
        s < min(len(ref), len(alt)) - p
        and ref[len(ref) - 1 - s] == alt[len(alt) - 1 - s]
    ):
        s += 1
    return p, s


def protein_diff(ref: str, alt: str) -> str:
    """Return a compact, human-readable protein difference.

    This is intentionally simpler than full HGVS: it reports the first changed
    residue and, for frameshifts, the number of novel residues before a stop.
    The TSV also contains both complete proteins, so downstream users can do a
    more formal annotation if needed.
    """
    if ref == alt:
        return ""
    p, s = _common_prefix_suffix(ref, alt)
    ref_mid = ref[p : len(ref) - s if s else len(ref)]
    alt_mid = alt[p : len(alt) - s if s else len(alt)]
    pos = p + 1
    if len(ref_mid) == 1 and len(alt_mid) == 1:
        return f"p.{aa3(ref_mid)}{pos}{aa3(alt_mid)}"
    if ref_mid and not alt_mid:
        if len(ref_mid) == 1:
            return f"p.{aa3(ref_mid)}{pos}del"
        return f"p.{aa3(ref_mid[0])}{pos}_{aa3(ref_mid[-1])}{pos + len(ref_mid) - 1}del"
    if alt_mid and not ref_mid:
        if len(alt_mid) == 1:
            return f"p.{pos}_{pos + 1}ins{aa3(alt_mid)}"
        return f"p.{pos}_{pos + 1}ins{len(alt_mid)}aa"
    # Both sides changed: describe the first residue and the net length.
    if ref_mid and alt_mid:
        if len(ref_mid) == len(alt_mid) == 1:
            return f"p.{aa3(ref_mid)}{pos}{aa3(alt_mid)}"
        return (
            f"p.{aa3(ref_mid[0])}{pos}_"
            f"{aa3(ref_mid[-1])}{pos + len(ref_mid) - 1}"
            f"delins{len(alt_mid)}aa"
        )
    return ""


def describe_frameshift(ref: str, alt: str) -> str:
    """Describe the novel C-terminal tail after a frameshift."""
    p, _ = _common_prefix_suffix(ref, alt)
    tail = alt[p:]
    if not tail:
        return ""
    aa = aa3(ref[p]) if p < len(ref) else "Ter"
    novel = len(tail.rstrip("*"))
    stop = tail.endswith("*")
    return f"p.{aa}{p + 1}fs*{novel}" if stop else f"p.{aa}{p + 1}fs*?"


def is_stop_codon(codon: str) -> bool:
    return codon.upper() in STOP_CODONS


def wrap(seq: str, width: int = 60) -> Iterable[str]:
    for i in range(0, len(seq), width):
        yield seq[i : i + width]


_SAFE_ID = re.compile(r"[^A-Za-z0-9_.:-]+")


def safe_id(value: str) -> str:
    """Make a value safe for a FASTA/TSV identifier."""
    return _SAFE_ID.sub("_", str(value))
