"""Variant normalization, allele application and mutated-genome blocks.

The central idea is to keep every edit in forward-strand genomic coordinates
(0-based, half-open) and to build a list of reference/alternate blocks.  A
transcript is then reconstructed by intersecting its exon/CDS intervals with
those blocks; this automatically handles deletions, insertions, replacements,
negative-strand genes, and variants that span exon junctions.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .fasta import FastaReader
from .sequtil import reverse_complement
from .vcf import ChosenAllele, VCFRecord


@dataclass
class RefCheck:
    chrom: str
    pos: int
    ref: str
    genome_ref: str
    match: bool
    offset: Optional[int] = None
    offset_strand: str = ""
    context: str = ""
    note: str = ""

    def as_dict(self) -> Dict[str, object]:
        return {
            "ref_match": self.match,
            "genome_ref": self.genome_ref,
            "ref_offset": "" if self.offset is None else self.offset,
            "ref_offset_strand": self.offset_strand,
            "genome_context": self.context,
            "ref_note": self.note,
        }


@dataclass
class Edit:
    chrom: str
    start0: int
    end0: int
    alt: str
    ref: str
    record: VCFRecord
    allele_index: int
    chosen: ChosenAllele
    kind: str = "complex"
    symbolic: bool = False
    unresolved: bool = False
    skipped: bool = False
    note: str = ""
    affected_ref: str = ""
    ref_check: Optional[RefCheck] = None

    @property
    def pos(self) -> int:
        return self.start0 + 1

    @property
    def length_delta(self) -> int:
        return len(self.alt) - (self.end0 - self.start0)

    @property
    def is_insertion(self) -> bool:
        return self.start0 == self.end0 and bool(self.alt)

    @property
    def is_deletion(self) -> bool:
        return self.end0 > self.start0 and self.alt == ""

    @property
    def span(self) -> int:
        return max(1, self.end0 - self.start0)

    @property
    def key(self) -> str:
        return (
            f"{self.chrom}:{self.record.pos}:{self.record.ref}:"
            f"{self.record.alt_string}#{self.allele_index}"
        )

    @property
    def changed_bases(self) -> int:
        """Approximate number of bases changed, after trimming shared anchors."""
        ref = self.affected_ref or self.ref
        alt = self.alt
        p = 0
        while p < min(len(ref), len(alt)) and ref[p] == alt[p]:
            p += 1
        s = 0
        while (
            s < min(len(ref), len(alt)) - p
            and ref[len(ref) - 1 - s] == alt[len(alt) - 1 - s]
        ):
            s += 1
        return max(len(ref) - p - s, len(alt) - p - s)

    def minimal_change(self) -> Tuple[int, str, str]:
        """Return ``(start0, ref, alt)`` after removing shared prefix/suffix."""
        ref = self.affected_ref or self.ref
        alt = self.alt
        start = self.start0
        p = 0
        while p < min(len(ref), len(alt)) and ref[p] == alt[p]:
            p += 1
        ref = ref[p:]
        alt = alt[p:]
        start += p
        s = 0
        while s < min(len(ref), len(alt)) and ref[-1 - s] == alt[-1 - s]:
            s += 1
        if s:
            ref = ref[:-s]
            alt = alt[:-s]
        return start, ref, alt

    def as_dict(self) -> Dict[str, object]:
        mstart, mref, malt = self.minimal_change()
        out: Dict[str, object] = {
            "chrom": self.chrom,
            "pos": self.record.pos,
            "id": self.record.id,
            "ref": self.record.ref,
            "alt_all": self.record.alt_string,
            "alt_selected": self.chosen.alt,
            "allele_index": self.allele_index,
            "qual": "" if self.record.qual is None else self.record.qual,
            "filter": self.record.filter,
            "gt": self.chosen.gt,
            "phased": self.chosen.phased,
            "svtype": self.record.info_value("SVTYPE") or "",
            "end": self.record.info_value("END") or "",
            "svlen": self.record.info_value("SVLEN") or "",
            "kind": self.kind,
            "symbolic": self.symbolic,
            "unresolved": self.unresolved,
            "skipped": self.skipped,
            "edit_start0": self.start0,
            "edit_end0": self.end0,
            "edit_alt": self.alt,
            "minimal_start0": mstart,
            "minimal_ref": mref,
            "minimal_alt": malt,
            "note": self.note,
        }
        if self.ref_check is not None:
            out.update(self.ref_check.as_dict())
        return out


def ref_matches(actual: str, expected: str) -> bool:
    if len(actual) != len(expected):
        return False
    for a, b in zip(actual.upper(), expected.upper()):
        if b != "N" and a != b:
            return False
    return True


def check_ref(
    genome: FastaReader,
    chrom: str,
    pos: int,
    ref: str,
    window: int = 50,
    context_flank: int = 10,
) -> RefCheck:
    """Validate VCF REF against the genome and quantify any offset.

    The offset is reported in bases relative to the VCF POS (negative means the
    REF was found upstream).  Both forward and reverse-complement occurrences
    are searched so that a strand/annotation mistake is visible in the report.
    """
    ref = ref.upper()
    if not genome.has_chrom(chrom):
        return RefCheck(
            chrom, pos, ref, "", False,
            note=f"chromosome {chrom!r} not present in genome FASTA",
        )
    try:
        chrom_len = genome.length(chrom)
    except KeyError:
        chrom_len = 0
    start0 = max(0, pos - 1)
    actual = genome.fetch(chrom, start0, start0 + len(ref))
    left = max(0, start0 - window)
    right = min(chrom_len, start0 + len(ref) + window)
    seq = genome.fetch(chrom, left, right)
    ctx_left = max(0, start0 - context_flank)
    ctx_right = min(chrom_len, start0 + len(ref) + context_flank)
    context = genome.fetch(chrom, ctx_left, ctx_right)
    if ref_matches(actual, ref):
        return RefCheck(chrom, pos, ref, actual, True, 0, "+", context)

    offset: Optional[int] = None
    strand = ""
    if len(ref) >= 2:
        # Searching for a single base in a 100 bp window is meaningless: every
        # base will occur nearby by chance and produce a spurious offset.
        idx = seq.find(ref)
        if idx >= 0:
            offset = (left + idx) - start0
            strand = "+"
        else:
            rc = reverse_complement(ref)
            idx = seq.find(rc)
            if idx >= 0:
                offset = (left + idx) - start0
                strand = "-"
    note = ""
    if offset is not None:
        note = (
            f"REF not found at POS; found {offset:+d} bp "
            f"({strand} strand). Check assembly/contig, 0/1-based offsets, "
            f"or transcript-vs-genome coordinates."
        )
    elif len(ref) == 1:
        note = (
            "single-base REF mismatch; an offset search is not meaningful. "
            "Check the assembly, chromosome naming, or the base at POS."
        )
    else:
        note = (
            "REF not found at POS or in the +/- window. Check the genome "
            "assembly, chromosome naming, or whether the VCF is transcript-based."
        )
    return RefCheck(chrom, pos, ref, actual, False, offset, strand, context, note)


def _sv_end(record: VCFRecord) -> Optional[int]:
    end = record.info_int("END")
    if end is not None:
        return end
    svlen = record.info_int("SVLEN")
    if svlen is not None:
        return record.pos + abs(svlen)
    return None


def _sv_sequence(record: VCFRecord, key: str) -> str:
    value = record.info_value(key)
    if not value:
        return ""
    value = value.split(",")[0].strip()
    if value.lower() in ("nan", "none", "."):
        return ""
    return value.upper()


def make_edit(
    record: VCFRecord,
    chosen: ChosenAllele,
    genome: FastaReader,
    ref_window: int = 50,
    ref_check: Optional[RefCheck] = None,
) -> Edit:
    """Convert one selected VCF allele into a genomic Edit."""
    ref = record.ref.upper()
    alt = chosen.alt.upper()
    start0 = record.pos - 1
    end0 = start0 + len(ref)
    symbolic = alt.startswith("<") or "[" in alt or "]" in alt
    if ref_check is None:
        ref_check = check_ref(
            genome, record.chrom, record.pos, ref, window=ref_window
        )
    if not ref_check.match:
        return Edit(
            record.chrom, start0, end0, alt, ref, record, chosen.index, chosen,
            kind="ref_mismatch", symbolic=symbolic, unresolved=True,
            note=ref_check.note, affected_ref=ref, ref_check=ref_check,
        )
    if symbolic:
        return _make_symbolic_edit(record, chosen, genome, ref_check)
    if alt in ("*", "<*>", "<NON_REF>"):
        return Edit(
            record.chrom, start0, end0, "", ref, record, chosen.index, chosen,
            kind="spanning_deletion", unresolved=True, note="spanning/unsupported allele",
            affected_ref=ref, ref_check=ref_check,
        )
    kind = "complex"
    if len(ref) == len(alt):
        kind = "snv" if len(ref) == 1 else "mnv"
    elif len(ref) < len(alt):
        kind = "ins"
    elif len(ref) > len(alt):
        kind = "del"
    return Edit(
        record.chrom, start0, end0, alt, ref, record, chosen.index, chosen,
        kind=kind, symbolic=False, affected_ref=ref, ref_check=ref_check,
    )


def _make_symbolic_edit(
    record: VCFRecord,
    chosen: ChosenAllele,
    genome: FastaReader,
    ref_check: RefCheck,
) -> Edit:
    alt = chosen.alt.upper()
    svtype = (record.info_value("SVTYPE") or alt.strip("<>").split(":")[0]).upper()
    start0 = record.pos - 1
    end = _sv_end(record)
    anchor = record.ref.upper()
    base_note = ""
    if not genome.has_chrom(record.chrom):
        return Edit(
            record.chrom, start0, start0 + len(anchor), "", record.ref, record,
            chosen.index, chosen, kind=f"symbolic_{svtype.lower()}", symbolic=True,
            unresolved=True, note="chromosome not in genome", affected_ref=anchor,
            ref_check=ref_check,
        )
    chrom_len = genome.length(record.chrom)

    if svtype == "DEL":
        if end is None or end <= record.pos:
            return Edit(
                record.chrom, start0, start0 + len(anchor), "", record.ref, record,
                chosen.index, chosen, kind="symbolic_del", symbolic=True,
                unresolved=True, note="<DEL> without usable END/SVLEN",
                affected_ref=anchor, ref_check=ref_check,
            )
        end = min(end, chrom_len)
        deleted = genome.fetch(record.chrom, record.pos, end)
        return Edit(
            record.chrom, record.pos, end, "", record.ref, record,
            chosen.index, chosen, kind="symbolic_del", symbolic=True,
            affected_ref=deleted, note=base_note, ref_check=ref_check,
        )

    if svtype == "INS":
        seq = _sv_sequence(record, "LEFT_SVINSSEQ")
        source = "LEFT_SVINSSEQ"
        if not seq:
            seq = _sv_sequence(record, "RIGHT_SVINSSEQ")
            source = "RIGHT_SVINSSEQ"
        if not seq:
            return Edit(
                record.chrom, start0, start0 + len(anchor), "", record.ref, record,
                chosen.index, chosen, kind="symbolic_ins", symbolic=True,
                unresolved=True,
                note="<INS> without LEFT_SVINSSEQ/RIGHT_SVINSSEQ; cannot reconstruct",
                affected_ref=anchor, ref_check=ref_check,
            )
        svlen = record.info_int("SVLEN")
        note = f"inserted sequence from {source}"
        if svlen is not None and len(seq) != abs(svlen):
            note += f"; SVLEN={svlen}, sequence length={len(seq)}"
        return Edit(
            record.chrom, record.pos, record.pos, seq, record.ref, record,
            chosen.index, chosen, kind="symbolic_ins", symbolic=True,
            affected_ref="", note=note, ref_check=ref_check,
        )

    if svtype in ("DUP", "DUP:TANDEM", "TANDEM:DUP"):
        if end is None or end <= record.pos:
            return Edit(
                record.chrom, start0, start0 + len(anchor), "", record.ref, record,
                chosen.index, chosen, kind="symbolic_dup", symbolic=True,
                unresolved=True, note="<DUP> without usable END/SVLEN",
                affected_ref=anchor, ref_check=ref_check,
            )
        end = min(end, chrom_len)
        dup_seq = genome.fetch(record.chrom, record.pos, end)
        return Edit(
            record.chrom, end, end, dup_seq, record.ref, record,
            chosen.index, chosen, kind="symbolic_dup", symbolic=True,
            affected_ref="", note="tandem duplication applied after END",
            ref_check=ref_check,
        )

    if svtype in ("INV", "INV3", "INV5"):
        if end is None or end <= record.pos:
            return Edit(
                record.chrom, start0, start0 + len(anchor), "", record.ref, record,
                chosen.index, chosen, kind="symbolic_inv", symbolic=True,
                unresolved=True, note="<INV> without usable END/SVLEN",
                affected_ref=anchor, ref_check=ref_check,
            )
        end = min(end, chrom_len)
        original = genome.fetch(record.chrom, record.pos, end)
        return Edit(
            record.chrom, record.pos, end, reverse_complement(original),
            record.ref, record, chosen.index, chosen, kind="symbolic_inv",
            symbolic=True, affected_ref=original,
            note="inverted interval POS+1..END", ref_check=ref_check,
        )

    return Edit(
        record.chrom, start0, start0 + len(anchor), "", record.ref, record,
        chosen.index, chosen, kind=f"symbolic_{svtype.lower()}", symbolic=True,
        unresolved=True,
        note=f"symbolic <{svtype}> cannot be represented as a single linear sequence",
        affected_ref=anchor, ref_check=ref_check,
    )


def edits_conflict(a: Edit, b: Edit) -> bool:
    """True when two edits cannot be applied to one consensus sequence."""
    if a.start0 == a.end0 and b.start0 == b.end0:
        # Multiple insertions at one point are allowed; order is VCF order.
        return False
    if a.start0 == a.end0:
        return b.start0 < a.start0 < b.end0
    if b.start0 == b.end0:
        return a.start0 < b.start0 < a.end0
    return a.start0 < b.end0 and b.start0 < a.end0


def resolve_edits(
    edits: Sequence[Edit], policy: str = "skip"
) -> Tuple[List[Edit], List[Edit]]:
    """Sort edits and remove/fail on overlaps according to policy."""
    ordered = sorted(
        edits,
        key=lambda e: (
            e.start0,
            e.end0,
            e.record.pos,
            e.allele_index,
        ),
    )
    accepted: List[Edit] = []
    rejected: List[Edit] = []
    # Only edits that can still overlap the current one need to be inspected.
    # This turns the common case from O(n^2) into O(n * local_overlap).
    active: List[Edit] = []
    for edit in ordered:
        if edit.start0 == edit.end0:
            active = [
                a
                for a in active
                if a.end0 > edit.start0
                or (a.start0 == a.end0 == edit.start0)
            ]
        else:
            active = [a for a in active if a.end0 > edit.start0]
        conflict = None
        for prev in active:
            if edits_conflict(prev, edit):
                conflict = prev
                break
        if conflict is not None:
            edit.skipped = True
            edit.note = (
                (edit.note + "; " if edit.note else "")
                + f"overlaps {conflict.key}; not applied"
            )
            if policy == "error":
                raise ValueError(edit.note)
            rejected.append(edit)
            continue
        accepted.append(edit)
        active.append(edit)
    return accepted, rejected


@dataclass
class Block:
    ref_start0: int
    ref_end0: int
    alt: Optional[str]  # None = reference sequence

    @property
    def is_insertion(self) -> bool:
        return self.ref_start0 == self.ref_end0


class MutatedRegion:
    """A reference interval with a set of edits applied, exposed as blocks."""

    def __init__(
        self,
        genome: FastaReader,
        chrom: str,
        edits: Sequence[Edit],
        start0: int,
        end0: int,
    ):
        self.genome = genome
        self.chrom = chrom
        self.edits = sorted(edits, key=lambda e: (e.start0, e.end0))
        if self.edits:
            start0 = min([start0] + [e.start0 for e in self.edits])
            end0 = max([end0] + [e.end0 for e in self.edits])
        self.start0 = max(0, start0)
        self.end0 = max(self.start0, end0)
        self.blocks: List[Block] = []
        cursor = self.start0
        for edit in self.edits:
            if edit.start0 > cursor:
                self.blocks.append(Block(cursor, edit.start0, None))
            if edit.start0 < cursor:
                # resolve_edits should have removed these; be defensive.
                raise ValueError(f"overlapping edit at {edit.chrom}:{edit.pos}")
            self.blocks.append(Block(edit.start0, edit.end0, edit.alt))
            cursor = edit.end0
        if cursor < self.end0:
            self.blocks.append(Block(cursor, self.end0, None))
        if not self.blocks:
            self.blocks.append(Block(self.start0, self.end0, None))

    def _slice_block(self, block: Block, a: int, b: int) -> str:
        inter_start = max(a, block.ref_start0)
        inter_end = min(b, block.ref_end0)
        if inter_end <= inter_start:
            return ""
        if block.alt is None:
            return self.genome.fetch(self.chrom, inter_start, inter_end)
        if inter_start == block.ref_start0 and inter_end == block.ref_end0:
            # The whole replacement is inside the requested interval.  Return
            # the entire ALT even when len(ALT) != len(REF); otherwise an
            # insertion written as REF=T/ALT=TC would lose the inserted base.
            return block.alt
        offset = inter_start - block.ref_start0
        return block.alt[offset : offset + (inter_end - inter_start)]

    def fetch(self, start0: int, end0: int) -> str:
        """Return the mutated plus-strand sequence for a reference interval."""
        if end0 <= start0:
            return ""
        out: List[str] = []
        for block in self.blocks:
            if block.is_insertion:
                if start0 <= block.ref_start0 <= end0 and block.alt:
                    out.append(block.alt)
                continue
            if block.ref_end0 <= start0 or block.ref_start0 >= end0:
                continue
            out.append(self._slice_block(block, start0, end0))
        return "".join(out)

    def extract_segment(self, start0: int, end0: int, strand: str) -> str:
        """Return a transcript-oriented sequence for one exon/CDS interval."""
        seq = self.fetch(start0, end0)
        return reverse_complement(seq) if strand == "-" else seq

    def boundary_insertions(self, start0: int, end0: int) -> int:
        """Count insertions anchored exactly at the interval boundaries."""
        n = 0
        for block in self.blocks:
            if block.is_insertion and block.ref_start0 in (start0, end0):
                n += 1
        return n


@dataclass
class EditIndex:
    """Edits grouped by chromosome and sorted for interval queries."""

    edits_by_chrom: Dict[str, List[Edit]] = field(default_factory=dict)

    def add(self, edit: Edit) -> None:
        self.edits_by_chrom.setdefault(edit.chrom, []).append(edit)

    def finalize(self) -> None:
        for edits in self.edits_by_chrom.values():
            edits.sort(key=lambda e: (e.start0, e.end0, e.record.pos))

    def query(self, chrom: str, start0: int, end0: int) -> List[Edit]:
        edits = self.edits_by_chrom.get(chrom, [])
        if not edits:
            return []
        # Linear scan is fine for transcript-sized intervals; the list is sorted
        # and most chromosomes have far fewer edits than the genome-wide total.
        out = []
        for edit in edits:
            if edit.start0 > end0:
                break
            if edit.end0 < start0:
                continue
            if edit.start0 == edit.end0:
                if start0 <= edit.start0 <= end0:
                    out.append(edit)
            elif edit.end0 > start0 and edit.start0 < end0:
                out.append(edit)
        return out
