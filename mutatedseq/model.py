"""Annotation data structures.

Coordinates in GTF/Segment objects are 1-based, inclusive, exactly as they
appear in Ensembl GTF.  Conversion to Python's 0-based half-open coordinates is
done through ``start0``/``end0``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence


@dataclass
class Segment:
    chrom: str
    start: int
    end: int
    strand: str
    rank: int = 0
    frame: Optional[int] = None
    exon_id: str = ""
    exon_number: Optional[int] = None
    source_feature: str = "CDS"

    @property
    def start0(self) -> int:
        return self.start - 1

    @property
    def end0(self) -> int:
        return self.end

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def contains0(self, pos0: int) -> bool:
        return self.start0 <= pos0 < self.end0

    def overlaps0(self, start0: int, end0: int) -> bool:
        return self.start0 < end0 and start0 < self.end0

    def overlap_len0(self, start0: int, end0: int) -> int:
        return max(0, min(self.end0, end0) - max(self.start0, start0))

    def as_tuple(self) -> tuple:
        return (self.chrom, self.start, self.end, self.strand)


@dataclass
class TranscriptModel:
    transcript_id: str
    gene_id: str
    chrom: str
    strand: str
    start: int
    end: int
    biotype: str = ""
    gene_name: str = ""
    gene_biotype: str = ""
    transcript_name: str = ""
    source: str = ""
    tags: List[str] = field(default_factory=list)
    exons: List[Segment] = field(default_factory=list)
    cds: List[Segment] = field(default_factory=list)
    stops: List[Segment] = field(default_factory=list)
    utr5: List[Segment] = field(default_factory=list)
    utr3: List[Segment] = field(default_factory=list)

    @property
    def start0(self) -> int:
        return self.start - 1

    @property
    def end0(self) -> int:
        return self.end

    @property
    def has_cds(self) -> bool:
        return bool(self.cds)

    @property
    def cds_start(self) -> Optional[int]:
        if not self.cds:
            return None
        return min(s.start for s in self.cds)

    @property
    def cds_end(self) -> Optional[int]:
        if not self.cds:
            return None
        return max(s.end for s in self.cds)

    @property
    def first_phase(self) -> Optional[int]:
        return self.cds[0].frame if self.cds else None

    @property
    def tag_set(self) -> set:
        return set(self.tags)

    @property
    def is_partial_start(self) -> bool:
        tags = self.tag_set
        return "mRNA_start_NF" in tags or "cds_start_NF" in tags

    @property
    def is_partial_end(self) -> bool:
        tags = self.tag_set
        return "mRNA_end_NF" in tags or "cds_end_NF" in tags

    @property
    def is_nmd(self) -> bool:
        return self.biotype == "nonsense_mediated_decay" or "NMD" in self.tag_set

    @property
    def is_nsd(self) -> bool:
        return self.biotype == "non_stop_decay"

    @property
    def is_immune_segment(self) -> bool:
        return self.biotype.endswith("_gene") and self.biotype.startswith(
            ("IG_", "TR_")
        )

    @property
    def location_key(self) -> tuple:
        return (self.chrom, self.start0, self.end0, self.transcript_id)

    def exon_rank(self, seg: Segment) -> int:
        for i, exon in enumerate(self.exons, 1):
            if seg.start >= exon.start and seg.end <= exon.end:
                return i
        return 0

    def to_summary_dict(self) -> Dict[str, object]:
        return {
            "transcript_id": self.transcript_id,
            "gene_id": self.gene_id,
            "gene_name": self.gene_name,
            "chrom": self.chrom,
            "strand": self.strand,
            "start": self.start,
            "end": self.end,
            "biotype": self.biotype,
            "gene_biotype": self.gene_biotype,
            "transcript_name": self.transcript_name,
            "exon_count": len(self.exons),
            "cds_count": len(self.cds),
            "tags": ",".join(self.tags),
            "is_partial_start": self.is_partial_start,
            "is_partial_end": self.is_partial_end,
            "is_nmd": self.is_nmd,
            "is_nsd": self.is_nsd,
            "is_immune_segment": self.is_immune_segment,
        }


def sort_segments(segments: Sequence[Segment], strand: str) -> List[Segment]:
    """Sort genomic segments into transcript (5' -> 3') order."""
    if strand == "-":
        return sorted(segments, key=lambda s: (-s.start, -s.end))
    return sorted(segments, key=lambda s: (s.start, s.end))


def subtract_segments(
    exon: Segment, cds_segments: Sequence[Segment]
) -> List[Segment]:
    """Return exon pieces not covered by any CDS segment (genomic order)."""
    pieces = [(exon.start, exon.end)]
    for cds in cds_segments:
        new_pieces = []
        for start, end in pieces:
            if cds.end < start or cds.start > end:
                new_pieces.append((start, end))
                continue
            if cds.start > start:
                new_pieces.append((start, cds.start - 1))
            if cds.end < end:
                new_pieces.append((cds.end + 1, end))
        pieces = new_pieces
    out = []
    for i, (start, end) in enumerate(sorted(pieces)):
        out.append(
            Segment(
                exon.chrom,
                start,
                end,
                exon.strand,
                rank=i,
                exon_id=exon.exon_id,
                exon_number=exon.exon_number,
                source_feature="UTR",
            )
        )
    return out


def derive_utrs(tx: TranscriptModel) -> None:
    """Fill tx.utr5/tx.utr3 when the GTF does not provide explicit UTR rows."""
    if not tx.cds or not tx.exons:
        return
    cds_start = tx.cds_start
    cds_end = tx.cds_end
    assert cds_start is not None and cds_end is not None
    utr5: List[Segment] = []
    utr3: List[Segment] = []
    coding = list(tx.cds) + list(tx.stops)
    for exon in tx.exons:
        for piece in subtract_segments(exon, coding):
            # Transcript-order position relative to the CDS.
            if tx.strand == "+":
                if piece.end < cds_start:
                    utr5.append(piece)
                elif piece.start > cds_end:
                    utr3.append(piece)
            else:
                if piece.start > cds_end:
                    utr5.append(piece)
                elif piece.end < cds_start:
                    utr3.append(piece)
    tx.utr5 = sort_segments(utr5, tx.strand)
    tx.utr3 = sort_segments(utr3, tx.strand)
