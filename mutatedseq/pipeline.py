"""Reference reconstruction, mutation application and effect annotation."""

from __future__ import annotations

import bisect
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from .fasta import FastaReader, read_fasta
from .gtf import Annotation
from .model import Segment, TranscriptModel, derive_utrs, sort_segments
from .sequtil import (
    STOP_CODONS,
    describe_frameshift,
    is_stop_codon,
    protein_diff,
    reverse_complement,
    translate,
)
from .variants import Edit, EditIndex, MutatedRegion


@dataclass
class RunOptions:
    extend_mode: str = "utr"  # none | utr | genomic | utr+genomic
    extend_bp: int = 5000
    per_variant: str = "auto"  # auto | full | fast | off
    flank: int = 5000
    sequences_in_tsv: bool = True
    all_transcripts: bool = False
    include_noncoding: bool = False
    only_biotypes: Optional[Set[str]] = None
    only_transcripts: Optional[Set[str]] = None
    max_transcripts: Optional[int] = None
    min_protein_change: bool = False
    verbose: bool = True


class GenomeSource:
    def __init__(self, genome: FastaReader):
        self.genome = genome

    def extract(self, chrom: str, start0: int, end0: int, strand: str) -> str:
        seq = self.genome.fetch(chrom, start0, end0)
        return reverse_complement(seq) if strand == "-" else seq


class RegionSource:
    def __init__(self, region: MutatedRegion):
        self.region = region

    def extract(self, chrom: str, start0: int, end0: int, strand: str) -> str:
        return self.region.extract_segment(start0, end0, strand)


@dataclass
class ReferenceBundle:
    tx: TranscriptModel
    pad5: int = 0
    stop_segments: List[Segment] = field(default_factory=list)
    cds3_segments: List[Segment] = field(default_factory=list)
    utr3_segments: List[Segment] = field(default_factory=list)
    genomic_tail_cds: Optional[Segment] = None
    genomic_tail_exon: Optional[Segment] = None
    ref_cdna: str = ""
    ref_cds: str = ""
    ref_cds_ext: str = ""
    ref_protein: str = ""          # annotated CDS translation (primary)
    ref_protein_annot: str = ""    # alias for ref_protein
    ref_protein_ext: str = ""      # CDS + 3' extension (readthrough)
    ref_protein_stop: bool = False
    ref_protein_ext_stop: bool = False
    ref_cdna_ensembl: str = ""
    ref_cds_ensembl: str = ""
    ref_protein_ensembl: str = ""
    cdna_matches: Optional[bool] = None
    cds_matches: Optional[bool] = None
    protein_matches: Optional[bool] = None
    reference_class: str = ""
    stop_reason: str = ""
    notes: List[str] = field(default_factory=list)
    stats: Dict[str, object] = field(default_factory=dict)


@dataclass
class MutatedBundle:
    alt_cdna: str = ""
    alt_cds: str = ""
    alt_cds_ext: str = ""
    alt_protein: str = ""          # annotated CDS translation (primary)
    alt_protein_annot: str = ""    # alias for alt_protein
    alt_protein_ext: str = ""      # CDS + 3' extension (readthrough)
    alt_protein_stop: bool = False
    alt_protein_ext_stop: bool = False
    n_edits: int = 0
    notes: List[str] = field(default_factory=list)


def _concat(source, chrom: str, segments: Sequence[Segment]) -> str:
    return "".join(source.extract(chrom, s.start0, s.end0, s.strand) for s in segments)


def _clip_segment(seg: Segment, max_len: int, from_end: bool) -> Segment:
    """Return the first/last max_len bases of a segment in transcript order."""
    if max_len <= 0 or seg.length <= max_len:
        return seg
    if from_end:
        if seg.strand == "+":
            return Segment(seg.chrom, seg.end - max_len + 1, seg.end, seg.strand)
        return Segment(seg.chrom, seg.start, seg.start + max_len - 1, seg.strand)
    if seg.strand == "+":
        return Segment(seg.chrom, seg.start, seg.start + max_len - 1, seg.strand)
    return Segment(seg.chrom, seg.end - max_len + 1, seg.end, seg.strand)


def _segment_after_cds(tx: TranscriptModel, last_cds_end: int, length: int) -> Segment:
    """Genomic segment immediately 3' of the CDS in transcript orientation."""
    if tx.strand == "+":
        return Segment(tx.chrom, last_cds_end + 1, last_cds_end + length, "+")
    return Segment(tx.chrom, last_cds_end - length, last_cds_end - 1, "-")


def _segment_after_exon(tx: TranscriptModel, length: int) -> Optional[Segment]:
    if not tx.exons:
        return None
    if tx.strand == "+":
        last = max(tx.exons, key=lambda e: e.end)
        return Segment(tx.chrom, last.end + 1, last.end + length, "+")
    last = min(tx.exons, key=lambda e: e.start)
    return Segment(tx.chrom, last.start - length, last.start - 1, "-")


def build_reference_bundle(
    tx: TranscriptModel,
    genome: FastaReader,
    options: RunOptions,
    ref_cdna: str = "",
    ref_cds: str = "",
    ref_protein: str = "",
) -> ReferenceBundle:
    """Reconstruct the reference cDNA/CDS/protein and classify the model."""
    b = ReferenceBundle(tx=tx)
    b.ref_cdna_ensembl = ref_cdna or ""
    b.ref_cds_ensembl = ref_cds or ""
    b.ref_protein_ensembl = ref_protein or ""
    source = GenomeSource(genome)
    try:
        b.ref_cdna = _concat(source, tx.chrom, tx.exons)
    except KeyError:
        b.notes.append(f"chromosome {tx.chrom} missing from genome")
        b.ref_cdna = ""

    if ref_cdna:
        b.cdna_matches = b.ref_cdna == ref_cdna
        if not b.cdna_matches:
            b.notes.append(
                f"cDNA differs from Ensembl ({len(b.ref_cdna)} vs {len(ref_cdna)} bp)"
            )

    if not tx.cds:
        b.reference_class = "no_cds"
        b.stop_reason = "no CDS annotation"
        return b

    cds_core = _concat(source, tx.chrom, tx.cds)
    cds_len = sum(s.length for s in tx.cds)
    b.stats["cds_core_len"] = cds_len

    # Determine whether the CDS already contains the stop codon.
    stop_in_cds = False
    stop_segments: List[Segment] = []
    if tx.stops:
        # Ensembl-style GTF: stop_codon is a separate row adjacent to the CDS.
        overlaps = any(
            st.overlaps0(c.start0, c.end0) for st in tx.stops for c in tx.cds
        )
        if overlaps:
            stop_in_cds = True
        else:
            stop_segments = sort_segments(tx.stops, tx.strand)
    else:
        if len(cds_core) % 3 == 0 and cds_core[-3:] in STOP_CODONS:
            stop_in_cds = True
    b.stop_segments = stop_segments

    # 5' N padding for partial models: the Ensembl CDS FASTA is authoritative
    # when available; otherwise the first CDS phase gives the partial codon.
    pad5 = 0
    if ref_cds:
        pad5 = len(ref_cds) - len(ref_cds.lstrip("N"))
    elif tx.is_partial_start and tx.first_phase is not None:
        pad5 = (3 - tx.first_phase) % 3
    b.pad5 = pad5

    core = "N" * pad5 + cds_core
    if stop_segments:
        core += _concat(source, tx.chrom, stop_segments)
    b.stats["core_len"] = len(core)

    # Some partial models have 1-3 bases beyond the GTF CDS/stop (a partial or
    # complete stop codon).  Infer those bases from the genome immediately 3'
    # of the annotated coding end, but only when they make the reconstruction
    # match Ensembl's CDS FASTA.
    cds3_segments: List[Segment] = []
    extra = len(ref_cds) - len(core) if ref_cds else 0
    if 0 < extra <= 3:
        plus_end = max(
            [c.end for c in tx.cds] + [s.end for s in stop_segments]
        )
        minus_end = min(
            [c.start for c in tx.cds] + [s.start for s in stop_segments]
        )
        anchor_end = plus_end if tx.strand == "+" else minus_end
        candidate = _segment_after_cds(tx, anchor_end, extra)
        candidate_seq = source.extract(
            candidate.chrom, candidate.start0, candidate.end0, candidate.strand
        )
        if not ref_cds or (core + candidate_seq) == ref_cds:
            cds3_segments = [candidate]
            b.stats["inferred_3p_extension"] = extra
        else:
            b.notes.append(
                "could not reproduce the 3' end of the Ensembl CDS from the "
                f"genome (extra={extra} bp)"
            )
    elif extra < 0:
        b.stats["ref_cds_shorter_than_gtf"] = abs(extra)
    elif extra > 3:
        b.notes.append(
            f"Ensembl CDS is {extra} bp longer than the GTF reconstruction; "
            "not inferred automatically"
        )
    b.cds3_segments = cds3_segments
    if cds3_segments:
        core += _concat(source, tx.chrom, cds3_segments)

    b.ref_cds = core
    if ref_cds:
        b.cds_matches = core == ref_cds
        if not b.cds_matches:
            b.notes.append(
                f"CDS differs from Ensembl: reconstruction={len(core)} bp, "
                f"Ensembl={len(ref_cds)} bp (pad5={pad5}, "
                f"core={len(cds_core)}, extra={extra})"
            )

    # 3' UTR / downstream extension.
    utr3 = list(tx.utr3)
    if not utr3 and options.extend_mode in ("utr", "utr+genomic"):
        # Derive from exons if the GTF did not provide explicit UTR rows.
        derive_utrs(tx)
        utr3 = list(tx.utr3)
    b.utr3_segments = utr3

    if tx.exons:
        tail = _segment_after_exon(tx, options.extend_bp)
        if tail is not None:
            chrom_len = genome.length(tx.chrom)
            if tail.start < 1 or tail.end > chrom_len:
                tail = Segment(
                    tail.chrom,
                    max(1, tail.start),
                    min(chrom_len, tail.end),
                    tail.strand,
                )
            if tail.end >= tail.start:
                b.genomic_tail_exon = tail
    if stop_segments:
        plus_end = max([c.end for c in tx.cds] + [s.end for s in stop_segments])
        minus_end = min([c.start for c in tx.cds] + [s.start for s in stop_segments])
    else:
        plus_end = max(c.end for c in tx.cds)
        minus_end = min(c.start for c in tx.cds)
    anchor_end = plus_end if tx.strand == "+" else minus_end
    tail = _segment_after_cds(tx, anchor_end, options.extend_bp)
    chrom_len = genome.length(tx.chrom)
    if tail.start < 1 or tail.end > chrom_len:
        tail = Segment(
            tail.chrom, max(1, tail.start), min(chrom_len, tail.end), tail.strand
        )
    if tail.end >= tail.start:
        b.genomic_tail_cds = tail

    b.ref_cds_ext = b.ref_cds + _extension_sequence(
        source, b, options, use_exon_tail=options.extend_mode != "genomic"
    )
    b.ref_protein, b.ref_protein_stop = translate(b.ref_cds, stop_at_stop=True)
    b.ref_protein_annot = b.ref_protein
    b.ref_protein_ext, b.ref_protein_ext_stop = translate(
        b.ref_cds_ext, stop_at_stop=True
    )
    if ref_protein:
        b.protein_matches = b.ref_protein_annot.rstrip("*") == ref_protein.rstrip("*")
        if not b.protein_matches:
            b.notes.append(
                f"protein differs from Ensembl ({len(b.ref_protein_annot)} vs "
                f"{len(ref_protein)} aa)"
            )

    # Reference model classification; this is the answer to "why no stop?".
    # Partial-model tags are reported even when a stop is present, because the
    # 5' N padding and the 3' incompleteness are independent of stop status.
    has_stop = b.ref_protein_stop
    if tx.is_partial_start and tx.is_partial_end:
        b.reference_class = "partial_5_3"
        b.stop_reason = "both 5' and 3' ends incomplete (mRNA_start_NF/mRNA_end_NF)"
    elif tx.is_partial_start:
        b.reference_class = "partial_5_complete_3" if has_stop else "partial_5"
        b.stop_reason = (
            "mRNA_start_NF/cds_start_NF: 5' incomplete"
            + (" (stop present)" if has_stop else " and no stop")
        )
    elif tx.is_partial_end:
        b.reference_class = "partial_3_with_stop" if has_stop else "partial_3"
        b.stop_reason = (
            "mRNA_end_NF/cds_end_NF: 3' incomplete"
            + (" (stop present)" if has_stop else " and no stop")
        )
    elif tx.is_nsd:
        b.reference_class = "non_stop_decay"
        b.stop_reason = "biotype non_stop_decay (no stop; degraded by NSD)"
    elif tx.is_immune_segment:
        b.reference_class = "immune_segment"
        b.stop_reason = f"biotype {tx.biotype}; V/D/J gene segment, not a full protein"
    elif has_stop:
        b.reference_class = "complete"
        b.stop_reason = "annotated stop codon present"
    elif tx.is_nmd:
        b.reference_class = "nmd_partial"
        b.stop_reason = "nonsense_mediated_decay model without an annotated stop"
    else:
        b.reference_class = "no_stop_other"
        b.stop_reason = "no in-frame stop in the available CDS/extension"
    b.stats["has_stop"] = b.ref_protein_stop
    b.stats["leading_n"] = pad5
    b.stats["stop_in_cds"] = stop_in_cds
    return b


def _extension_sequence(
    source, bundle: ReferenceBundle, options: RunOptions, use_exon_tail: bool
) -> str:
    if options.extend_mode == "none" or options.extend_bp <= 0:
        return ""
    tx = bundle.tx
    remaining = options.extend_bp
    pieces: List[str] = []
    if options.extend_mode in ("utr", "utr+genomic"):
        for seg in bundle.utr3_segments:
            if remaining <= 0:
                break
            clipped = _clip_segment(seg, remaining, from_end=False)
            seq = source.extract(clipped.chrom, clipped.start0, clipped.end0, clipped.strand)
            if len(seq) > remaining:
                seq = seq[:remaining]
            pieces.append(seq)
            remaining -= len(seq)
    if options.extend_mode in ("genomic", "utr+genomic") and remaining > 0:
        tail = bundle.genomic_tail_exon if use_exon_tail else bundle.genomic_tail_cds
        if tail is not None:
            seq = source.extract(tail.chrom, tail.start0, tail.end0, tail.strand)
            # seq is already transcript-oriented; keep the first `remaining`
            # bases of the 3' tail for both strands.
            if len(seq) > remaining:
                seq = seq[:remaining]
            pieces.append(seq)
    return "".join(pieces)


def mutate_bundle(
    bundle: ReferenceBundle,
    region: MutatedRegion,
    options: RunOptions,
) -> MutatedBundle:
    tx = bundle.tx
    source = RegionSource(region)
    out = MutatedBundle(n_edits=len(region.edits))
    try:
        out.alt_cdna = _concat(source, tx.chrom, tx.exons)
    except KeyError:
        out.notes.append(f"chromosome {tx.chrom} missing from genome")
    cds_core = _concat(source, tx.chrom, tx.cds)
    stop = _concat(source, tx.chrom, bundle.stop_segments)
    ext3 = _concat(source, tx.chrom, bundle.cds3_segments)
    out.alt_cds = "N" * bundle.pad5 + cds_core + stop + ext3
    out.alt_cds_ext = out.alt_cds + _extension_sequence(
        source, bundle, options, use_exon_tail=options.extend_mode != "genomic"
    )
    out.alt_protein, out.alt_protein_stop = translate(
        out.alt_cds, stop_at_stop=True
    )
    out.alt_protein_annot = out.alt_protein
    out.alt_protein_ext, out.alt_protein_ext_stop = translate(
        out.alt_cds_ext, stop_at_stop=True
    )
    return out


def _first_codon_diff(
    ref_cds: str, alt_cds: str
) -> Tuple[Optional[int], str, str, str, str]:
    """Return first differing codon and amino-acid change."""
    n = min(len(ref_cds), len(alt_cds))
    for i in range(0, n - 2, 3):
        rc = ref_cds[i : i + 3]
        ac = alt_cds[i : i + 3]
        if rc != ac:
            ra = translate(rc, stop_at_stop=False)[0]
            aa = translate(ac, stop_at_stop=False)[0]
            return i // 3 + 1, rc, ac, ra, aa
    if len(ref_cds) != len(alt_cds):
        i = (n // 3) * 3
        rc = ref_cds[i : i + 3]
        ac = alt_cds[i : i + 3]
        ra = translate(rc, stop_at_stop=False)[0]
        aa = translate(ac, stop_at_stop=False)[0]
        return i // 3 + 1, rc, ac, ra, aa
    return None, "", "", "", ""


@dataclass
class VariantLocation:
    region: str = ""
    exon_rank: int = 0
    in_cds: bool = False
    in_stop: bool = False
    in_utr5: bool = False
    in_utr3: bool = False
    splice: str = ""
    distance: Optional[int] = None
    note: str = ""


def locate_edit(
    tx: TranscriptModel,
    edit: Edit,
    flank: int = 5000,
    bundle: Optional[ReferenceBundle] = None,
) -> VariantLocation:
    loc = VariantLocation()
    a, b = edit.start0, edit.end0
    cds_segments = list(tx.cds)
    stop_segments = list(tx.stops)
    if bundle is not None:
        cds_segments += list(bundle.cds3_segments)
        stop_segments = list(bundle.stop_segments)
    if edit.is_insertion:
        p = a
        in_exon = any(
            e.start0 <= p <= e.end0 for e in tx.exons
        )
        in_cds = any(c.start0 <= p <= c.end0 for c in cds_segments)
        in_stop = any(s.start0 <= p <= s.end0 for s in stop_segments)
        in_utr5 = any(u.start0 <= p <= u.end0 for u in tx.utr5)
        in_utr3 = any(u.start0 <= p <= u.end0 for u in tx.utr3)
    else:
        in_exon = any(e.overlaps0(a, b) for e in tx.exons)
        in_cds = any(c.overlaps0(a, b) for c in cds_segments)
        in_stop = any(s.overlaps0(a, b) for s in stop_segments)
        in_utr5 = any(u.overlaps0(a, b) for u in tx.utr5)
        in_utr3 = any(u.overlaps0(a, b) for u in tx.utr3)

    loc.in_stop = in_stop
    # The stop codon is coding for consequence purposes.
    loc.in_cds = in_cds or in_stop
    loc.in_utr5 = in_utr5 and not loc.in_cds
    loc.in_utr3 = in_utr3 and not loc.in_cds

    # Splice proximity: 2 bp into the intron, 3 bp into the exon.
    for i, exon in enumerate(tx.exons, 1):
        # In transcript orientation the 5' exon end is the acceptor (unless it
        # is the first exon) and the 3' end is the donor (unless it is the last
        # exon).  For a minus-strand gene the genomic coordinates are swapped.
        boundaries = []
        if tx.strand == "+":
            if i > 1:
                boundaries.append((exon.start0, "acceptor"))
            if i < len(tx.exons):
                boundaries.append((exon.end0, "donor"))
        else:
            if i > 1:
                boundaries.append((exon.end0, "acceptor"))
            if i < len(tx.exons):
                boundaries.append((exon.start0, "donor"))
        for boundary, kind in boundaries:
            if edit.is_insertion:
                dist = abs(a - boundary)
                in_region = a >= boundary - 2 and a <= boundary + 2
            else:
                dist = min(abs(a - boundary), abs(b - boundary))
                in_region = (a <= boundary + 2 and b >= boundary - 2)
            if in_region and dist <= 5:
                loc.splice = kind
                loc.note = f"within {dist} bp of a splice {kind} site"
                break

    if in_exon:
        for i, exon in enumerate(tx.exons, 1):
            if (edit.is_insertion and exon.start0 <= a <= exon.end0) or (
                not edit.is_insertion and exon.overlaps0(a, b)
            ):
                loc.exon_rank = i
                break
        if loc.in_cds:
            loc.region = "exon_cds"
        elif loc.in_utr5:
            loc.region = "exon_5utr"
        elif loc.in_utr3:
            loc.region = "exon_3utr"
        else:
            loc.region = "exon_noncoding"
    elif tx.start0 <= (a if edit.is_insertion else a) <= tx.end0:
        loc.region = "intron"
    elif edit.is_insertion:
        if a < tx.start0:
            loc.region = "upstream"
            loc.distance = tx.start0 - a
        else:
            loc.region = "downstream"
            loc.distance = a - tx.end0
    else:
        if b <= tx.start0:
            loc.region = "upstream"
            loc.distance = tx.start0 - b
        elif a >= tx.end0:
            loc.region = "downstream"
            loc.distance = a - tx.end0
        else:
            loc.region = "intron"
    return loc


def _classify_snv(
    tx: TranscriptModel,
    bundle: ReferenceBundle,
    edit: Edit,
) -> Dict[str, object]:
    """Fast codon annotation for a SNV in the CDS."""
    out: Dict[str, object] = {}
    if len(edit.ref) != 1 or len(edit.alt) != 1:
        return out
    # Build the CDS index map lazily for this transcript.
    segments = list(tx.cds) + list(bundle.stop_segments) + list(bundle.cds3_segments)
    offset = bundle.pad5
    idx: Optional[int] = None
    for seg in segments:
        if seg.start0 <= edit.start0 < seg.end0:
            if tx.strand == "+":
                idx = offset + (edit.start0 - seg.start0)
            else:
                idx = offset + (seg.end0 - 1 - edit.start0)
            break
        offset += seg.length
    if idx is None:
        return out
    codon_start = (idx // 3) * 3
    codon_ref = bundle.ref_cds[codon_start : codon_start + 3]
    if len(codon_ref) < 3:
        return out
    alt_base = edit.alt if tx.strand == "+" else reverse_complement(edit.alt)
    within = idx % 3
    codon_alt = codon_ref[:within] + alt_base + codon_ref[within + 1 :]
    aa_ref, _ = translate(codon_ref)
    aa_alt, _ = translate(codon_alt)
    out.update(
        {
            "cds_pos": idx + 1,
            "protein_pos": idx // 3 + 1,
            "codon_ref": codon_ref,
            "codon_alt": codon_alt,
            "aa_ref": aa_ref,
            "aa_alt": aa_alt,
        }
    )
    if aa_ref == aa_alt:
        out["consequence"] = "synonymous"
    elif aa_alt == "*":
        out["consequence"] = "nonsense"
    elif idx == 0 and codon_alt not in ("ATG", "CTG", "TTG", "GTG"):
        out["consequence"] = "start_lost"
    elif aa_ref == "*" and aa_alt != "*":
        out["consequence"] = "stop_lost"
    else:
        out["consequence"] = "missense"
    return out


def _single_variant(
    tx: TranscriptModel,
    bundle: ReferenceBundle,
    edit: Edit,
    genome: FastaReader,
    options: RunOptions,
) -> Tuple[MutatedBundle, Dict[str, object]]:
    """Apply only one edit to the transcript and return its effect."""
    start0 = min(tx.start0, edit.start0)
    end0 = max(tx.end0, edit.end0)
    region = MutatedRegion(genome, tx.chrom, [edit], start0, end0)
    mb = mutate_bundle(bundle, region, options)
    extra: Dict[str, object] = {}
    pos, codon_ref, codon_alt, aa_ref, aa_alt = _first_codon_diff(
        bundle.ref_cds, mb.alt_cds
    )
    if pos is not None:
        extra["protein_pos"] = pos
        extra["cds_pos"] = (pos - 1) * 3 + 1
        extra["codon_ref"] = codon_ref
        extra["codon_alt"] = codon_alt
        extra["aa_ref"] = aa_ref
        extra["aa_alt"] = aa_alt
    extra["protein_diff"] = protein_diff(bundle.ref_protein, mb.alt_protein)
    extra["protein_ext_diff"] = protein_diff(
        bundle.ref_protein_ext, mb.alt_protein_ext
    )
    if bundle.ref_protein != mb.alt_protein:
        import difflib

        sm = difflib.SequenceMatcher(
            None, bundle.ref_protein, mb.alt_protein, autojunk=False
        )
        # Report the first aligned difference for the TSV; the complete
        # proteins are in the transcript report.
        for tag, i1, i2, j1, j2 in sm.get_opcodes():
            if tag != "equal":
                extra["aa_diff_span"] = f"{i1 + 1}-{i2}->{j1 + 1}-{j2}"
                break
    return mb, extra


def annotate_edit(
    tx: TranscriptModel,
    bundle: ReferenceBundle,
    edit: Edit,
    genome: FastaReader,
    options: RunOptions,
) -> Dict[str, object]:
    loc = locate_edit(tx, edit, options.flank, bundle=bundle)
    result: Dict[str, object] = {
        "region": loc.region,
        "exon_rank": loc.exon_rank,
        "splice": loc.splice,
        "variant_distance": "" if loc.distance is None else loc.distance,
        "location_note": loc.note,
        "consequence": "",
        "cds_pos": "",
        "protein_pos": "",
        "codon_ref": "",
        "codon_alt": "",
        "aa_ref": "",
        "aa_alt": "",
        "protein_diff": "",
        "stop_gained": False,
        "stop_lost": False,
        "start_lost": False,
    }
    if edit.unresolved:
        result["consequence"] = "unresolved"
        return result
    if loc.in_cds:
        if options.per_variant == "fast" and edit.kind == "snv":
            result.update(_classify_snv(tx, bundle, edit))
        elif options.per_variant == "fast":
            delta = edit.length_delta
            if delta % 3:
                result["consequence"] = "frameshift"
            elif delta < 0:
                result["consequence"] = "inframe_deletion"
            elif delta > 0:
                result["consequence"] = "inframe_insertion"
            else:
                result["consequence"] = "coding_variant"
        elif options.per_variant == "full":
            try:
                mb, extra = _single_variant(tx, bundle, edit, genome, options)
            except Exception as exc:  # keep the report useful on edge cases
                result["consequence"] = "effect_error"
                result["location_note"] = f"{loc.note}; {exc}"
                return result
            result.update(extra)
            # Use the CDS length change, not the genomic length change.  A
            # deletion spanning an intron can be in-frame even when its VCF
            # REF/ALT length difference is not divisible by three.
            cds_delta = len(mb.alt_cds) - len(bundle.ref_cds)
            result["cds_length_delta"] = cds_delta
            # Derive a consequence from the single-variant comparison.
            if bundle.ref_protein == mb.alt_protein:
                # A coding indel can be in-frame and synonymous only if the
                # amino-acid sequence is unchanged; still distinguish the DNA
                # change from a true synonymous SNV.
                if edit.kind == "snv":
                    result["consequence"] = "synonymous"
                elif cds_delta == 0:
                    result["consequence"] = "synonymous"
                else:
                    result["consequence"] = "coding_inframe"
            else:
                aa_ref = str(result.get("aa_ref", ""))
                aa_alt = str(result.get("aa_alt", ""))
                if cds_delta % 3:
                    result["consequence"] = "frameshift"
                elif aa_alt == "*":
                    result["consequence"] = "nonsense"
                elif aa_ref == "*" and aa_alt != "*":
                    result["consequence"] = "stop_lost"
                elif (
                    result.get("protein_pos") == 1
                    and str(result.get("codon_ref", ""))
                    in ("ATG", "CTG", "TTG", "GTG")
                    and str(result.get("codon_alt", ""))
                    not in ("ATG", "CTG", "TTG", "GTG")
                ):
                    result["consequence"] = "start_lost"
                elif cds_delta < 0:
                    result["consequence"] = "inframe_deletion"
                elif cds_delta > 0:
                    result["consequence"] = "inframe_insertion"
                elif edit.kind in ("snv", "mnv"):
                    result["consequence"] = "missense"
                else:
                    result["consequence"] = "protein_altering"
            # Stop/start flags are orthogonal to the primary consequence: a
            # frameshift can also introduce a premature stop.  Do not overwrite
            # the primary consequence with "stop_gained".
            if bundle.ref_protein_stop and not mb.alt_protein_stop:
                result["stop_lost"] = True
            if (not bundle.ref_protein_stop and mb.alt_protein_stop) or (
                mb.alt_protein_stop
                and len(mb.alt_protein_annot) < len(bundle.ref_protein_annot)
                and not result["stop_lost"]
            ):
                result["stop_gained"] = True
            if result["consequence"] in ("start_lost", "nonsense"):
                if result["consequence"] == "start_lost":
                    result["start_lost"] = True
                if result["consequence"] == "nonsense":
                    result["stop_gained"] = True
            if result["consequence"] == "stop_lost":
                result["stop_lost"] = True
        else:
            result["consequence"] = "coding_variant"
    elif loc.splice:
        result["consequence"] = f"splice_{loc.splice}"
    elif loc.in_utr5:
        result["consequence"] = "5_prime_UTR"
    elif loc.in_utr3:
        result["consequence"] = "3_prime_UTR"
    elif loc.region == "intron":
        result["consequence"] = "intron"
    elif loc.region == "upstream":
        result["consequence"] = "upstream_gene"
    elif loc.region == "downstream":
        result["consequence"] = "downstream_gene"
    elif loc.region == "exon_noncoding":
        result["consequence"] = "exonic_noncoding"
    else:
        result["consequence"] = "intergenic"
    if result["consequence"] == "nonsense":
        result["stop_gained"] = True
    elif result["consequence"] == "stop_lost":
        result["stop_lost"] = True
    elif result["consequence"] == "start_lost":
        result["start_lost"] = True
    return result


def transcript_aggregate(
    bundle: ReferenceBundle,
    mb: MutatedBundle,
    variants: Sequence[Edit],
    consequences: Sequence[str],
    options: RunOptions,
) -> Dict[str, object]:
    tx = bundle.tx
    ref = bundle.ref_protein
    alt = mb.alt_protein
    aggregate: Dict[str, object] = {
        "n_variants": len(variants),
        "protein_changed": ref != alt,
        "protein_diff": protein_diff(ref, alt),
        "protein_length_ref": len(ref),
        "protein_length_alt": len(alt),
        "protein_ext_changed": bundle.ref_protein_ext != mb.alt_protein_ext,
        "protein_ext_diff": protein_diff(bundle.ref_protein_ext, mb.alt_protein_ext),
        "protein_ext_length_ref": len(bundle.ref_protein_ext),
        "protein_ext_length_alt": len(mb.alt_protein_ext),
        "protein_stop_ref": bundle.ref_protein_stop,
        "protein_stop_alt": mb.alt_protein_stop,
        "protein_stop_ext_ref": bundle.ref_protein_ext_stop,
        "protein_stop_ext_alt": mb.alt_protein_ext_stop,
        "stop_lost": bool(bundle.ref_protein_stop and not mb.alt_protein_stop),
        "stop_gained": bool(
            (not bundle.ref_protein_stop and mb.alt_protein_stop)
            or (mb.alt_protein_stop and len(mb.alt_protein_annot) < len(bundle.ref_protein_annot))
        ),
        "start_lost": False,
        "frameshift": False,
        "consequences": ",".join(sorted(set(c for c in consequences if c))),
    }
    consequence_set = set(consequences)
    # Per-variant consequences are more informative than the aggregate alone
    # when several variants hit the same transcript (e.g. a stop-lost plus a
    # separate premature stop).
    if "stop_lost" in consequence_set:
        aggregate["stop_lost"] = True
    if {"nonsense", "stop_gained"} & consequence_set:
        aggregate["stop_gained"] = True
    if "start_lost" in consequence_set:
        aggregate["start_lost"] = True
    if "frameshift" in consequence_set:
        aggregate["frameshift"] = True
    # Start lost: the first codon was a recognized start and is no longer one.
    # Do not use the translated amino acid alone: CTG/TTG/GTG are alternative
    # starts and should not be reported as start-lost merely because they
    # translate to L/V rather than M.
    start_codons = {"ATG", "CTG", "TTG", "GTG"}
    if (
        bundle.ref_cds[:3] in start_codons
        and mb.alt_cds[:3] not in start_codons
    ):
        aggregate["start_lost"] = True
    # Frameshift: use the net CDS length change.  This correctly ignores
    # intronic parts of a deletion and is accurate for single edits.
    cds_delta = len(mb.alt_cds) - len(bundle.ref_cds)
    aggregate["cds_length_delta"] = cds_delta
    if cds_delta % 3 != 0:
        aggregate["frameshift"] = True
    return aggregate


class TranscriptIndex:
    """Bin index for finding transcripts overlapping a variant.

    A transcript can start far before a variant and still contain it, so a
    simple bisect on transcript starts is not sufficient.  Fixed-size bins
    (16 kb) keep the index small and make both short variants and large SVs
    cheap to query.
    """

    BIN_SIZE = 16384

    def __init__(self, annotation: Annotation, transcripts: Optional[Set[str]] = None):
        self.by_chrom: Dict[str, List[TranscriptModel]] = {}
        self.bins: Dict[str, Dict[int, List[TranscriptModel]]] = {}
        for tx in annotation.iter_transcripts():
            if transcripts is not None and tx.transcript_id not in transcripts:
                continue
            self.by_chrom.setdefault(tx.chrom, []).append(tx)
            chrom_bins = self.bins.setdefault(tx.chrom, {})
            b0 = tx.start0 // self.BIN_SIZE
            b1 = max(b0, (tx.end0 - 1) // self.BIN_SIZE)
            for b in range(b0, b1 + 1):
                chrom_bins.setdefault(b, []).append(tx)

    def query(self, chrom: str, start0: int, end0: int, flank: int = 0) -> List[TranscriptModel]:
        if chrom not in self.bins:
            return []
        a = max(0, start0 - flank)
        b = max(a, end0 + flank)
        if b == a:
            b = a  # insertion point
        b0 = a // self.BIN_SIZE
        b1 = max(b0, (b - 1) // self.BIN_SIZE) if b > a else b0
        chrom_bins = self.bins[chrom]
        seen: Set[str] = set()
        out: List[TranscriptModel] = []
        for bin_id in range(b0, b1 + 1):
            for tx in chrom_bins.get(bin_id, []):
                if tx.transcript_id in seen:
                    continue
                # For a point (insertion) b == a; use an inclusive test.
                if tx.start0 <= b and tx.end0 >= a:
                    seen.add(tx.transcript_id)
                    out.append(tx)
        return out


def load_reference_sequences(
    path: str, wanted: Set[str]
) -> Dict[str, str]:
    """Load versionless FASTA IDs, but only the wanted transcripts."""
    out: Dict[str, str] = {}
    if not path:
        return out
    from .fasta import iter_fasta

    for name, seq in iter_fasta(path):
        base = name.split(".")[0]
        if base in wanted:
            out[base] = seq
    return out


def load_reference_proteins(
    path: str, wanted: Set[str]
) -> Dict[str, str]:
    """Load Ensembl peptide FASTA, keyed by transcript ID.

    Ensembl pep.all.fa headers are ``>ENSP... ... transcript:ENST...``; the
    first field is a protein ID, so a naive first-field lookup would miss every
    transcript.  Fall back to the first field for simple test FASTA files.
    """
    import re
    from .fasta import iter_fasta_records

    out: Dict[str, str] = {}
    if not path:
        return out
    pattern = re.compile(r"transcript:(\S+)")
    for header, seq in iter_fasta_records(path):
        first = header.split()[0].split(".")[0]
        match = pattern.search(header)
        key = match.group(1).split(".")[0] if match else first
        if key in wanted:
            out[key] = seq
    return out
