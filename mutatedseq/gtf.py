"""Minimal GTF/GFF3 parser and transcript model builder.

Only the features needed for transcript/CDS reconstruction are retained:
gene, transcript/mRNA, exon, CDS, stop_codon and UTR rows.  Ensembl GTF uses
1-based inclusive coordinates; that convention is preserved in Segment.
"""

from __future__ import annotations

import gzip
import re
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .model import Segment, TranscriptModel, derive_utrs, sort_segments


_ATTR_GTF = re.compile(r"([A-Za-z0-9_.:-]+)\s+\"([^\"]*)\"")
_ATTR_GFF = re.compile(r"([A-Za-z0-9_.:-]+)=([^;]+)")


def open_maybe_gzip(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def parse_attributes(text: str) -> Dict[str, str]:
    """Parse GTF (key "value") or GFF3 (key=value) attributes."""
    attrs: Dict[str, str] = {}
    matches = _ATTR_GTF.findall(text)
    if not matches:
        matches = [
            (key, value.strip())
            for key, value in _ATTR_GFF.findall(text)
        ]
    for key, value in matches:
        # Duplicate keys such as tag are joined; callers split on comma.
        if key in attrs:
            attrs[key] = attrs[key] + "," + value
        else:
            attrs[key] = value
    return attrs


def _frame(value: str) -> Optional[int]:
    value = value.strip()
    if value in ("0", "1", "2"):
        return int(value)
    return None


def _exon_number(value: str) -> Optional[int]:
    value = value.strip()
    if value.isdigit():
        return int(value)
    return None


class Annotation:
    """In-memory annotation.

    ``transcripts`` maps transcript_id (versionless, i.e. ENST0000...) to a
    TranscriptModel.  Versioned IDs are accepted in lookups through
    ``get_transcript``.
    """

    def __init__(self) -> None:
        self.transcripts: Dict[str, TranscriptModel] = {}
        self.gene_names: Dict[str, str] = {}
        self.gene_biotypes: Dict[str, str] = {}

    def __len__(self) -> int:
        return len(self.transcripts)

    def get_transcript(self, transcript_id: str) -> Optional[TranscriptModel]:
        if not transcript_id:
            return None
        if transcript_id in self.transcripts:
            return self.transcripts[transcript_id]
        base = transcript_id.split(".")[0]
        return self.transcripts.get(base)

    def iter_transcripts(self) -> Iterable[TranscriptModel]:
        return self.transcripts.values()

    def by_chrom(self) -> Dict[str, List[TranscriptModel]]:
        out: Dict[str, List[TranscriptModel]] = defaultdict(list)
        for tx in self.transcripts.values():
            out[tx.chrom].append(tx)
        for tx_list in out.values():
            tx_list.sort(key=lambda t: (t.start0, t.end0, t.transcript_id))
        return out


def parse_gtf(path: str, verbose: bool = True) -> Annotation:
    """Parse gene/transcript/exon/CDS/stop_codon/UTR rows into an Annotation."""
    ann = Annotation()

    # transcript-level metadata, keyed by transcript_id
    tx_meta: Dict[str, Dict[str, object]] = {}
    tx_exons: Dict[str, List[Segment]] = defaultdict(list)
    tx_cds: Dict[str, List[Segment]] = defaultdict(list)
    tx_stops: Dict[str, List[Segment]] = defaultdict(list)
    tx_utr5: Dict[str, List[Segment]] = defaultdict(list)
    tx_utr3: Dict[str, List[Segment]] = defaultdict(list)
    tx_chrom: Dict[str, str] = {}
    tx_strand: Dict[str, str] = {}
    tx_gene_id: Dict[str, str] = {}

    line_no = 0
    with open_maybe_gzip(path) as handle:
        for line in handle:
            line_no += 1
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) < 9:
                continue
            chrom, source, feature, start_s, end_s, score, strand, frame_s, attr_s = (
                fields[:9]
            )
            if feature not in {
                "gene",
                "transcript",
                "mRNA",
                "exon",
                "CDS",
                "stop_codon",
                "five_prime_utr",
                "three_prime_utr",
            }:
                continue
            try:
                start = int(start_s)
                end = int(end_s)
            except ValueError:
                continue
            attrs = parse_attributes(attr_s)
            gene_id = attrs.get("gene_id", "")
            transcript_id = attrs.get("transcript_id", "")
            if feature == "gene":
                if gene_id:
                    if attrs.get("gene_name"):
                        ann.gene_names[gene_id] = attrs["gene_name"]
                    if attrs.get("gene_biotype"):
                        ann.gene_biotypes[gene_id] = attrs["gene_biotype"]
                continue
            if not transcript_id:
                continue
            if feature in ("transcript", "mRNA"):
                tx_meta[transcript_id] = {
                    "gene_id": gene_id,
                    "gene_name": attrs.get("gene_name", ""),
                    "gene_biotype": attrs.get("gene_biotype", ""),
                    "transcript_name": attrs.get("transcript_name", ""),
                    "biotype": attrs.get("transcript_biotype", "")
                    or attrs.get("biotype", ""),
                    "source": source,
                    "tags": [t for t in attrs.get("tag", "").split(",") if t],
                    "start": start,
                    "end": end,
                    "chrom": chrom,
                    "strand": strand,
                }
                tx_chrom[transcript_id] = chrom
                tx_strand[transcript_id] = strand
                continue
            if feature in ("five_prime_utr", "three_prime_utr"):
                seg = Segment(
                    chrom, start, end, strand,
                    exon_id=attrs.get("exon_id", ""),
                    exon_number=_exon_number(attrs.get("exon_number", "")),
                    source_feature=feature,
                )
                if feature == "five_prime_utr":
                    tx_utr5[transcript_id].append(seg)
                else:
                    tx_utr3[transcript_id].append(seg)
                continue
            if not gene_id:
                # Some GFF3 files only carry Parent; keep going with a placeholder.
                gene_id = attrs.get("Parent", "").split(":")[0]
            if gene_id and transcript_id not in tx_gene_id:
                tx_gene_id[transcript_id] = gene_id
            seg = Segment(
                chrom,
                start,
                end,
                strand,
                frame=_frame(frame_s),
                exon_id=attrs.get("exon_id", ""),
                exon_number=_exon_number(attrs.get("exon_number", "")),
                source_feature=feature,
            )
            if feature == "exon":
                tx_exons[transcript_id].append(seg)
            elif feature == "CDS":
                tx_cds[transcript_id].append(seg)
            elif feature == "stop_codon":
                tx_stops[transcript_id].append(seg)
            if transcript_id not in tx_chrom:
                tx_chrom[transcript_id] = chrom
                tx_strand[transcript_id] = strand

    all_ids = set(tx_meta) | set(tx_exons) | set(tx_cds) | set(tx_stops)
    for transcript_id in all_ids:
        meta = tx_meta.get(transcript_id, {})
        chrom = str(meta.get("chrom") or tx_chrom.get(transcript_id, ""))
        strand = str(meta.get("strand") or tx_strand.get(transcript_id, "+"))
        if not chrom:
            continue
        exons = tx_exons.get(transcript_id, [])
        cds = tx_cds.get(transcript_id, [])
        stops = tx_stops.get(transcript_id, [])
        # Determine transcript span from exons when the transcript row is absent.
        start = int(meta.get("start") or (min(e.start for e in exons) if exons else 1))
        end = int(meta.get("end") or (max(e.end for e in exons) if exons else start))
        gene_id = str(meta.get("gene_id") or tx_gene_id.get(transcript_id, ""))
        tx = TranscriptModel(
            transcript_id=transcript_id.split(".")[0],
            gene_id=gene_id,
            chrom=chrom,
            strand=strand if strand in ("+", "-") else "+",
            start=start,
            end=end,
            biotype=str(meta.get("biotype") or ""),
            gene_name=str(meta.get("gene_name") or ann.gene_names.get(gene_id, "")),
            gene_biotype=str(meta.get("gene_biotype") or ann.gene_biotypes.get(gene_id, "")),
            transcript_name=str(meta.get("transcript_name") or ""),
            source=str(meta.get("source") or ""),
            tags=list(meta.get("tags") or []),  # type: ignore[arg-type]
        )
        tx.exons = _sort_exons(exons, tx.strand)
        tx.cds = _sort_cds(cds, tx.strand)
        tx.stops = sort_segments(stops, tx.strand)
        tx.utr5 = sort_segments(tx_utr5.get(transcript_id, []), tx.strand)
        tx.utr3 = sort_segments(tx_utr3.get(transcript_id, []), tx.strand)
        if tx.exons and (not tx.utr5 or not tx.utr3):
            derive_utrs(tx)
        ann.transcripts[tx.transcript_id] = tx

    if verbose:
        print(
            f"[mutatedseq] parsed {len(ann.transcripts)} transcripts, "
            f"{sum(1 for t in ann.transcripts.values() if t.cds)} with CDS "
            f"from {path}"
        )
    return ann


def _sort_exons(exons: Sequence[Segment], strand: str) -> List[Segment]:
    if not exons:
        return []
    numbered = [e for e in exons if e.exon_number is not None]
    if len(numbered) == len(exons):
        # exon_number is transcript order in Ensembl/GENCODE.
        return sorted(exons, key=lambda e: (e.exon_number, e.start))
    return sort_segments(exons, strand)


def _sort_cds(cds: Sequence[Segment], strand: str) -> List[Segment]:
    if not cds:
        return []
    numbered = [c for c in cds if c.exon_number is not None]
    if len(numbered) == len(cds):
        return sorted(cds, key=lambda c: (c.exon_number, c.start))
    return sort_segments(cds, strand)
