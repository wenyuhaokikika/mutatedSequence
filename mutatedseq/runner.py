"""End-to-end orchestration and report writing."""

from __future__ import annotations

import csv
import json
import os
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from .fasta import FastaReader, write_fasta
from .gtf import Annotation, parse_gtf
from .model import TranscriptModel
from .pipeline import (
    ReferenceBundle,
    RunOptions,
    TranscriptIndex,
    annotate_edit,
    build_reference_bundle,
    load_reference_proteins,
    load_reference_sequences,
    mutate_bundle,
    transcript_aggregate,
)
from .sequtil import reverse_complement, safe_id
from .vcf import ChosenAllele, VCFRecord, choose_allele, parse_vcf
from .variants import (
    Edit,
    MutatedRegion,
    RefCheck,
    check_ref,
    make_edit,
    resolve_edits,
)


@dataclass
class VariantResult:
    record: VCFRecord
    original_chrom: str
    chosen: Optional[ChosenAllele] = None
    edit: Optional[Edit] = None
    status: str = ""
    note: str = ""
    ref_check: Optional[RefCheck] = None


@dataclass
class RunStats:
    n_records: int = 0
    n_selected: int = 0
    n_applied: int = 0
    n_skipped_allele: int = 0
    n_ref_mismatch: int = 0
    n_unresolved: int = 0
    n_overlap: int = 0
    n_symbolic: int = 0
    n_transcripts_total: int = 0
    n_transcripts_processed: int = 0
    n_transcripts_mutated: int = 0
    n_transcripts_protein_changed: int = 0
    n_stop_lost: int = 0
    n_stop_gained: int = 0
    n_start_lost: int = 0
    n_frameshift: int = 0
    reference_classes: Counter = field(default_factory=Counter)
    consequences: Counter = field(default_factory=Counter)
    skip_reasons: Counter = field(default_factory=Counter)
    cds_matches: Counter = field(default_factory=Counter)
    cdna_matches: Counter = field(default_factory=Counter)
    protein_matches: Counter = field(default_factory=Counter)


def canonicalize_annotation(ann: Annotation, genome: FastaReader) -> List[TranscriptModel]:
    """Rename transcript chromosomes to the genome FASTA names; drop missing."""
    missing: List[TranscriptModel] = []
    for tx in list(ann.transcripts.values()):
        resolved = genome.resolve(tx.chrom)
        if resolved is None:
            missing.append(tx)
            del ann.transcripts[tx.transcript_id]
        else:
            tx.chrom = resolved
    return missing


def select_transcripts(
    ann: Annotation, options: RunOptions
) -> List[TranscriptModel]:
    out: List[TranscriptModel] = []
    wanted = options.only_transcripts
    for tx in ann.iter_transcripts():
        if not tx.cds and not getattr(options, "include_noncoding", False):
            continue
        if options.only_biotypes and tx.biotype not in options.only_biotypes:
            continue
        if wanted:
            if tx.transcript_id not in wanted and tx.transcript_id.split(".")[0] not in wanted:
                continue
        out.append(tx)
    out.sort(key=lambda t: (t.chrom, t.start0, t.transcript_id))
    if options.max_transcripts:
        out = out[: options.max_transcripts]
    return out


def _variant_base_row(vr: VariantResult) -> Dict[str, object]:
    rec = vr.record
    chosen = vr.chosen
    edit = vr.edit
    row: Dict[str, object] = {
        "chrom": rec.chrom,
        "chrom_vcf": vr.original_chrom,
        "pos": rec.pos,
        "id": rec.id,
        "ref": rec.ref,
        "alt_all": rec.alt_string,
        "alt_selected": chosen.alt if chosen else "",
        "allele_index": chosen.index if chosen else "",
        "gt": chosen.gt if chosen else "",
        "phased": chosen.phased if chosen else "",
        "qual": "" if rec.qual is None else rec.qual,
        "filter": rec.filter,
        "svtype": rec.info_value("SVTYPE") or "",
        "end": rec.info_value("END") or "",
        "svlen": rec.info_value("SVLEN") or "",
        "status": vr.status,
        "variant_note": vr.note,
        "edit_kind": edit.kind if edit else "",
        "edit_start0": edit.start0 if edit else "",
        "edit_end0": edit.end0 if edit else "",
        "edit_length_delta": edit.length_delta if edit else "",
        "symbolic": edit.symbolic if edit else "",
        "unresolved": edit.unresolved if edit else "",
        "skipped": edit.skipped if edit else "",
        "edit_note": edit.note if edit else "",
        "ref_match": "",
        "genome_ref": "",
        "ref_offset": "",
        "ref_offset_strand": "",
        "genome_context": "",
        "ref_note": "",
    }
    rc = vr.ref_check
    if rc is None and edit is not None:
        rc = edit.ref_check
    if rc is not None:
        row.update(rc.as_dict())
    return row


def _transcript_headers(tx: TranscriptModel) -> Dict[str, object]:
    return {
        "transcript_id": tx.transcript_id,
        "gene_id": tx.gene_id,
        "gene_name": tx.gene_name,
        "transcript_biotype": tx.biotype,
        "gene_biotype": tx.gene_biotype,
        "strand": tx.strand,
        "tx_start": tx.start,
        "tx_end": tx.end,
        "exon_count": len(tx.exons),
        "cds_count": len(tx.cds),
        "tags": ",".join(tx.tags),
    }


def _seq(value: str, include: bool) -> str:
    return value if include else ""


def _fasta_header(tx: TranscriptModel, bundle: ReferenceBundle, kind: str) -> str:
    return (
        f"{safe_id(tx.transcript_id)}|{safe_id(tx.gene_name or tx.gene_id)}|"
        f"{safe_id(tx.chrom)}:{tx.start}-{tx.end}({tx.strand})|{kind}|"
        f"{bundle.reference_class}"
    )


def _write_tsv(path: str, rows: Sequence[Dict[str, object]], columns: Sequence[str]) -> None:
    with open(path, "wt", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(columns), delimiter="\t",
            extrasaction="ignore", lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in columns})


def _bool(value: object) -> str:
    if value is True:
        return "true"
    if value is False:
        return "false"
    return str(value)


def _run_summary_text(stats: RunStats, missing: Sequence[TranscriptModel]) -> List[str]:
    lines = [
        "mutatedSequence run summary",
        "===========================",
        f"VCF records parsed: {stats.n_records}",
        f"alleles selected: {stats.n_selected}",
        f"edits applied: {stats.n_applied}",
        f"allele-skipped records: {stats.n_skipped_allele}",
        f"REF mismatches: {stats.n_ref_mismatch}",
        f"unresolved alleles (symbolic/REF mismatch): {stats.n_unresolved}",
        f"overlap-rejected edits: {stats.n_overlap}",
        f"symbolic edits: {stats.n_symbolic}",
        f"transcripts in annotation (CDS): {stats.n_transcripts_total}",
        f"transcripts processed: {stats.n_transcripts_processed}",
        f"transcripts with edits: {stats.n_transcripts_mutated}",
        f"transcripts with changed protein: {stats.n_transcripts_protein_changed}",
        f"stop_lost transcripts: {stats.n_stop_lost}",
        f"stop_gained transcripts: {stats.n_stop_gained}",
        f"start_lost transcripts: {stats.n_start_lost}",
        f"frameshift transcripts: {stats.n_frameshift}",
        "",
        "Reference model classes:",
    ]
    for key, value in stats.reference_classes.most_common():
        lines.append(f"  {key}: {value}")
    lines.extend(["", "Variant consequences:"])
    for key, value in stats.consequences.most_common():
        lines.append(f"  {key}: {value}")
    lines.extend(["", "Selection skip reasons:"])
    for key, value in stats.skip_reasons.most_common():
        lines.append(f"  {key}: {value}")
    lines.extend(["", "Reference validation (true/false/unknown):"])
    for name, counter in (
        ("CDS vs Ensembl", stats.cds_matches),
        ("cDNA vs Ensembl", stats.cdna_matches),
        ("protein vs Ensembl", stats.protein_matches),
    ):
        lines.append(
            f"  {name}: "
            + ", ".join(f"{k}={v}" for k, v in counter.most_common())
        )
    if missing:
        lines.extend(["", f"Transcripts dropped (chromosome absent from FASTA): {len(missing)}"])
        for tx in missing[:20]:
            lines.append(f"  {tx.transcript_id} {tx.chrom}:{tx.start}-{tx.end}")
        if len(missing) > 20:
            lines.append(f"  ... and {len(missing) - 20} more")
    return lines


def run(
    genome_path: str,
    gtf_path: str,
    vcf_path: str,
    outdir: str,
    options: RunOptions,
    cds_path: Optional[str] = None,
    cdna_path: Optional[str] = None,
    pep_path: Optional[str] = None,
    sample: Optional[str] = None,
    include_filtered: bool = False,
    filter_values: Optional[Sequence[str]] = None,
    min_qual: Optional[float] = None,
    allele_mode: str = "alt",
    alt_index: int = 1,
    allow_missing_gt: bool = False,
    include_hom_ref: bool = False,
    overlap_policy: str = "skip",
    ref_window: int = 50,
    strict_ref: bool = False,
    write_mutated_genome: bool = False,
    write_normalized_vcf: bool = True,
) -> RunStats:
    stats = RunStats()

    if options.verbose:
        print("[mutatedseq] loading genome ...")
    genome = FastaReader(genome_path)
    if options.verbose:
        print("[mutatedseq] parsing annotation ...")
    ann = parse_gtf(gtf_path, verbose=options.verbose)
    missing = canonicalize_annotation(ann, genome)
    transcripts = select_transcripts(ann, options)
    stats.n_transcripts_total = len(transcripts)
    if options.verbose:
        print(
            f"[mutatedseq] {len(transcripts)} transcripts selected "
            f"({len(missing)} dropped: chromosome missing from FASTA)"
        )

    wanted = {tx.transcript_id for tx in transcripts}
    ref_cdna = load_reference_sequences(cdna_path, wanted) if cdna_path else {}
    ref_cds = load_reference_sequences(cds_path, wanted) if cds_path else {}
    ref_pep = load_reference_proteins(pep_path, wanted) if pep_path else {}
    if options.verbose and (ref_cdna or ref_cds or ref_pep):
        print(
            f"[mutatedseq] reference FASTA loaded: cDNA={len(ref_cdna)}, "
            f"CDS={len(ref_cds)}, protein={len(ref_pep)}"
        )

    if options.verbose:
        print("[mutatedseq] parsing VCF ...")
    records, sample_names, selected_sample = parse_vcf(
        vcf_path,
        sample=sample,
        include_filtered=include_filtered,
        filter_values=filter_values,
        min_qual=min_qual,
    )
    stats.n_records = len(records)
    if options.per_variant == "auto":
        options.per_variant = "full" if len(records) <= 20000 else "fast"
    if options.verbose:
        print(
            f"[mutatedseq] per-variant effect mode: {options.per_variant} "
            f"({len(records)} records)"
        )
    if sample and sample_names and sample not in sample_names:
        raise ValueError(
            f"sample {sample!r} not in VCF; available: {', '.join(sample_names)}"
        )
    if options.verbose and sample_names:
        print(
            f"[mutatedseq] VCF samples: {', '.join(sample_names)}; "
            f"using {selected_sample!r}"
        )

    tx_index = TranscriptIndex(
        ann, transcripts={tx.transcript_id for tx in transcripts}
    )

    results: List[VariantResult] = []
    accepted_edits: List[Edit] = []
    rejected_edits: List[Edit] = []
    vr_by_edit: Dict[int, VariantResult] = {}
    ref_warn_limit = 20
    ref_warned = 0
    for record in records:
        original_chrom = record.chrom
        resolved = genome.resolve(record.chrom)
        if resolved is not None:
            record.chrom = resolved
        # Check REF once per record, before allele selection, so that even a
        # genotype-skipped record is validated.  In strict mode the first
        # mismatch aborts; otherwise the variant is reported and not applied.
        ref_check = check_ref(
            genome, record.chrom, record.pos, record.ref, window=ref_window
        )
        if not ref_check.match:
            stats.n_ref_mismatch += 1
            if strict_ref:
                raise ValueError(
                    f"REF mismatch at {record.chrom}:{record.pos}: "
                    f"VCF REF={record.ref}, genome={ref_check.genome_ref or 'N/A'}. "
                    f"{ref_check.note} Aborting because --strict is set."
                )
            if ref_warned < ref_warn_limit:
                print(
                    f"[mutatedseq] WARNING: REF mismatch at "
                    f"{record.chrom}:{record.pos}: VCF REF={record.ref}, "
                    f"genome={ref_check.genome_ref or 'N/A'}; variant not applied "
                    f"({ref_check.note}). Use --strict to abort.",
                    file=sys.stderr,
                )
            ref_warned += 1
        chosen = choose_allele(
            record,
            selected_sample,
            mode=allele_mode,
            alt_index=alt_index,
            allow_missing_gt=allow_missing_gt,
            include_hom_ref=include_hom_ref,
        )
        vr = VariantResult(
            record=record,
            original_chrom=original_chrom,
            chosen=chosen,
            ref_check=ref_check,
        )
        if chosen.skip:
            vr.status = chosen.reason
            stats.n_skipped_allele += 1
            stats.skip_reasons[chosen.reason] += 1
            results.append(vr)
            continue
        stats.n_selected += 1
        edit = make_edit(
            record, chosen, genome, ref_window=ref_window, ref_check=ref_check
        )
        edit.chrom = resolved or record.chrom
        vr.edit = edit
        vr_by_edit[id(edit)] = vr
        if edit.unresolved:
            vr.status = "unresolved"
            vr.note = edit.note or "unresolved allele"
            stats.n_unresolved += 1
            if edit.symbolic:
                stats.n_symbolic += 1
            results.append(vr)
            continue
        if edit.symbolic:
            stats.n_symbolic += 1
        vr.status = "selected"
        accepted_edits.append(edit)
        results.append(vr)

    if ref_warned > ref_warn_limit:
        print(
            f"[mutatedseq] WARNING: {ref_warned} REF mismatches; only the first "
            f"{ref_warn_limit} were printed. See variants.tsv / diagnostics.txt.",
            file=sys.stderr,
        )
    # Create output directories only after the strict REF validation has
    # passed, so a failed strict run does not leave empty directories.
    os.makedirs(outdir, exist_ok=True)
    sequences_dir = os.path.join(outdir, "sequences")
    os.makedirs(sequences_dir, exist_ok=True)

    # Resolve overlaps separately per chromosome.
    by_chrom: Dict[str, List[Edit]] = defaultdict(list)
    for edit in accepted_edits:
        by_chrom[edit.chrom].append(edit)
    resolved_edits: List[Edit] = []
    for chrom, edits in by_chrom.items():
        good, bad = resolve_edits(edits, policy=overlap_policy)
        resolved_edits.extend(good)
        rejected_edits.extend(bad)
    for edit in rejected_edits:
        stats.n_overlap += 1
        # Mark the corresponding result.
        for vr in results:
            if vr.edit is edit:
                vr.status = "overlap_skipped"
                vr.note = edit.note
                break
    stats.n_applied = len(resolved_edits)

    # Accumulate edits per transcript by querying the transcript index once per
    # edit.  This is much faster than scanning all edits for every transcript.
    edits_by_tx: Dict[str, List[Edit]] = defaultdict(list)
    tx_by_edit: Dict[str, List[TranscriptModel]] = defaultdict(list)
    for edit in resolved_edits:
        txs = tx_index.query(edit.chrom, edit.start0, edit.end0, flank=0)
        for tx in txs:
            edits_by_tx[tx.transcript_id].append(edit)
            tx_by_edit[edit.key].append(tx)

    affected_ids = set(edits_by_tx)
    if options.all_transcripts:
        processed = transcripts
    else:
        processed = [tx for tx in transcripts if tx.transcript_id in affected_ids]
    stats.n_transcripts_processed = len(processed)
    if options.verbose:
        print(
            f"[mutatedseq] {len(processed)} transcripts to report "
            f"({len(affected_ids)} with overlapping edits)"
        )

    effect_rows: List[Dict[str, object]] = []
    transcript_rows: List[Dict[str, object]] = []
    reported_edit_ids: Set[int] = set()
    fasta: Dict[str, List[Tuple[str, str]]] = defaultdict(list)
    for fasta_name in (
        "transcripts.ref.fa", "transcripts.alt.fa",
        "cds.ref.fa", "cds.alt.fa", "cds_ext.ref.fa", "cds_ext.alt.fa",
        "proteins.ref.fa", "proteins.alt.fa", "proteins.changed.fa",
        "proteins.extended.ref.fa", "proteins.extended.alt.fa",
        "proteins.extended.changed.fa",
        "genomic_regions.plus.ref.fa", "genomic_regions.plus.alt.fa",
        "genomic_regions.transcript.ref.fa",
        "genomic_regions.transcript.alt.fa",
    ):
        fasta[fasta_name] = []
    report_columns = [
        # variant
        "chrom", "chrom_vcf", "pos", "id", "ref", "alt_all", "alt_selected",
        "allele_index", "gt", "phased", "qual", "filter", "svtype", "end",
        "svlen", "status", "variant_note", "edit_kind", "edit_start0",
        "edit_end0", "edit_length_delta", "symbolic", "unresolved", "skipped",
        "edit_note", "ref_match", "genome_ref", "ref_offset",
        "ref_offset_strand", "genome_context", "ref_note",
        # transcript
        "transcript_id", "gene_id", "gene_name", "transcript_biotype",
        "gene_biotype", "strand", "tx_start", "tx_end", "exon_count",
        "cds_count", "tags", "exon_rank", "region", "splice", "consequence",
        "cds_length_delta", "cds_pos", "protein_pos", "codon_ref", "codon_alt",
        "aa_ref", "aa_alt", "stop_gained", "stop_lost", "start_lost",
        "protein_diff_single", "protein_ext_diff_single",
        "protein_diff_aggregate", "protein_ext_diff_aggregate",
        # sequences
        "transcript_ref", "transcript_alt", "cds_ref", "cds_alt",
        "cds_ext_ref", "cds_ext_alt", "protein_ref", "protein_alt",
        "protein_ext_ref", "protein_ext_alt",
        # model / status
        "reference_class", "reference_stop_reason", "transcript_status",
        "n_variants", "variant_keys", "transcript_notes",
    ]
    transcript_columns = [
        "transcript_id", "gene_id", "gene_name", "chrom", "strand", "start",
        "end", "transcript_biotype", "gene_biotype", "tags", "exon_count",
        "cds_count", "reference_class", "reference_stop_reason",
        "n_variants", "variant_keys", "consequences", "cdna_length_ref",
        "cdna_length_alt", "cds_length_ref", "cds_length_alt",
        "cds_ext_length_ref", "cds_ext_length_alt", "protein_length_ref",
        "protein_length_alt", "protein_ext_length_ref",
        "protein_ext_length_alt", "protein_stop_ref", "protein_stop_alt",
        "protein_stop_ext_ref", "protein_stop_ext_alt", "protein_changed",
        "protein_ext_changed", "stop_lost", "stop_gained", "start_lost",
        "frameshift", "protein_diff", "protein_ext_diff",
        "cdna_matches_ensembl", "cds_matches_ensembl", "protein_matches_ensembl",
        "cdna_ref", "cdna_alt", "cds_ref", "cds_alt", "cds_ext_ref",
        "cds_ext_alt", "protein_ref", "protein_alt", "protein_ext_ref",
        "protein_ext_alt", "transcript_status", "notes",
    ]

    for tx in processed:
        tx_edits = edits_by_tx.get(tx.transcript_id, [])
        bundle = build_reference_bundle(
            tx,
            genome,
            options,
            ref_cdna=ref_cdna.get(tx.transcript_id, ""),
            ref_cds=ref_cds.get(tx.transcript_id, ""),
            ref_protein=ref_pep.get(tx.transcript_id, ""),
        )
        stats.reference_classes[bundle.reference_class] += 1
        stats.cds_matches[str(bundle.cds_matches)] += 1
        stats.cdna_matches[str(bundle.cdna_matches)] += 1
        stats.protein_matches[str(bundle.protein_matches)] += 1

        if tx_edits:
            region = MutatedRegion(
                genome, tx.chrom, tx_edits, tx.start0, tx.end0
            )
            mb = mutate_bundle(bundle, region, options)
            stats.n_transcripts_mutated += 1
        else:
            # No overlapping edit: use the reference sequences directly.
            from .pipeline import MutatedBundle

            mb = MutatedBundle(
                alt_cdna=bundle.ref_cdna,
                alt_cds=bundle.ref_cds,
                alt_cds_ext=bundle.ref_cds_ext,
                alt_protein=bundle.ref_protein,
                alt_protein_annot=bundle.ref_protein_annot,
                alt_protein_ext=bundle.ref_protein_ext,
                alt_protein_stop=bundle.ref_protein_stop,
                alt_protein_ext_stop=bundle.ref_protein_ext_stop,
            )
            region = MutatedRegion(genome, tx.chrom, [], tx.start0, tx.end0)

        consequences: List[str] = []
        single_effects: List[Dict[str, object]] = []
        for edit in tx_edits:
            effect = annotate_edit(tx, bundle, edit, genome, options)
            consequences.append(str(effect.get("consequence", "")))
            single_effects.append(effect)

        aggregate = transcript_aggregate(
            bundle, mb, tx_edits, consequences, options
        )
        if aggregate["protein_changed"]:
            stats.n_transcripts_protein_changed += 1
        if aggregate["stop_lost"]:
            stats.n_stop_lost += 1
        if aggregate["stop_gained"]:
            stats.n_stop_gained += 1
        if aggregate["start_lost"]:
            stats.n_start_lost += 1
        if aggregate["frameshift"]:
            stats.n_frameshift += 1
        for consequence in consequences:
            if consequence:
                stats.consequences[consequence] += 1

        # Per-variant report rows.
        for i, (edit, effect) in enumerate(zip(tx_edits, single_effects)):
            row = _variant_base_row(vr_by_edit[id(edit)])
            reported_edit_ids.add(id(edit))
            row.update(_transcript_headers(tx))
            row.update(
                {
                    "region": effect.get("region", ""),
                    "splice": effect.get("splice", ""),
                    "consequence": effect.get("consequence", ""),
                    "cds_length_delta": effect.get("cds_length_delta", ""),
                    "cds_pos": effect.get("cds_pos", ""),
                    "protein_pos": effect.get("protein_pos", ""),
                    "codon_ref": effect.get("codon_ref", ""),
                    "codon_alt": effect.get("codon_alt", ""),
                    "aa_ref": effect.get("aa_ref", ""),
                    "aa_alt": effect.get("aa_alt", ""),
                    "stop_gained": _bool(effect.get("stop_gained", False)),
                    "stop_lost": _bool(effect.get("stop_lost", False)),
                    "start_lost": _bool(effect.get("start_lost", False)),
                    "protein_diff_single": effect.get("protein_diff", ""),
                    "protein_ext_diff_single": effect.get("protein_ext_diff", ""),
                    "protein_diff_aggregate": aggregate["protein_diff"],
                    "protein_ext_diff_aggregate": aggregate["protein_ext_diff"],
                    "transcript_ref": _seq(bundle.ref_cdna, options.sequences_in_tsv),
                    "transcript_alt": _seq(mb.alt_cdna, options.sequences_in_tsv),
                    "cds_ref": _seq(bundle.ref_cds, options.sequences_in_tsv),
                    "cds_alt": _seq(mb.alt_cds, options.sequences_in_tsv),
                    "cds_ext_ref": _seq(bundle.ref_cds_ext, options.sequences_in_tsv),
                    "cds_ext_alt": _seq(mb.alt_cds_ext, options.sequences_in_tsv),
                    "protein_ref": _seq(bundle.ref_protein, options.sequences_in_tsv),
                    "protein_alt": _seq(mb.alt_protein, options.sequences_in_tsv),
                    "protein_ext_ref": _seq(
                        bundle.ref_protein_ext, options.sequences_in_tsv
                    ),
                    "protein_ext_alt": _seq(
                        mb.alt_protein_ext, options.sequences_in_tsv
                    ),
                    "reference_class": bundle.reference_class,
                    "reference_stop_reason": bundle.stop_reason,
                    "transcript_status": "mutated",
                    "n_variants": len(tx_edits),
                    "variant_keys": ";".join(e.key for e in tx_edits),
                    "transcript_notes": "; ".join(bundle.notes + mb.notes),
                }
            )
            effect_rows.append(row)

        # Transcript-level aggregate row.
        tr = {
            **_transcript_headers(tx),
            "start": tx.start,
            "end": tx.end,
            "reference_class": bundle.reference_class,
            "reference_stop_reason": bundle.stop_reason,
            "n_variants": len(tx_edits),
            "variant_keys": ";".join(e.key for e in tx_edits),
            "consequences": aggregate["consequences"],
            "cdna_length_ref": len(bundle.ref_cdna),
            "cdna_length_alt": len(mb.alt_cdna),
            "cds_length_ref": len(bundle.ref_cds),
            "cds_length_alt": len(mb.alt_cds),
            "cds_ext_length_ref": len(bundle.ref_cds_ext),
            "cds_ext_length_alt": len(mb.alt_cds_ext),
            "protein_length_ref": aggregate["protein_length_ref"],
            "protein_length_alt": aggregate["protein_length_alt"],
            "protein_ext_length_ref": aggregate["protein_ext_length_ref"],
            "protein_ext_length_alt": aggregate["protein_ext_length_alt"],
            "protein_stop_ref": _bool(aggregate["protein_stop_ref"]),
            "protein_stop_alt": _bool(aggregate["protein_stop_alt"]),
            "protein_stop_ext_ref": _bool(aggregate["protein_stop_ext_ref"]),
            "protein_stop_ext_alt": _bool(aggregate["protein_stop_ext_alt"]),
            "protein_changed": _bool(aggregate["protein_changed"]),
            "protein_ext_changed": _bool(aggregate["protein_ext_changed"]),
            "stop_lost": _bool(aggregate["stop_lost"]),
            "stop_gained": _bool(aggregate["stop_gained"]),
            "start_lost": _bool(aggregate["start_lost"]),
            "frameshift": _bool(aggregate["frameshift"]),
            "protein_diff": aggregate["protein_diff"],
            "protein_ext_diff": aggregate["protein_ext_diff"],
            "cdna_matches_ensembl": _bool(bundle.cdna_matches),
            "cds_matches_ensembl": _bool(bundle.cds_matches),
            "protein_matches_ensembl": _bool(bundle.protein_matches),
            "cdna_ref": _seq(bundle.ref_cdna, options.sequences_in_tsv),
            "cdna_alt": _seq(mb.alt_cdna, options.sequences_in_tsv),
            "cds_ref": _seq(bundle.ref_cds, options.sequences_in_tsv),
            "cds_alt": _seq(mb.alt_cds, options.sequences_in_tsv),
            "cds_ext_ref": _seq(bundle.ref_cds_ext, options.sequences_in_tsv),
            "cds_ext_alt": _seq(mb.alt_cds_ext, options.sequences_in_tsv),
            "protein_ref": _seq(bundle.ref_protein, options.sequences_in_tsv),
            "protein_alt": _seq(mb.alt_protein, options.sequences_in_tsv),
            "protein_ext_ref": _seq(
                bundle.ref_protein_ext, options.sequences_in_tsv
            ),
            "protein_ext_alt": _seq(
                mb.alt_protein_ext, options.sequences_in_tsv
            ),
            "transcript_status": "mutated" if tx_edits else "unchanged",
            "notes": "; ".join(bundle.notes + mb.notes),
        }
        transcript_rows.append(tr)

        # FASTA outputs.
        header = _fasta_header(tx, bundle, "transcript")
        fasta["transcripts.ref.fa"].append((header, bundle.ref_cdna))
        fasta["transcripts.alt.fa"].append((header, mb.alt_cdna))
        fasta["cds.ref.fa"].append((header, bundle.ref_cds))
        fasta["cds.alt.fa"].append((header, mb.alt_cds))
        fasta["cds_ext.ref.fa"].append((header, bundle.ref_cds_ext))
        fasta["cds_ext.alt.fa"].append((header, mb.alt_cds_ext))
        fasta["proteins.ref.fa"].append((header, bundle.ref_protein))
        fasta["proteins.alt.fa"].append((header, mb.alt_protein))
        fasta["proteins.extended.ref.fa"].append(
            (header, bundle.ref_protein_ext)
        )
        fasta["proteins.extended.alt.fa"].append(
            (header, mb.alt_protein_ext)
        )
        if mb.alt_protein != bundle.ref_protein:
            fasta["proteins.changed.fa"].append((header, mb.alt_protein))
        if mb.alt_protein_ext != bundle.ref_protein_ext:
            fasta["proteins.extended.changed.fa"].append(
                (header, mb.alt_protein_ext)
            )

        # Genomic region (plus and transcript orientation).
        try:
            genomic_alt = region.fetch(tx.start0, tx.end0)
            genomic_ref = genome.fetch(tx.chrom, tx.start0, tx.end0)
            g_header = _fasta_header(tx, bundle, "genomic")
            fasta["genomic_regions.plus.ref.fa"].append((g_header, genomic_ref))
            fasta["genomic_regions.plus.alt.fa"].append((g_header, genomic_alt))
            fasta["genomic_regions.transcript.ref.fa"].append(
                (
                    g_header,
                    reverse_complement(genomic_ref)
                    if tx.strand == "-"
                    else genomic_ref,
                )
            )
            fasta["genomic_regions.transcript.alt.fa"].append(
                (
                    g_header,
                    reverse_complement(genomic_alt)
                    if tx.strand == "-"
                    else genomic_alt,
                )
            )
        except KeyError:
            pass

        if options.verbose and stats.n_transcripts_processed and len(transcript_rows) % 1000 == 0:
            print(f"[mutatedseq] processed {len(transcript_rows)} transcripts ...")

    # Variants with no overlapping transcript still get a report row.
    for vr in results:
        if vr.edit is not None and id(vr.edit) in reported_edit_ids:
            continue
        row = _variant_base_row(vr)
        if vr.status == "unresolved":
            consequence = "unresolved"
        elif vr.status == "overlap_skipped":
            consequence = "overlap_skipped"
        elif vr.edit is None:
            consequence = vr.status or "not_selected"
        else:
            consequence = "no_overlapping_transcript"
        row.update(
            {
                "transcript_status": "no_overlapping_transcript",
                "consequence": consequence,
                "transcript_notes": "variant not applied to any selected transcript",
            }
        )
        effect_rows.append(row)
        if consequence:
            stats.consequences[consequence] += 1

    # Write TSVs.
    _write_tsv(os.path.join(outdir, "mutation_report.tsv"), effect_rows, report_columns)
    _write_tsv(os.path.join(outdir, "transcripts.tsv"), transcript_rows, transcript_columns)

    variant_columns = [
        "chrom", "chrom_vcf", "pos", "id", "ref", "alt_all", "alt_selected",
        "allele_index", "gt", "phased", "qual", "filter", "svtype", "end",
        "svlen", "status", "variant_note", "edit_kind", "edit_start0",
        "edit_end0", "edit_length_delta", "symbolic", "unresolved", "skipped",
        "edit_note", "ref_match", "genome_ref", "ref_offset",
        "ref_offset_strand", "genome_context", "ref_note", "n_transcripts",
        "genes", "transcripts",
    ]
    variant_rows = []
    for vr in results:
        row = _variant_base_row(vr)
        txs = tx_by_edit.get(vr.edit.key, []) if vr.edit else []
        row["n_transcripts"] = len(txs)
        row["genes"] = ",".join(sorted({t.gene_name or t.gene_id for t in txs}))
        row["transcripts"] = ",".join(sorted({t.transcript_id for t in txs}))
        variant_rows.append(row)
    _write_tsv(os.path.join(outdir, "variants.tsv"), variant_rows, variant_columns)

    # Write FASTA files.
    for name, records in fasta.items():
        write_fasta(os.path.join(sequences_dir, name), records)

    # Optional whole-chromosome mutated FASTA.
    if write_mutated_genome:
        for chrom, edits in by_chrom.items():
            good = [e for e in resolved_edits if e.chrom == chrom]
            if not good:
                continue
            region = MutatedRegion(genome, chrom, good, 0, genome.length(chrom))
            seq = region.fetch(0, genome.length(chrom))
            write_fasta(
                os.path.join(sequences_dir, f"genome.{safe_id(chrom)}.alt.fa"),
                [(f"{chrom}|mutated_genome", seq)],
            )

    if write_normalized_vcf:
        _write_normalized_vcf(
            os.path.join(outdir, "normalized.vcf"),
            genome,
            results,
            resolved_edits,
        )

    # Summary / diagnostics.
    summary = {
        "n_records": stats.n_records,
        "n_selected": stats.n_selected,
        "n_applied": stats.n_applied,
        "n_skipped_allele": stats.n_skipped_allele,
        "n_ref_mismatch": stats.n_ref_mismatch,
        "n_unresolved": stats.n_unresolved,
        "n_overlap": stats.n_overlap,
        "n_symbolic": stats.n_symbolic,
        "n_transcripts_total": stats.n_transcripts_total,
        "n_transcripts_processed": stats.n_transcripts_processed,
        "n_transcripts_mutated": stats.n_transcripts_mutated,
        "n_transcripts_protein_changed": stats.n_transcripts_protein_changed,
        "n_stop_lost": stats.n_stop_lost,
        "n_stop_gained": stats.n_stop_gained,
        "n_start_lost": stats.n_start_lost,
        "n_frameshift": stats.n_frameshift,
        "reference_classes": dict(stats.reference_classes),
        "consequences": dict(stats.consequences),
        "skip_reasons": dict(stats.skip_reasons),
        "cds_matches": dict(stats.cds_matches),
        "cdna_matches": dict(stats.cdna_matches),
        "protein_matches": dict(stats.protein_matches),
        "sample": selected_sample,
        "allele_mode": allele_mode,
        "alt_index": alt_index,
        "overlap_policy": overlap_policy,
        "extend_mode": options.extend_mode,
        "extend_bp": options.extend_bp,
    }
    with open(os.path.join(outdir, "run_summary.json"), "wt", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    with open(os.path.join(outdir, "diagnostics.txt"), "wt", encoding="utf-8") as handle:
        handle.write("\n".join(_run_summary_text(stats, missing)) + "\n")

    if options.verbose:
        print(
            f"[mutatedseq] done: {stats.n_applied} edits applied, "
            f"{stats.n_transcripts_mutated} transcripts mutated, "
            f"{stats.n_transcripts_protein_changed} changed proteins"
        )
        print(f"[mutatedseq] reports written to {outdir}")
    return stats


def _write_normalized_vcf(
    path: str,
    genome: FastaReader,
    results: Sequence[VariantResult],
    applied: Sequence[Edit],
) -> None:
    """Write explicit REF/ALT records that bcftools consensus can consume."""
    applied_ids = {id(e) for e in applied}
    used_chroms = sorted({e.chrom for e in applied})
    with open(path, "wt", encoding="utf-8", newline="\n") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##source=mutatedSequence\n")
        handle.write(
            '##INFO=<ID=ORIG_ALT,Number=.,Type=String,Description="Original ALT alleles">\n'
        )
        handle.write(
            '##INFO=<ID=ALLELE_INDEX,Number=1,Type=Integer,Description="Selected ALT index">\n'
        )
        handle.write(
            '##INFO=<ID=GT_SELECTED,Number=1,Type=String,Description="Selected genotype">\n'
        )
        handle.write(
            '##INFO=<ID=ORIG_GT,Number=1,Type=String,Description="Original genotype">\n'
        )
        handle.write(
            '##INFO=<ID=EDIT_KIND,Number=1,Type=String,Description="mutatedSequence edit class">\n'
        )
        for chrom in used_chroms:
            try:
                length = genome.length(chrom)
                handle.write(f"##contig=<ID={chrom},length={length}>\n")
            except KeyError:
                pass
        handle.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        )
        for vr in results:
            edit = vr.edit
            if edit is None or not edit.symbolic:
                if edit is None or id(edit) not in applied_ids:
                    continue
                pos = edit.record.pos
                ref = edit.record.ref
                alt = edit.chosen.alt
                if not alt:
                    continue
            else:
                if id(edit) not in applied_ids:
                    continue
                if edit.kind == "symbolic_del":
                    pos = edit.record.pos
                    ref = genome.fetch(edit.chrom, pos - 1, edit.end0)
                    alt = ref[:1]
                elif edit.kind == "symbolic_ins":
                    pos = edit.record.pos
                    anchor = genome.fetch(edit.chrom, pos - 1, pos)
                    ref = anchor
                    alt = anchor + edit.alt
                elif edit.kind == "symbolic_dup":
                    pos = edit.end0  # insertion after the duplicated interval
                    anchor = genome.fetch(edit.chrom, pos - 1, pos)
                    ref = anchor
                    alt = anchor + edit.alt
                elif edit.kind == "symbolic_inv":
                    pos = edit.record.pos
                    anchor = genome.fetch(edit.chrom, pos - 1, pos)
                    segment = genome.fetch(edit.chrom, pos, edit.end0)
                    ref = anchor + segment
                    alt = anchor + reverse_complement(segment)
                else:
                    continue
            if not ref or not alt:
                continue
            gt = vr.chosen.gt if vr.chosen else ""
            # The normalized VCF has one ALT per record.  Remap the genotype
            # so bcftools consensus/unphased tools apply the selected ALT.
            if gt and ("/" in gt or "|" in gt):
                gt_out = "1|1" if "|" in gt else "1/1"
            elif gt:
                gt_out = "1"
            else:
                gt_out = "1/1"
            info = (
                f"ORIG_ALT={vr.record.alt_string};"
                f"ALLELE_INDEX={vr.chosen.index if vr.chosen else ''};"
                f"GT_SELECTED={gt};ORIG_GT={gt};EDIT_KIND={edit.kind}"
            )
            filt = vr.record.filter if vr.record.filter else "."
            qual = "." if vr.record.qual is None else str(vr.record.qual)
            handle.write(
                f"{edit.chrom}\t{pos}\t{vr.record.id}\t{ref}\t{alt}\t{qual}\t"
                f"{filt}\t{info}\tGT\t{gt_out}\n"
            )
