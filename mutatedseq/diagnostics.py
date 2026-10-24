"""Diagnostics for the two recurring questions:

1. Why does the VCF REF appear offset / why does the Ensembl CDS have extra
   leading characters?
2. Why do many transcripts not end in a stop codon?

The command is intentionally separate from the mutation pipeline so it can be
run on a full Ensembl release without producing huge sequence outputs.
"""

from __future__ import annotations

import csv
import os
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Sequence, Tuple

from .fasta import FastaReader
from .gtf import parse_gtf
from .pipeline import (
    GenomeSource,
    RunOptions,
    build_reference_bundle,
    load_reference_proteins,
    load_reference_sequences,
)
from .runner import canonicalize_annotation, select_transcripts
from .sequtil import STOP_CODONS
from .variants import check_ref
from .vcf import parse_vcf


def _common_prefix_suffix(a: str, b: str) -> Tuple[int, int]:
    p = 0
    while p < min(len(a), len(b)) and a[p] == b[p]:
        p += 1
    s = 0
    while s < min(len(a), len(b)) - p and a[-1 - s] == b[-1 - s]:
        s += 1
    return p, s


def diagnose(
    genome_path: str,
    gtf_path: str,
    outdir: str,
    cds_path: Optional[str] = None,
    cdna_path: Optional[str] = None,
    pep_path: Optional[str] = None,
    vcf_path: Optional[str] = None,
    max_transcripts: Optional[int] = 5000,
    ref_window: int = 50,
    verbose: bool = True,
) -> None:
    os.makedirs(outdir, exist_ok=True)
    if verbose:
        print("[mutatedseq] diagnostics: loading genome and annotation ...")
    genome = FastaReader(genome_path)
    ann = parse_gtf(gtf_path, verbose=verbose)
    missing = canonicalize_annotation(ann, genome)
    options = RunOptions(
        extend_mode="none",
        per_variant="off",
        all_transcripts=False,
        max_transcripts=None,
        verbose=False,
    )
    transcripts = select_transcripts(ann, options)
    if max_transcripts:
        transcripts = transcripts[:max_transcripts]
    wanted = {tx.transcript_id for tx in transcripts}
    ref_cdna = load_reference_sequences(cdna_path, wanted) if cdna_path else {}
    ref_cds = load_reference_sequences(cds_path, wanted) if cds_path else {}
    ref_pep = load_reference_proteins(pep_path, wanted) if pep_path else {}

    rows: List[Dict[str, object]] = []
    class_counts: Counter = Counter()
    cds_match_counts: Counter = Counter()
    cdna_match_counts: Counter = Counter()
    protein_match_counts: Counter = Counter()
    protein_match_by_class: Dict[str, List[int]] = defaultdict(lambda: [0, 0])
    diff_counts: Counter = Counter()
    raw_diff_counts: Counter = Counter()
    raw_prefix_counts: Counter = Counter()
    leading_n_counts: Counter = Counter()
    stop_counts: Counter = Counter()
    mod3_counts: Counter = Counter()
    examples: Dict[str, List[str]] = defaultdict(list)
    raw_examples: Dict[str, List[str]] = defaultdict(list)

    for tx in transcripts:
        bundle = build_reference_bundle(
            tx,
            genome,
            options,
            ref_cdna=ref_cdna.get(tx.transcript_id, ""),
            ref_cds=ref_cds.get(tx.transcript_id, ""),
            ref_protein=ref_pep.get(tx.transcript_id, ""),
        )
        ens = bundle.ref_cds_ensembl
        recon = bundle.ref_cds
        diff = len(ens) - len(recon) if ens else ""
        raw_len = sum(seg.length for seg in tx.cds)
        raw_seq = ""
        if ens:
            raw_seq = "".join(
                GenomeSource(genome).extract(
                    tx.chrom, seg.start0, seg.end0, tx.strand
                )
                for seg in tx.cds
            )
            raw_diff = len(ens) - len(raw_seq)
            raw_diff_counts[str(raw_diff)] += 1
            rp, rs = _common_prefix_suffix(raw_seq, ens)
            raw_prefix_counts[f"diff={raw_diff},prefix={rp},suffix={rs}"] += 1
            if raw_diff not in (0, 3) and len(raw_examples[str(raw_diff)]) < 8:
                raw_examples[str(raw_diff)].append(
                    f"{tx.transcript_id} ({tx.biotype}); raw={len(raw_seq)} "
                    f"ens={len(ens)} pad5={bundle.pad5} prefix={rp} suffix={rs} "
                    f"tags={','.join(tx.tags)}"
                )
        if ens:
            diff_counts[str(diff)] += 1
            leading_n = len(ens) - len(ens.lstrip("N"))
            leading_n_counts[str(leading_n)] += 1
            is_stop = len(ens) % 3 == 0 and ens[-3:] in STOP_CODONS
            stop_counts["stop" if is_stop else "no_stop"] += 1
            mod3_counts[str(len(ens) % 3)] += 1
            p, s = _common_prefix_suffix(recon, ens)
            if diff not in (0, 3) and len(examples[str(diff)]) < 8:
                examples[str(diff)].append(
                    f"{tx.transcript_id} ({tx.biotype}); "
                    f"recon={len(recon)} ens={len(ens)} pad5={bundle.pad5} "
                    f"prefix={p} suffix={s} tags={','.join(tx.tags)}"
                )
        class_counts[bundle.reference_class] += 1
        cds_match_counts[str(bundle.cds_matches)] += 1
        cdna_match_counts[str(bundle.cdna_matches)] += 1
        protein_match_counts[str(bundle.protein_matches)] += 1
        if bundle.protein_matches is not None:
            # [non-matching, matching]
            protein_match_by_class[bundle.reference_class][
                1 if bundle.protein_matches else 0
            ] += 1
        rows.append(
            {
                "transcript_id": tx.transcript_id,
                "gene_id": tx.gene_id,
                "gene_name": tx.gene_name,
                "chrom": tx.chrom,
                "strand": tx.strand,
                "biotype": tx.biotype,
                "tags": ",".join(tx.tags),
                "cds_len_reconstructed": len(recon),
                "cds_len_ensembl": len(ens),
                "length_diff": diff,
                "common_prefix": (p if ens else ""),
                "common_suffix": (s if ens else ""),
                "raw_gtf_cds_len": raw_len,
                "raw_length_diff": (len(ens) - len(raw_seq) if ens else ""),
                "raw_common_prefix": (rp if ens else ""),
                "raw_common_suffix": (rs if ens else ""),
                "leading_n_ensembl": (len(ens) - len(ens.lstrip("N")) if ens else ""),
                "ensembl_end_codon": (ens[-3:] if ens else ""),
                "ensembl_has_stop": (
                    (len(ens) % 3 == 0 and ens[-3:] in STOP_CODONS) if ens else ""
                ),
                "reference_class": bundle.reference_class,
                "reference_stop_reason": bundle.stop_reason,
                "cds_matches": bundle.cds_matches,
                "cdna_matches": bundle.cdna_matches,
                "protein_matches": bundle.protein_matches,
                "notes": "; ".join(bundle.notes),
            }
        )
        if verbose and len(rows) % 2000 == 0:
            print(f"[mutatedseq] diagnostics: {len(rows)} transcripts ...")

    # VCF REF validation.
    vcf_rows: List[Dict[str, object]] = []
    ref_match_counts: Counter = Counter()
    offset_counts: Counter = Counter()
    strand_counts: Counter = Counter()
    ref_examples: List[str] = []
    if vcf_path:
        records, samples, selected = parse_vcf(vcf_path, include_filtered=True)
        for rec in records:
            resolved = genome.resolve(rec.chrom)
            chrom = resolved or rec.chrom
            rc = check_ref(genome, chrom, rec.pos, rec.ref, window=ref_window)
            ref_match_counts["match" if rc.match else "mismatch"] += 1
            if not rc.match:
                offset_counts["none" if rc.offset is None else str(rc.offset)] += 1
                strand_counts[rc.offset_strand or "none"] += 1
                if len(ref_examples) < 20:
                    ref_examples.append(
                        f"{rec.chrom}:{rec.pos} REF={rec.ref} "
                        f"genome={rc.genome_ref} offset={rc.offset} "
                        f"strand={rc.offset_strand} filter={rec.filter}"
                    )
            vcf_rows.append(
                {
                    "chrom": rec.chrom,
                    "chrom_resolved": chrom,
                    "pos": rec.pos,
                    "ref": rec.ref,
                    "alt": rec.alt_string,
                    "filter": rec.filter,
                    "ref_match": rc.match,
                    "genome_ref": rc.genome_ref,
                    "ref_offset": "" if rc.offset is None else rc.offset,
                    "ref_offset_strand": rc.offset_strand,
                    "context": rc.context,
                    "note": rc.note,
                }
            )

    if rows:
        with open(
            os.path.join(outdir, "diagnostics_transcripts.tsv"),
            "wt", encoding="utf-8", newline="\n",
        ) as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(rows[0]), delimiter="\t",
                extrasaction="ignore", lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
    if vcf_rows:
        with open(
            os.path.join(outdir, "diagnostics_vcf.tsv"),
            "wt", encoding="utf-8", newline="\n",
        ) as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(vcf_rows[0]), delimiter="\t",
                extrasaction="ignore", lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(vcf_rows)

    lines = [
        "mutatedSequence diagnostics",
        "===========================",
        f"Genome: {genome_path}",
        f"GTF:    {gtf_path}",
        f"CDS:    {cds_path or '(not supplied)'}",
        f"cDNA:   {cdna_path or '(not supplied)'}",
        f"Protein:{pep_path or '(not supplied)'}",
        f"VCF:    {vcf_path or '(not supplied)'}",
        f"transcripts examined: {len(rows)}",
        f"transcripts dropped because the chromosome is absent: {len(missing)}",
        "",
        "1a. Raw GTF CDS (no stop, no N pad) versus the Ensembl CDS FASTA",
        f"   raw length difference (Ensembl - raw GTF CDS): {dict(raw_diff_counts)}",
        f"   raw common prefix/suffix patterns: {dict(raw_prefix_counts.most_common(12))}",
        f"   leading N count in Ensembl CDS: {dict(leading_n_counts)}",
        f"   Ensembl CDS terminal codon status: {dict(stop_counts)}",
        f"   Ensembl CDS length mod 3: {dict(mod3_counts)}",
    ]
    if raw_examples:
        lines.append("   examples of non-0/3 raw differences:")
        for key in sorted(raw_examples, key=lambda x: -abs(int(x))):
            lines.append(f"     diff={key}")
            for example in raw_examples[key]:
                lines.append(f"       {example}")
    lines.extend(
        [
            "",
            "   This is the direct answer to the 'extra bases at the front' question:",
            "   - A +3 difference with prefix=len(raw GTF CDS) is the stop codon at",
            "     the 3' end. Ensembl GTF CDS rows exclude the stop codon; cds.all.fa",
            "     includes it. The old createDB.py mode shows this as 'XTAA', not as",
            "     an extra genomic base at the 5' end.",
            "   - Leading Ns (1-2) are 5'-partial CDS models (mRNA_start_NF /",
            "     cds_start_NF). The Ns complete the first codon; the GTF CDS starts",
            "     at the first known base, so the raw genomic CDS is shorter at the",
            "     5' end. A +4/+5 difference is 1-2 leading Ns plus the +3 stop.",
            "   - The X in the old createDB.py mode is a placeholder created by",
            "     b.replace(a,'X'); it is not a genomic base.",
            "",
            "1b. Reconstruction validation (with stop and N pad added back)",
            f"   cDNA match counts: {dict(cdna_match_counts)}",
            f"   CDS  match counts: {dict(cds_match_counts)}",
            f"   protein match counts: {dict(protein_match_counts)}",
            "   protein match by class (nonmatch,match): "
            + str({k: tuple(v) for k, v in protein_match_by_class.items()}),
            "   (3'-incomplete models may differ at the last peptide residue; "
            "the CDS-level comparison is the authoritative check.)",
            f"   CDS length difference (Ensembl - reconstruction): {dict(diff_counts)}",
        ]
    )
    if examples:
        lines.append("   examples of non-0/3 reconstruction differences:")
        for key in sorted(examples, key=lambda x: -abs(int(x))):
            lines.append(f"     diff={key}")
            for example in examples[key]:
                lines.append(f"       {example}")
    lines.extend(
        [
            "",
            "2. Why transcripts do not end with a stop codon",
            f"   reference classes: {dict(class_counts)}",
        ]
    )
    for key, value in class_counts.most_common():
        if key != "complete":
            lines.append(f"   {key}: {value}")
    lines.extend(
        [
            "   Common explanations:",
            "   - mRNA_end_NF / cds_end_NF: 3' incomplete model; the stop codon is",
            "     outside the known transcript. This is the dominant reason in",
            "     Ensembl (about 18% of cds.all.fa records).",
            "   - mRNA_start_NF / cds_start_NF: 5' incomplete model; the CDS may be",
            "     partial and can also lack a stop.",
            "   - non_stop_decay: no stop by definition; the transcript is degraded.",
            "   - IG_*_gene / TR_*_gene V/D/J segments: not full-length proteins.",
            "   - nonsense_mediated_decay and pseudogene models may be partial.",
            "   - A terminal stop triplet in a CDS whose length is not a multiple",
            "     of 3 is not an in-frame stop.",
            "",
            "3. VCF REF validation",
        ]
    )
    if vcf_path:
        lines.extend(
            [
                f"   REF matches: {dict(ref_match_counts)}",
                f"   mismatch offsets (bases relative to POS): {dict(offset_counts)}",
                f"   offset strand: {dict(strand_counts)}",
            ]
        )
        if ref_examples:
            lines.append("   examples:")
            for example in ref_examples:
                lines.append(f"     {example}")
        lines.extend(
            [
                "   Remember:",
                "   - VCF REF/ALT are forward-strand. For a minus-strand transcript",
                "     the cDNA/CDS is reverse-complemented, so the alleles must be",
                "     reverse-complemented when applying them at transcript level.",
                "   - For an indel, POS is the anchor base; the changed bases start",
                "     at POS+1 for a pure insertion and the deleted base is POS+1",
                "     for REF=AT/ALT=A. Looking only at POS gives a 1-bp illusion.",
                "   - Soft-masked Ensembl DNA (lowercase repeats) must be compared",
                "     case-insensitively or upper-cased.",
                "   - Different assemblies/contig naming produce many mismatches, not",
                "     a consistent 1-3 bp shift. Check the offset distribution.",
                "",
            ]
        )
    else:
        lines.append("   (no VCF supplied)")
    lines.append("Per-transcript details: diagnostics_transcripts.tsv")
    if vcf_rows:
        lines.append("Per-VCF details: diagnostics_vcf.tsv")

    text = "\n".join(lines) + "\n"
    with open(os.path.join(outdir, "diagnostics.txt"), "wt", encoding="utf-8") as handle:
        handle.write(text)
    if verbose:
        print(text)
