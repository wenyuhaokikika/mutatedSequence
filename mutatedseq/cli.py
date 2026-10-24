"""Command-line interface for mutatedSequence."""

from __future__ import annotations

import argparse
import os
import sys
from typing import List, Optional

from . import __version__
from .diagnostics import diagnose
from .pipeline import RunOptions
from .runner import run


def _common_reference_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("-g", "--genome", required=True, help="genome FASTA (optionally .gz)")
    parser.add_argument("-a", "--gtf", required=True, help="Ensembl GTF/GFF3 annotation")
    parser.add_argument("--cds", help="Ensembl cds.all.fa[.gz] (recommended)")
    parser.add_argument("--cdna", help="Ensembl cdna.all.fa[.gz] (optional validation)")
    parser.add_argument("--pep", help="Ensembl pep.all.fa[.gz] (optional validation)")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mutatedseq",
        description=(
            "Generate allele-aware mutated transcript, CDS and protein sequences "
            "from a VCF and an Ensembl annotation."
        ),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    p_run = sub.add_parser("run", help="apply variants and write reports")
    _common_reference_args(p_run)
    p_run.add_argument("-v", "--vcf", required=True, help="input VCF[.gz]")
    p_run.add_argument("-o", "--out", required=True, help="output directory")
    p_run.add_argument("--sample", help="sample name (default: first VCF sample)")
    p_run.add_argument(
        "--allele-mode",
        choices=["alt", "genotype", "hom", "hap1", "hap2"],
        default="alt",
        help=(
            "which allele to apply: alt=requested ALT; genotype=first non-ref in GT; "
            "hom=only homozygous ALT; hap1/hap2=allele on a haplotype (GT order)"
        ),
    )
    p_run.add_argument(
        "--alt-index", type=int, default=1,
        help=(
            "1-based ALT index used by --allele-mode alt for multiallelic "
            "records (default: 1); biallelic records use their only ALT"
        ),
    )
    p_run.add_argument(
        "--allow-missing-gt", action="store_true",
        help="apply the ALT when GT is ./. or absent",
    )
    p_run.add_argument(
        "--include-hom-ref", action="store_true",
        help="apply ALT even for 0/0 genotypes (off by default)",
    )
    p_run.add_argument("--include-filtered", action="store_true",
                       help="do not require FILTER=PASS/.")
    p_run.add_argument(
        "--filter-values", default="",
        help="comma-separated FILTER values to accept (implies --include-filtered)",
    )
    p_run.add_argument("--min-qual", type=float, help="drop records below this QUAL")
    p_run.add_argument(
        "--overlap", choices=["skip", "error"], default="skip",
        help="policy for overlapping variants (default: skip and report)",
    )
    p_run.add_argument(
        "--strict", action="store_true",
        help=(
            "abort with exit code 2 if a parsed VCF REF does not match the "
            "genome (default: warn, report ref_match=false, and do not apply "
            "that variant)"
        ),
    )
    p_run.add_argument("--ref-window", type=int, default=50,
                       help="window for REF offset diagnostics (default: 50)")
    p_run.add_argument(
        "--extend-mode", choices=["none", "utr", "genomic", "utr+genomic"],
        default="utr",
        help=(
            "3' extension for readthrough/stop-loss translation: utr=spliced 3'UTR "
            "(default), genomic=old unspliced 5 kb behaviour, utr+genomic=both"
        ),
    )
    p_run.add_argument("--extend-bp", type=int, default=5000,
                       help="maximum 3' extension length (default: 5000)")
    p_run.add_argument(
        "--per-variant", choices=["auto", "full", "fast", "off"], default="auto",
        help=(
            "per-variant effect annotation: full=recompute each variant alone; "
            "fast=codon mapping for SNVs; off=location only. auto chooses full "
            "for <=20k records and fast otherwise"
        ),
    )
    p_run.add_argument("--all-transcripts", action="store_true",
                       help="report transcripts without overlapping variants too")
    p_run.add_argument("--include-noncoding", action="store_true",
                       help="include transcripts without a CDS")
    p_run.add_argument("--biotypes", default="",
                       help="comma-separated transcript biotypes to keep")
    p_run.add_argument("--transcripts", default="",
                       help="comma-separated transcript IDs to keep")
    p_run.add_argument("--max-transcripts", type=int,
                       help="limit the number of transcripts (testing)")
    p_run.add_argument("--no-sequences-in-tsv", action="store_true",
                       help="omit full sequences from TSV (FASTA outputs remain)")
    p_run.add_argument("--write-mutated-genome", action="store_true",
                       help="write whole-chromosome mutated FASTA (can be large)")
    p_run.add_argument("--no-normalized-vcf", action="store_true",
                       help="do not write normalized.vcf")
    p_run.add_argument("-q", "--quiet", action="store_true")

    p_diag = sub.add_parser(
        "diagnose",
        help="explain REF/CDS offsets and non-stop transcripts",
    )
    _common_reference_args(p_diag)
    p_diag.add_argument("-v", "--vcf", help="optional VCF for REF mismatch diagnostics")
    p_diag.add_argument("-o", "--out", required=True, help="output directory")
    p_diag.add_argument("--max-transcripts", type=int, default=5000,
                        help="number of transcripts to inspect (default: 5000; 0=all)")
    p_diag.add_argument("--ref-window", type=int, default=50)
    p_diag.add_argument("-q", "--quiet", action="store_true")

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "run":
        options = RunOptions(
            extend_mode=args.extend_mode,
            extend_bp=args.extend_bp,
            per_variant=args.per_variant,
            sequences_in_tsv=not args.no_sequences_in_tsv,
            all_transcripts=args.all_transcripts,
            include_noncoding=args.include_noncoding,
            only_biotypes=(
                {x.strip() for x in args.biotypes.split(",") if x.strip()}
                if args.biotypes else None
            ),
            only_transcripts=(
                {x.strip() for x in args.transcripts.split(",") if x.strip()}
                if args.transcripts else None
            ),
            max_transcripts=args.max_transcripts,
            verbose=not args.quiet,
        )
        filters = {x.strip() for x in args.filter_values.split(",") if x.strip()}
        if filters:
            args.include_filtered = True
        try:
            run(
                genome_path=args.genome,
                gtf_path=args.gtf,
                vcf_path=args.vcf,
                outdir=args.out,
                options=options,
                cds_path=args.cds,
                cdna_path=args.cdna,
                pep_path=args.pep,
                sample=args.sample,
                include_filtered=args.include_filtered,
                filter_values=filters or None,
                min_qual=args.min_qual,
                allele_mode=args.allele_mode,
                alt_index=args.alt_index,
                allow_missing_gt=args.allow_missing_gt,
                include_hom_ref=args.include_hom_ref,
                overlap_policy=args.overlap,
                ref_window=args.ref_window,
                strict_ref=args.strict,
                write_mutated_genome=args.write_mutated_genome,
                write_normalized_vcf=not args.no_normalized_vcf,
            )
        except (ValueError, FileNotFoundError, KeyError) as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 2
        return 0
    if args.command == "diagnose":
        try:
            diagnose(
                genome_path=args.genome,
                gtf_path=args.gtf,
                outdir=args.out,
                cds_path=args.cds,
                cdna_path=args.cdna,
                pep_path=args.pep,
                vcf_path=args.vcf,
                max_transcripts=(args.max_transcripts if args.max_transcripts else None),
                ref_window=args.ref_window,
                verbose=not args.quiet,
            )
        except (ValueError, FileNotFoundError, KeyError) as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 2
        return 0
    parser.error("unknown command")
    return 2
