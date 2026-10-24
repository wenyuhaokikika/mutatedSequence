"""Generate a small, self-contained synthetic dataset for mutatedSequence.

The dataset deliberately covers:
  * plus- and minus-strand genes
  * multi-exon genes and a deletion spanning an intron
  * SNV: start-lost, missense, synonymous, nonsense, stop-lost
  * 1-bp frameshift insertion/deletion and 3-bp in-frame deletion
  * 5' UTR, 3' UTR, intron and splice-donor variants
  * a 5'-partial model (mRNA_start_NF, 1 N pad)
  * a 3'-partial model with no stop codon (mRNA_end_NF)
  * a stop codon split across an exon junction
  * multiallelic and phased genotypes
  * symbolic <DEL>, <INS>, <DUP>, <INV>
  * REF mismatch, overlapping variants, 0/0, ./. and LowQual records

Run:
    python tests/make_test_data.py
"""

from __future__ import annotations

import os
import random
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from mutatedseq.sequtil import CODON_TABLE, reverse_complement, translate

DATA = os.path.join(HERE, "data")

AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
AA_CODON = {
    "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT",
    "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAA", "L": "CTG",
    "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT",
    "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT",
}
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def rand_dna(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def protein_seq(length: int) -> str:
    assert length >= 1
    return "M" + "".join(AA_ORDER[(i - 1) % 20] for i in range(1, length))


def sense_codons(n: int, offset: int = 0) -> str:
    return "".join(AA_CODON[AA_ORDER[(i + offset) % 20]] for i in range(n))


def rc(seq: str) -> str:
    return reverse_complement(seq)


@dataclass
class TestGene:
    name: str
    chrom: str
    strand: str
    exons: List[Tuple[int, int]]  # transcript order
    cds_segments: List[Tuple[int, int]]  # transcript order
    stop_segments: List[Tuple[int, int]]  # transcript order
    utr5: List[Tuple[int, int]]
    utr3: List[Tuple[int, int]]
    transcript_seq: str
    cds_seq: str  # includes 5' N pad (if any) and stop (if any)
    protein: str
    cds_core: str  # GTF CDS concatenation, without pad/stop
    pad5: int = 0
    biotype: str = "protein_coding"
    tags: List[str] = field(default_factory=list)
    cds_frames: List[int] = field(default_factory=list)

    @property
    def gene_id(self) -> str:
        return "GENE_" + self.name

    @property
    def tx_id(self) -> str:
        return "TX_" + self.name

    @property
    def start(self) -> int:
        return min(s for s, e in self.exons)

    @property
    def end(self) -> int:
        return max(e for s, e in self.exons)


GENES: List[TestGene] = []


def build_gene(
    name: str,
    chrom: str,
    strand: str,
    exons: Sequence[Tuple[int, int]],
    cds_segments: Sequence[Tuple[int, int]],
    stop_segments: Sequence[Tuple[int, int]],
    utr5: Sequence[Tuple[int, int]],
    utr3: Sequence[Tuple[int, int]],
    exon_contents: Sequence[str],
    cds_core: str,
    cds_seq: str,
    pad5: int = 0,
    biotype: str = "protein_coding",
    tags: Optional[Sequence[str]] = None,
    cds_frames: Optional[Sequence[int]] = None,
) -> TestGene:
    assert len(exons) == len(exon_contents)
    for (start, end), content in zip(exons, exon_contents):
        assert end - start + 1 == len(content), (name, start, end, len(content))
    # Transcript sequence is the concatenation of transcript-oriented exons.
    if strand == "+":
        transcript = "".join(exon_contents)
    else:
        transcript = "".join(rc(c) for c in exon_contents)
    protein = translate(cds_seq, stop_at_stop=True)[0]
    g = TestGene(
        name=name, chrom=chrom, strand=strand,
        exons=list(exons), cds_segments=list(cds_segments),
        stop_segments=list(stop_segments), utr5=list(utr5), utr3=list(utr3),
        transcript_seq=transcript, cds_seq=cds_seq, protein=protein,
        cds_core=cds_core, pad5=pad5, biotype=biotype, tags=list(tags or []),
        cds_frames=list(cds_frames or [0] * len(cds_segments)),
    )
    GENES.append(g)
    return g


def make_genes() -> None:
    # ------------------------------------------------------------------ A
    # plus strand, 3 exons, complete CDS + stop
    cds_a_core = "".join(AA_CODON[a] for a in protein_seq(233))
    utr5_a = rand_dna(100, 11)
    utr3_a = sense_codons(132, 3) + rand_dna(2, 12)
    stop = "TAA"
    exon1 = utr5_a + cds_a_core[:300]
    exon2 = cds_a_core[300:600]
    exon3 = cds_a_core[600:] + stop + utr3_a
    build_gene(
        "A", "chr1", "+",
        exons=[(501, 900), (1401, 1700), (2201, 2700)],
        cds_segments=[(601, 900), (1401, 1700), (2201, 2299)],
        stop_segments=[(2300, 2302)],
        utr5=[(501, 600)], utr3=[(2303, 2700)],
        exon_contents=[exon1, exon2, exon3],
        cds_core=cds_a_core, cds_seq=cds_a_core + stop,
    )

    # ------------------------------------------------------------------ B
    # minus strand, 2 exons, complete CDS + stop
    cds_b_core = "".join(AA_CODON[a] for a in protein_seq(199))
    utr5_b = rand_dna(100, 21)
    utr3_b = rand_dna(100, 22)
    # Exon 1 (higher coordinates) = RC(cds_exon1) + RC(utr5)
    exon_b1 = rc(cds_b_core[:300]) + rc(utr5_b)
    # Exon 2 (lower coordinates) = RC(utr3) + RC(stop) + RC(cds_exon2)
    exon_b2 = rc(utr3_b) + rc(stop) + rc(cds_b_core[300:])
    build_gene(
        "B", "chr1", "-",
        exons=[(3801, 4200), (2901, 3300)],
        cds_segments=[(3801, 4100), (3004, 3300)],
        stop_segments=[(3001, 3003)],
        utr5=[(4101, 4200)], utr3=[(2901, 3000)],
        exon_contents=[exon_b1, exon_b2],
        cds_core=cds_b_core, cds_seq=cds_b_core + stop,
    )

    # ------------------------------------------------------------------ C
    # plus strand, 5'-partial model: first CDS phase 2 -> one N pad
    cds_c_known = rand_dna(2, 31) + sense_codons(166, 4)
    utr3_c = rand_dna(97, 32)
    exon_c = cds_c_known + stop + utr3_c
    build_gene(
        "C", "chr1", "+",
        exons=[(5001, 5600)],
        cds_segments=[(5001, 5500)],
        stop_segments=[(5501, 5503)],
        utr5=[], utr3=[(5504, 5600)],
        exon_contents=[exon_c],
        cds_core=cds_c_known,
        cds_seq="N" + cds_c_known + stop,
        pad5=1,
        tags=["mRNA_start_NF"],
        cds_frames=[2],
    )

    # ------------------------------------------------------------------ D
    # plus strand, 3'-partial model: no stop codon
    cds_d_core = "".join(AA_CODON[a] for a in protein_seq(200))
    utr3_d = sense_codons(66, 5) + rand_dna(2, 41)
    exon_d = cds_d_core + utr3_d
    build_gene(
        "D", "chr1", "+",
        exons=[(6001, 6800)],
        cds_segments=[(6001, 6600)],
        stop_segments=[],
        utr5=[], utr3=[(6601, 6800)],
        exon_contents=[exon_d],
        cds_core=cds_d_core,
        cds_seq=cds_d_core,
        tags=["mRNA_end_NF"],
    )

    # ------------------------------------------------------------------ E
    # plus strand, stop codon split across an exon junction
    cds_e_core = "".join(AA_CODON[a] for a in protein_seq(100))
    utr5_e = rand_dna(100, 51)
    utr3_e = sense_codons(99, 6) + rand_dna(1, 52)
    exon_e1 = utr5_e + cds_e_core + stop[0]
    exon_e2 = stop[1:] + utr3_e
    build_gene(
        "E", "chr1", "+",
        exons=[(7001, 7401), (7601, 7900)],
        cds_segments=[(7101, 7400)],
        stop_segments=[(7401, 7401), (7601, 7602)],
        utr5=[(7001, 7100)], utr3=[(7603, 7900)],
        exon_contents=[exon_e1, exon_e2],
        cds_core=cds_e_core,
        cds_seq=cds_e_core + stop,
    )

    # ------------------------------------------------------------------ G
    # chr2, plus strand, one exon; used for multiallelic/SV tests
    cds_g_core = "".join(AA_CODON[a] for a in protein_seq(100))
    utr5_g = rand_dna(100, 61)
    utr3_g = sense_codons(66, 7) + rand_dna(2, 62)
    exon_g = utr5_g + cds_g_core + stop + utr3_g
    build_gene(
        "G", "chr2", "+",
        exons=[(501, 1103)],
        cds_segments=[(601, 900)],
        stop_segments=[(901, 903)],
        utr5=[(501, 600)], utr3=[(904, 1103)],
        exon_contents=[exon_g],
        cds_core=cds_g_core, cds_seq=cds_g_core + stop,
    )


def place(genome: Dict[str, List[str]], chrom: str, start: int, seq: str) -> None:
    start0 = start - 1
    genome[chrom][start0 : start0 + len(seq)] = list(seq)


def write_fasta(path: str, records: List[Tuple[str, str]], width: int = 60) -> None:
    with open(path, "wt", encoding="ascii", newline="\n") as handle:
        for name, seq in records:
            handle.write(">" + name + "\n")
            for i in range(0, len(seq), width):
                handle.write(seq[i : i + width] + "\n")


def write_gtf(path: str, genes: Sequence[TestGene]) -> None:
    lines = ["#!genome-build synthetic-test", "#!genome-version test"]
    for g in genes:
        attrs = (
            f'gene_id "{g.gene_id}"; gene_name "{g.name}"; '
            f'gene_biotype "{g.biotype}";'
        )
        lines.append(
            f"{g.chrom}\ttest\tgene\t{g.start}\t{g.end}\t.\t{g.strand}\t.\t{attrs}"
        )
        tag_text = "".join(f' tag "{tag}";' for tag in g.tags)
        tr_attrs = (
            f'gene_id "{g.gene_id}"; transcript_id "{g.tx_id}"; '
            f'gene_name "{g.name}"; transcript_name "{g.tx_id}"; '
            f'gene_biotype "{g.biotype}"; transcript_biotype "{g.biotype}";'
            + tag_text
        )
        lines.append(
            f"{g.chrom}\ttest\ttranscript\t{g.start}\t{g.end}\t.\t{g.strand}\t.\t{tr_attrs}"
        )
        for i, (start, end) in enumerate(g.exons, 1):
            attrs = (
                f'gene_id "{g.gene_id}"; transcript_id "{g.tx_id}"; '
                f'exon_number "{i}"; gene_name "{g.name}"; '
                f'transcript_biotype "{g.biotype}";' + tag_text
            )
            lines.append(
                f"{g.chrom}\ttest\texon\t{start}\t{end}\t.\t{g.strand}\t.\t{attrs}"
            )
        for (start, end), frame in zip(g.cds_segments, g.cds_frames):
            exon_num = next(
                i for i, (s, e) in enumerate(g.exons, 1) if start >= s and end <= e
            )
            attrs = (
                f'gene_id "{g.gene_id}"; transcript_id "{g.tx_id}"; '
                f'exon_number "{exon_num}"; gene_name "{g.name}"; '
                f'transcript_biotype "{g.biotype}";' + tag_text
            )
            lines.append(
                f"{g.chrom}\ttest\tCDS\t{start}\t{end}\t.\t{g.strand}\t{frame}\t{attrs}"
            )
        for start, end in g.stop_segments:
            exon_num = next(
                i for i, (s, e) in enumerate(g.exons, 1) if start >= s and end <= e
            )
            attrs = (
                f'gene_id "{g.gene_id}"; transcript_id "{g.tx_id}"; '
                f'exon_number "{exon_num}"; gene_name "{g.name}"; '
                f'transcript_biotype "{g.biotype}";' + tag_text
            )
            lines.append(
                f"{g.chrom}\ttest\tstop_codon\t{start}\t{end}\t.\t{g.strand}\t0\t{attrs}"
            )
        for start, end in g.utr5:
            attrs = (
                f'gene_id "{g.gene_id}"; transcript_id "{g.tx_id}"; '
                f'gene_name "{g.name}"; transcript_biotype "{g.biotype}";' + tag_text
            )
            lines.append(
                f"{g.chrom}\ttest\tfive_prime_utr\t{start}\t{end}\t.\t{g.strand}\t.\t{attrs}"
            )
        for start, end in g.utr3:
            attrs = (
                f'gene_id "{g.gene_id}"; transcript_id "{g.tx_id}"; '
                f'gene_name "{g.name}"; transcript_biotype "{g.biotype}";' + tag_text
            )
            lines.append(
                f"{g.chrom}\ttest\tthree_prime_utr\t{start}\t{end}\t.\t{g.strand}\t.\t{attrs}"
            )
    with open(path, "wt", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(lines) + "\n")


class VCFBuilder:
    def __init__(self) -> None:
        self.rows: List[str] = []

    def add(
        self,
        chrom: str,
        pos: int,
        vid: str,
        ref: str,
        alt: str,
        gt1: str = "0/1",
        gt2: str = "0/0",
        info: str = ".",
        filt: str = "PASS",
        qual: str = "50",
    ) -> None:
        self.rows.append(
            f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t{qual}\t{filt}\t{info}\t"
            f"GT:GQ:DP\t{gt1}:50:30\t{gt2}:50:30"
        )

    def write(self, path: str, contigs: Dict[str, int]) -> None:
        lines = [
            "##fileformat=VCFv4.2",
            "##source=mutatedSequence-synthetic-test",
        ]
        for name, length in contigs.items():
            lines.append(f"##contig=<ID={name},length={length}>")
        lines.extend(
            [
                '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
                '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">',
                '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="SV length">',
                '##INFO=<ID=LEFT_SVINSSEQ,Number=1,Type=String,Description="Left insertion sequence">',
                '##INFO=<ID=RIGHT_SVINSSEQ,Number=1,Type=String,Description="Right insertion sequence">',
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
            ]
        )
        with open(path, "wt", encoding="utf-8", newline="\n") as handle:
            handle.write("\n".join(lines + self.rows) + "\n")


def cds_pos(g: TestGene, idx: int, genome: Dict[str, List[str]]) -> int:
    """0-based genomic position for a GTF-CDS index in transcript order."""
    offset = 0
    for start, end in g.cds_segments:
        length = end - start + 1
        if idx < offset + length:
            rel = idx - offset
            if g.strand == "+":
                return (start - 1) + rel
            return (end - 1) - rel
        offset += length
    raise IndexError(f"{g.name}: CDS index {idx} out of range")


def other_base(base: str) -> str:
    return {"A": "C", "C": "G", "G": "T", "T": "A"}[base.upper()]


def tx_alt_to_forward(g: TestGene, alt_tx: str) -> str:
    return alt_tx if g.strand == "+" else rc(alt_tx)


def pick_missense_alt(g: TestGene, idx: int) -> str:
    """Pick a transcript-oriented base that changes the amino acid at idx."""
    cds_index = g.pad5 + idx
    codon_start = (cds_index // 3) * 3
    codon = g.cds_seq[codon_start : codon_start + 3]
    within = cds_index % 3
    ref_base = codon[within]
    aa_ref = translate(codon, stop_at_stop=False)[0]
    for alt in "ACGT":
        if alt == ref_base:
            continue
        alt_codon = codon[:within] + alt + codon[within + 1 :]
        aa_alt = translate(alt_codon, stop_at_stop=False)[0]
        if aa_alt != aa_ref and aa_alt != "*":
            return alt
    return next(b for b in "ACGT" if b != ref_base)


def add_snv_cds(
    vcf: VCFBuilder, g: TestGene, genome: Dict[str, List[str]],
    idx: int, alt_tx: str, vid: str, gt1: str = "0/1", gt2: str = "0/0",
) -> None:
    g0 = cds_pos(g, idx, genome)
    ref = genome[g.chrom][g0]
    alt = tx_alt_to_forward(g, alt_tx)
    assert alt != ref, (g.name, idx, ref, alt)
    vcf.add(g.chrom.replace("chr", ""), g0 + 1, vid, ref, alt, gt1=gt1, gt2=gt2)


def add_deletion_forward(
    vcf: VCFBuilder, g: TestGene, genome: Dict[str, List[str]],
    idx: int, n: int, vid: str, gt1: str = "0/1", gt2: str = "0/0",
) -> None:
    g0 = cds_pos(g, idx, genome)
    assert g0 >= 1
    ref = "".join(genome[g.chrom][g0 - 1 : g0 + n])
    alt = genome[g.chrom][g0 - 1]
    vcf.add(g.chrom.replace("chr", ""), g0, vid, ref, alt, gt1=gt1, gt2=gt2)


def add_insertion_forward(
    vcf: VCFBuilder, g: TestGene, genome: Dict[str, List[str]],
    idx: int, ins: str, vid: str, gt1: str = "0/1", gt2: str = "0/0",
) -> None:
    g0 = cds_pos(g, idx, genome)
    ref = genome[g.chrom][g0]
    vcf.add(g.chrom.replace("chr", ""), g0 + 1, vid, ref, ref + ins, gt1=gt1, gt2=gt2)


def add_symbolic_del(
    vcf: VCFBuilder, g: TestGene, genome: Dict[str, List[str]],
    start0: int, end0: int, vid: str, gt1: str = "0/1", gt2: str = "0/0",
) -> None:
    # Delete 0-based [start0, end0); POS is the anchor before start0.
    pos = start0
    ref = genome[g.chrom][start0 - 1]
    length = end0 - start0
    info = f"SVTYPE=DEL;END={end0};SVLEN=-{length}"
    vcf.add(g.chrom.replace("chr", ""), pos, vid, ref, "<DEL>", gt1=gt1, gt2=gt2, info=info)


def add_symbolic_ins(
    vcf: VCFBuilder, g: TestGene, genome: Dict[str, List[str]],
    g0: int, seq: str, vid: str, gt1: str = "0/1", gt2: str = "0/0",
) -> None:
    # Insert before 0-based g0; POS is the base immediately before g0.
    pos = g0
    ref = genome[g.chrom][g0 - 1]
    info = f"SVTYPE=INS;END={pos};SVLEN={len(seq)};LEFT_SVINSSEQ={seq}"
    vcf.add(g.chrom.replace("chr", ""), pos, vid, ref, "<INS>", gt1=gt1, gt2=gt2, info=info)


def add_symbolic_dup(
    vcf: VCFBuilder, g: TestGene, genome: Dict[str, List[str]],
    start0: int, end0: int, vid: str, gt1: str = "0/1", gt2: str = "0/0",
) -> None:
    pos = start0  # anchor before duplicated interval
    ref = genome[g.chrom][start0 - 1]
    length = end0 - start0
    info = f"SVTYPE=DUP;END={end0};SVLEN={length}"
    vcf.add(g.chrom.replace("chr", ""), pos, vid, ref, "<DUP>", gt1=gt1, gt2=gt2, info=info)


def add_symbolic_inv(
    vcf: VCFBuilder, g: TestGene, genome: Dict[str, List[str]],
    start0: int, end0: int, vid: str, gt1: str = "0/1", gt2: str = "0/0",
) -> None:
    pos = start0
    ref = genome[g.chrom][start0 - 1]
    length = end0 - start0
    info = f"SVTYPE=INV;END={end0};SVLEN={length}"
    vcf.add(g.chrom.replace("chr", ""), pos, vid, ref, "<INV>", gt1=gt1, gt2=gt2, info=info)


def make_vcf(genes: Dict[str, TestGene], genome: Dict[str, List[str]], path: str) -> None:
    vcf = VCFBuilder()
    A, B, C, D, E, G = (genes[x] for x in "ABCDEG")

    # Gene A (plus): coding effects.
    add_snv_cds(vcf, A, genome, 1, "A", "A_start_lost", gt1="1/1")          # ATG -> AAG
    add_snv_cds(vcf, A, genome, 3 * 10, "A", "A_missense")                  # CTG -> ATG, L->M
    add_snv_cds(vcf, A, genome, 3 * 20 + 2, "C", "A_synonymous")            # TAT -> TAC, Y->Y
    # Codon 29 is AAA (K); AAA -> TAA is a one-base nonsense change.
    add_snv_cds(vcf, A, genome, 3 * 29, "T", "A_nonsense", gt1="1/1")
    add_insertion_forward(vcf, A, genome, 3 * 40, "C", "A_frameshift_ins")
    add_deletion_forward(vcf, A, genome, 3 * 50, 3, "A_inframe_del", gt1="1/1")
    # Two-base MNV at codon 70: CTG (L) -> GAG (E).
    g0 = cds_pos(A, 3 * 70, genome)
    vcf.add("1", g0 + 1, "A_mnv", "".join(genome["chr1"][g0 : g0 + 2]), "GA")
    # Stop lost: stop codon TAA at 2300..2302, change the first base T->C.
    vcf.add("1", 2300, "A_stop_lost", genome["chr1"][2299], "C")
    # Splice donor at the first base of intron 2 (1701), and an intronic SNV
    # at 2000; both are outside the large deletion below.
    vcf.add("1", 1701, "A_splice_donor", genome["chr1"][1700], other_base(genome["chr1"][1700]))
    vcf.add("1", 550, "A_utr5", genome["chr1"][549], other_base(genome["chr1"][549]))
    vcf.add("1", 2400, "A_utr3", genome["chr1"][2399], other_base(genome["chr1"][2399]), gt1="1/1")
    vcf.add("1", 2000, "A_intron", genome["chr1"][1999], other_base(genome["chr1"][1999]))
    # A deletion spanning exon1, intron1 and exon2: 650 genomic bp, 150 CDS bp.
    add_symbolic_del(vcf, A, genome, 850, 1500, "A_del_exon_span")
    # Overlapping deletion + SNV: the deletion wins, the SNV is skipped.
    add_deletion_forward(vcf, A, genome, 3 * 60, 3, "A_overlap_del", gt1="1/1")
    g0 = cds_pos(A, 3 * 60, genome)
    snv_pos = g0 + 1  # second base of the deleted codon
    snv_ref = genome["chr1"][snv_pos]
    vcf.add("1", snv_pos + 1, "A_overlap_snv", snv_ref, other_base(snv_ref), gt1="1/1")
    # REF mismatch and genotype edge cases in intergenic space.
    bad_pos = 450
    bad_ref = genome["chr1"][bad_pos - 1]
    vcf.add("1", bad_pos, "A_bad_ref", other_base(bad_ref), "A")
    vcf.add("1", 470, "A_hom_ref", genome["chr1"][469], other_base(genome["chr1"][469]), gt1="0/0")
    vcf.add("1", 480, "A_missing_gt", genome["chr1"][479], other_base(genome["chr1"][479]), gt1="./.")
    vcf.add("1", 490, "A_lowqual", genome["chr1"][489], other_base(genome["chr1"][489]), gt1="0/1", filt="LowQual")

    # Gene B (minus): the ALT is the reverse complement of the transcript ALT.
    add_snv_cds(vcf, B, genome, 3 * 10, "A", "B_missense", gt1="1/1")       # transcript CTG -> ATG
    add_deletion_forward(vcf, B, genome, 3 * 50, 1, "B_frameshift_del")
    vcf.add("1", 3301, "B_splice", genome["chr1"][3300], other_base(genome["chr1"][3300]))
    vcf.add("1", 2950, "B_utr3", genome["chr1"][2949], other_base(genome["chr1"][2949]), gt1="1/1")
    # Minus-strand symbolic SVs: 9 bp deletion (in-frame) and 3 bp insertion.
    add_symbolic_del(
        vcf, B, genome, cds_pos(B, 330, genome),
        cds_pos(B, 330, genome) + 9, "B_symbolic_del",
    )
    add_symbolic_ins(
        vcf, B, genome, cds_pos(B, 100, genome),
        "GCT", "B_symbolic_ins",
    )

    # Gene C: 5'-partial.
    add_snv_cds(vcf, C, genome, 30, pick_missense_alt(C, 30), "C_missense")
    add_insertion_forward(vcf, C, genome, 60, "G", "C_frameshift_ins")

    # Gene D: no stop; create a stop and a frameshift. Codon 29 is AAA (K).
    add_snv_cds(vcf, D, genome, 3 * 29, "T", "D_stop_gained", gt1="1/1")
    add_insertion_forward(vcf, D, genome, 3 * 40, "T", "D_frameshift_ins")

    # Gene E: split stop codon across the exon junction.
    vcf.add("1", 7401, "E_stop_lost_split", genome["chr1"][7400], "C")

    # Gene G (chr2): multiallelic, phased, symbolic SVs.
    # Codon 94 is CAA (Q).  At its first base, ALT1=A gives AAA (K, missense)
    # and ALT2=T gives TAA (stop).  Both alleles are valid at the same REF base.
    g0 = cds_pos(G, 3 * 94, genome)
    ref = genome["chr2"][g0]
    assert ref == "C", (g0, ref)
    vcf.add("2", g0 + 1, "G_multiallelic", ref, "A,T", gt1="1/2", gt2="0/1")
    add_snv_cds(vcf, G, genome, 3 * 95, "A", "G_phased", gt1="1|0", gt2="0|1")
    # ATG -> CTG: CTG is an alternative start, not a start-lost event.
    add_snv_cds(vcf, G, genome, 0, "C", "G_alt_start")
    add_symbolic_del(vcf, G, genome, cds_pos(G, 3 * 10, genome), cds_pos(G, 3 * 10, genome) + 9, "G_symbolic_del")
    add_symbolic_ins(vcf, G, genome, cds_pos(G, 3 * 30, genome), "GCT", "G_symbolic_ins")
    add_symbolic_dup(vcf, G, genome, cds_pos(G, 3 * 50, genome), cds_pos(G, 3 * 50, genome) + 30, "G_symbolic_dup", gt1="1/1")
    add_symbolic_inv(vcf, G, genome, cds_pos(G, 3 * 70, genome), cds_pos(G, 3 * 70, genome) + 30, "G_symbolic_inv")
    # Explicit insertion in the 3' UTR.
    vcf.add("2", 1500, "G_explicit_ins", genome["chr2"][1499], genome["chr2"][1499] + "TT")
    # More edge cases on chr2.
    bad_pos = 2000
    bad_ref = genome["chr2"][bad_pos - 1]
    vcf.add("2", bad_pos, "G_bad_ref", other_base(bad_ref), "A")
    vcf.add("2", 2100, "G_lowqual", genome["chr2"][2099], other_base(genome["chr2"][2099]), filt="LowQual")
    vcf.add("2", 2200, "G_homref", genome["chr2"][2199], other_base(genome["chr2"][2199]), gt1="0/0")
    vcf.add("2", 2300, "G_missing_gt", genome["chr2"][2299], other_base(genome["chr2"][2299]), gt1="./.")

    vcf.write(path, {"1": len(genome["chr1"]), "2": len(genome["chr2"])})


def main() -> None:
    os.makedirs(DATA, exist_ok=True)
    genome = {
        "chr1": list(rand_dna(8000, 1)),
        "chr2": list(rand_dna(4000, 2)),
    }
    make_genes()
    # The exon contents are not stored; rebuild them here explicitly from the
    # same sequence rules.  (Keeping the generator deterministic.)
    a, b, c, d, e, g = (x for x in GENES)
    place(genome, "chr1", 501, a.transcript_seq[:400])
    place(genome, "chr1", 1401, a.transcript_seq[400:700])
    place(genome, "chr1", 2201, a.transcript_seq[700:])
    # Gene B: transcript order exon1 (3801..4200), exon2 (2901..3300).
    b_exon1_len = 4200 - 3801 + 1
    b_exon1 = rc(b.transcript_seq[:b_exon1_len])
    b_exon2 = rc(b.transcript_seq[b_exon1_len:])
    place(genome, "chr1", 3801, b_exon1)
    place(genome, "chr1", 2901, b_exon2)
    # Single-exon genes.
    place(genome, "chr1", 5001, c.transcript_seq)
    place(genome, "chr1", 6001, d.transcript_seq)
    # Gene E: split across two exons.
    e_exon1_len = 7401 - 7001 + 1
    place(genome, "chr1", 7001, e.transcript_seq[:e_exon1_len])
    place(genome, "chr1", 7601, e.transcript_seq[e_exon1_len:])
    place(genome, "chr2", 501, g.transcript_seq)

    genome_seq = {k: "".join(v) for k, v in genome.items()}
    write_fasta(
        os.path.join(DATA, "test.fa"),
        [(name, seq) for name, seq in genome_seq.items()],
    )
    write_gtf(os.path.join(DATA, "test.gtf"), GENES)
    write_fasta(
        os.path.join(DATA, "test.cdna.fa"),
        [(f"{x.tx_id}.1", x.transcript_seq) for x in GENES],
    )
    write_fasta(
        os.path.join(DATA, "test.cds.fa"),
        [(f"{x.tx_id}.1", x.cds_seq) for x in GENES],
    )
    write_fasta(
        os.path.join(DATA, "test.pep.fa"),
        [(f"{x.tx_id}.1", x.protein) for x in GENES],
    )
    make_vcf({x.name: x for x in GENES}, genome, os.path.join(DATA, "test.vcf"))
    print(f"Wrote test data to {DATA}")
    for x in GENES:
        print(
            f"  {x.tx_id}: {x.chrom}:{x.start}-{x.end}({x.strand}) "
            f"cDNA={len(x.transcript_seq)} CDS={len(x.cds_seq)} "
            f"protein={len(x.protein)} class_tags={x.tags}"
        )


if __name__ == "__main__":
    main()
