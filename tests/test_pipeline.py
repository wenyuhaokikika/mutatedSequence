"""End-to-end tests for mutatedSequence using the synthetic dataset."""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "tmp_out")
OUT_HAP2 = os.path.join(HERE, "tmp_out_hap2")
OUT_HAP1 = os.path.join(HERE, "tmp_out_hap1")
OUT_DIAG = os.path.join(HERE, "tmp_out_diag")
OUT_STRICT = os.path.join(HERE, "tmp_out_strict")

sys.path.insert(0, ROOT)


def run_cli(args, check=True):
    return subprocess.run(
        [sys.executable, "-m", "mutatedseq", *args],
        cwd=ROOT,
        check=check,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )


def read_tsv(path):
    with open(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path):
    out = {}
    name = None
    parts = []
    with open(path, "rt", encoding="ascii") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    out[name] = "".join(parts)
                name = line[1:].split("|")[0].split(".")[0]
                parts = []
            else:
                parts.append(line.strip())
        if name is not None:
            out[name] = "".join(parts)
    return out


class PipelineTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not os.path.exists(os.path.join(DATA, "test.vcf")):
            subprocess.run(
                [sys.executable, os.path.join(HERE, "make_test_data.py")],
                cwd=ROOT, check=True,
            )
        for d in (OUT, OUT_HAP2, OUT_HAP1, OUT_DIAG, OUT_STRICT):
            shutil.rmtree(d, ignore_errors=True)
        common = [
            "-g", os.path.join(DATA, "test.fa"),
            "-a", os.path.join(DATA, "test.gtf"),
            "--cds", os.path.join(DATA, "test.cds.fa"),
            "--cdna", os.path.join(DATA, "test.cdna.fa"),
            "--pep", os.path.join(DATA, "test.pep.fa"),
        ]
        vcf = os.path.join(DATA, "test.vcf")
        cls.common = common
        cls.default_run = run_cli(
            ["run", *common, "-v", vcf, "-o", OUT, "--sample", "S1"]
        )
        run_cli([
            "run", *common, "-v", vcf, "-o", OUT_HAP2,
            "--sample", "S1", "--allele-mode", "hap2",
        ])
        run_cli([
            "run", *common, "-v", vcf, "-o", OUT_HAP1,
            "--sample", "S1", "--allele-mode", "hap1",
        ])
        run_cli([
            "diagnose", *common, "-v", vcf, "-o", OUT_DIAG,
            "--max-transcripts", "0",
        ])

        cls.report = read_tsv(os.path.join(OUT, "mutation_report.tsv"))
        cls.transcripts = read_tsv(os.path.join(OUT, "transcripts.tsv"))
        cls.variants = read_tsv(os.path.join(OUT, "variants.tsv"))
        cls.by_id = {}
        for row in cls.report:
            cls.by_id.setdefault(row["id"], []).append(row)
        cls.tx = {row["transcript_id"]: row for row in cls.transcripts}

    def effect(self, vid, tid=None):
        rows = self.by_id[vid]
        if tid is not None:
            rows = [r for r in rows if r["transcript_id"] == tid]
        self.assertEqual(len(rows), 1, (vid, len(rows)))
        return rows[0]

    def test_reference_sequences_match_design(self):
        ref_tx = read_fasta(os.path.join(OUT, "sequences", "transcripts.ref.fa"))
        ref_cds = read_fasta(os.path.join(OUT, "sequences", "cds.ref.fa"))
        ref_prot = read_fasta(os.path.join(OUT, "sequences", "proteins.ref.fa"))
        # The test data has 6 transcripts with deterministic lengths.
        expected = {
            "TX_A": (1200, 702, 233),
            "TX_B": (800, 600, 199),
            "TX_C": (600, 504, 167),
            "TX_D": (800, 600, 200),
            "TX_E": (701, 303, 100),
            "TX_G": (603, 303, 100),
        }
        design_tx = read_fasta(os.path.join(DATA, "test.cdna.fa"))
        design_cds = read_fasta(os.path.join(DATA, "test.cds.fa"))
        design_prot = read_fasta(os.path.join(DATA, "test.pep.fa"))
        for tid, (cdna_len, cds_len, prot_len) in expected.items():
            self.assertEqual(len(ref_tx[tid]), cdna_len, tid)
            self.assertEqual(len(ref_cds[tid]), cds_len, tid)
            self.assertEqual(len(ref_prot[tid]), prot_len, tid)
            self.assertEqual(ref_tx[tid], design_tx[tid], tid)
            self.assertEqual(ref_cds[tid], design_cds[tid], tid)
            self.assertEqual(ref_prot[tid], design_prot[tid], tid)
        self.assertTrue(ref_prot["TX_C"].startswith("X"))
        self.assertEqual(self.tx["TX_C"]["reference_class"], "partial_5_complete_3")
        self.assertEqual(self.tx["TX_D"]["reference_class"], "partial_3")
        self.assertEqual(self.tx["TX_D"]["protein_stop_ref"], "false")

    def test_snv_effects(self):
        nonsense = self.effect("A_nonsense")
        self.assertEqual(nonsense["consequence"], "nonsense")
        self.assertEqual(nonsense["aa_ref"], "K")
        self.assertEqual(nonsense["aa_alt"], "*")
        self.assertEqual(nonsense["stop_gained"], "true")
        syn = self.effect("A_synonymous")
        self.assertEqual(syn["consequence"], "synonymous")
        self.assertEqual(syn["protein_diff_single"], "")
        start = self.effect("A_start_lost")
        self.assertEqual(start["consequence"], "start_lost")
        self.assertEqual(start["start_lost"], "true")
        stop = self.effect("A_stop_lost")
        self.assertEqual(stop["consequence"], "stop_lost")
        self.assertEqual(stop["stop_lost"], "true")
        self.assertEqual(stop["codon_ref"], "TAA")
        missense = self.effect("A_missense")
        self.assertEqual(missense["consequence"], "missense")
        self.assertEqual(missense["aa_ref"], "L")
        self.assertEqual(missense["aa_alt"], "M")

    def test_frameshift_and_inframe(self):
        fs = self.effect("A_frameshift_ins")
        self.assertEqual(fs["consequence"], "frameshift")
        self.assertEqual(fs["cds_length_delta"], "1")
        inframe = self.effect("A_inframe_del")
        self.assertEqual(inframe["consequence"], "inframe_deletion")
        self.assertEqual(inframe["cds_length_delta"], "-3")
        # A deletion spanning an intron removes 150 CDS bases: in-frame.
        span = self.effect("A_del_exon_span")
        self.assertEqual(span["consequence"], "inframe_deletion")
        self.assertEqual(span["cds_length_delta"], "-150")

    def test_negative_strand_alleles(self):
        # VCF REF/ALT are forward-strand; the transcript codon is CTG -> ATG.
        missense = self.effect("B_missense")
        self.assertEqual(missense["consequence"], "missense")
        self.assertEqual(missense["codon_ref"], "CTG")
        self.assertEqual(missense["codon_alt"], "ATG")
        self.assertEqual(missense["aa_ref"], "L")
        self.assertEqual(missense["aa_alt"], "M")
        fs = self.effect("B_frameshift_del")
        self.assertEqual(fs["consequence"], "frameshift")
        self.assertEqual(fs["cds_length_delta"], "-1")
        splice = self.effect("B_splice")
        self.assertEqual(splice["consequence"], "splice_acceptor")
        mnv = self.effect("A_mnv")
        self.assertEqual(mnv["consequence"], "missense")
        self.assertEqual(mnv["codon_ref"], "CTG")
        self.assertEqual(mnv["codon_alt"], "GAG")
        self.assertEqual(mnv["aa_ref"], "L")
        self.assertEqual(mnv["aa_alt"], "E")
        # Minus-strand symbolic SVs must also be applied in transcript orientation.
        bdel = self.effect("B_symbolic_del")
        self.assertEqual(bdel["consequence"], "inframe_deletion")
        self.assertEqual(bdel["cds_length_delta"], "-9")
        bins = self.effect("B_symbolic_ins")
        self.assertEqual(bins["consequence"], "inframe_insertion")
        self.assertEqual(bins["cds_length_delta"], "3")

    def test_partial_models_and_split_stop(self):
        c_mis = self.effect("C_missense")
        self.assertIn(c_mis["consequence"], ("missense", "synonymous"))
        c_fs = self.effect("C_frameshift_ins")
        self.assertEqual(c_fs["consequence"], "frameshift")
        d_stop = self.effect("D_stop_gained")
        self.assertEqual(d_stop["consequence"], "nonsense")
        self.assertEqual(d_stop["stop_gained"], "true")
        d_fs = self.effect("D_frameshift_ins")
        self.assertEqual(d_fs["consequence"], "frameshift")
        e_stop = self.effect("E_stop_lost_split")
        self.assertEqual(e_stop["consequence"], "stop_lost")
        self.assertEqual(e_stop["stop_lost"], "true")
        self.assertGreater(
            int(self.tx["TX_E"]["protein_length_alt"]),
            int(self.tx["TX_E"]["protein_length_ref"]),
        )

    def test_multiallelic_and_phased(self):
        alt1 = self.effect("G_multiallelic")
        self.assertEqual(alt1["alt_selected"], "A")
        self.assertEqual(alt1["codon_ref"], "CAA")
        self.assertEqual(alt1["codon_alt"], "AAA")
        self.assertEqual(alt1["consequence"], "missense")
        # hap2 picks ALT2=T -> TAA stop.
        hap2_rows = read_tsv(os.path.join(OUT_HAP2, "mutation_report.tsv"))
        row = next(
            r for r in hap2_rows
            if r["id"] == "G_multiallelic" and r["transcript_id"] == "TX_G"
        )
        self.assertEqual(row["alt_selected"], "T")
        self.assertEqual(row["consequence"], "nonsense")
        self.assertEqual(row["aa_alt"], "*")
        alt_start = self.effect("G_alt_start")
        self.assertEqual(alt_start["consequence"], "missense")
        self.assertEqual(alt_start["start_lost"], "false")
        self.assertEqual(alt_start["codon_ref"], "ATG")
        self.assertEqual(alt_start["codon_alt"], "CTG")
        # G_phased is 1|0: hap1 applies, hap2 does not.
        hap1_rows = read_tsv(os.path.join(OUT_HAP1, "mutation_report.tsv"))
        self.assertTrue(
            any(r["id"] == "G_phased" and r["transcript_id"] == "TX_G" for r in hap1_rows)
        )
        self.assertFalse(
            any(
                r["id"] == "G_phased" and r["transcript_id"] == "TX_G"
                for r in hap2_rows
            )
        )

    def test_symbolic_variants(self):
        cases = {
            "G_symbolic_del": ("inframe_deletion", "-9"),
            "G_symbolic_ins": ("inframe_insertion", "3"),
            "G_symbolic_dup": ("inframe_insertion", "30"),
            "B_symbolic_del": ("inframe_deletion", "-9"),
            "B_symbolic_ins": ("inframe_insertion", "3"),
        }
        for vid, (consequence, delta) in cases.items():
            row = self.effect(vid)
            self.assertEqual(row["consequence"], consequence, vid)
            self.assertEqual(row["cds_length_delta"], delta, vid)
        inv = self.effect("G_symbolic_inv")
        self.assertIn(inv["consequence"], ("protein_altering", "inframe_deletion", "missense"))

    def test_ref_mismatch_overlap_and_genotypes(self):
        # Default behaviour: warning + ref_match=false, not a fatal error.
        self.assertIn("WARNING: REF mismatch", self.default_run.stdout)
        bad = next(r for r in self.variants if r["id"] == "A_bad_ref")
        self.assertEqual(bad["status"], "unresolved")
        self.assertEqual(bad["ref_match"], "False")
        overlap = next(r for r in self.variants if r["id"] == "A_overlap_snv")
        self.assertEqual(overlap["status"], "overlap_skipped")
        hom = next(r for r in self.variants if r["id"] == "A_hom_ref")
        self.assertEqual(hom["status"], "hom_ref")
        # The overlapping deletion itself is applied.
        row = self.effect("A_overlap_del")
        self.assertEqual(row["consequence"], "inframe_deletion")
        # --strict must abort on the first REF mismatch with exit code 2.
        strict = run_cli(
            ["run", *self.common, "-v", os.path.join(DATA, "test.vcf"),
             "-o", OUT_STRICT, "--sample", "S1", "--strict"],
            check=False,
        )
        self.assertEqual(strict.returncode, 2)
        self.assertIn("REF mismatch", strict.stdout)
        self.assertIn("--strict", strict.stdout)

    def test_diagnostics_raw_offsets(self):
        with open(os.path.join(OUT_DIAG, "diagnostics.txt"), "rt", encoding="utf-8") as handle:
            text = handle.read()
        self.assertIn("raw length difference", text)
        self.assertIn("'3'", text)
        self.assertIn("'4'", text)
        self.assertIn("leading N count", text)
        self.assertIn("mismatch offsets", text)
        diag = read_tsv(os.path.join(OUT_DIAG, "diagnostics_transcripts.tsv"))
        by_id = {r["transcript_id"]: r for r in diag}
        self.assertEqual(by_id["TX_A"]["raw_length_diff"], "3")
        self.assertEqual(by_id["TX_C"]["raw_length_diff"], "4")
        self.assertEqual(by_id["TX_C"]["leading_n_ensembl"], "1")
        self.assertEqual(by_id["TX_D"]["raw_length_diff"], "0")

    def test_outputs_exist(self):
        for name in (
            "mutation_report.tsv", "transcripts.tsv", "variants.tsv",
            "normalized.vcf", "run_summary.json", "diagnostics.txt",
        ):
            self.assertTrue(os.path.exists(os.path.join(OUT, name)), name)
        for name in (
            "transcripts.alt.fa", "cds.alt.fa", "proteins.alt.fa",
            "proteins.changed.fa", "genomic_regions.plus.alt.fa",
            "genomic_regions.transcript.alt.fa",
        ):
            self.assertTrue(
                os.path.exists(os.path.join(OUT, "sequences", name)), name
            )
        # normalized.vcf should not contain applied symbolic ALT alleles.
        with open(os.path.join(OUT, "normalized.vcf"), "rt", encoding="utf-8") as handle:
            text = handle.read()
        self.assertIn("G_symbolic_del", text)
        self.assertNotIn("\t<DEL>", text.split("#CHROM")[1])


if __name__ == "__main__":
    unittest.main(verbosity=2)
