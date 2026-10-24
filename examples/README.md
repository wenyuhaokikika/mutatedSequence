# Example output

[English](../README.md) | [简体中文](../README.zh-CN.md)

These files were generated from the synthetic test dataset:

```bash
python tests/make_test_data.py
python -m mutatedseq run \
  -g tests/data/test.fa \
  -a tests/data/test.gtf \
  -v tests/data/test.vcf \
  -o examples/demo_out \
  --cds tests/data/test.cds.fa \
  --cdna tests/data/test.cdna.fa \
  --pep tests/data/test.pep.fa \
  --sample S1
```

- `mutation_report.tsv` — one row per variant × transcript, including original
  VCF REF/ALT, selected allele, REF check, codon/AA change, and full
  before/after cDNA, CDS and protein sequences.
- `transcripts.tsv` — one row per transcript with the aggregate effect.
- `variants.tsv` — one row per VCF record / selected allele.
- `run_summary.json` / `diagnostics.txt` — run statistics.

The sequences are synthetic and are not biological sequences. The purpose is
to make every column inspectable and to demonstrate the different scenarios in
[the main README](../README.md#test-scenarios).
