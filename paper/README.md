# Manuscript draft

This directory contains a bioRxiv-style preprint draft for mutatedSequence.

| File | Purpose |
|---|---|
| `manuscript.md` | Editable source of the manuscript |
| `manuscript.pdf` | Print-ready draft generated from the Markdown |
| `manuscript.docx` | Editable Word draft exported with Pandoc |
| `manuscript.html` | HTML intermediate used to build the PDF |
| `cover-letter.md` | Draft cover letter for bioRxiv |
| `SUBMISSION_CHECKLIST.md` | Items to resolve before posting |
| `figures/` | Figure 1–3 as 300 dpi PNG and vector PDF |
| `make_figures.py` | Regenerates the figures from the values reported in the manuscript |
| `build_pdf.py` | Converts the Markdown to HTML and PDF using Microsoft Edge |

## Rebuild

```bash
python paper/make_figures.py
python paper/build_pdf.py
# optional Word export
pandoc paper/manuscript.md -o paper/manuscript.docx --resource-path=paper
```

`build_pdf.py` expects Microsoft Edge at the default Windows location. If it is
not installed, open `paper/manuscript.html` in a browser and print to PDF.

## What is real and what is a placeholder

All quantitative results in the draft come from the repository:

- 1,050/1,050 chromosome 21 CDS nucleotide sequences exact against Ensembl 109;
- 952/1,050 annotated proteins exact, with 98 partial-3′/partial-5′+3′ model
  differences;
- the raw CDS length-difference distribution (+3: 831, 0: 136, +5: 40, +4: 29,
  +1: 8, +2: 6);
- 121,766 genome-wide CDS records, 99,812 terminal stops, 21,954 without a
  terminal stop, 18,559 `mRNA_end_NF`, and 758 out-of-frame stop triplets;
- the demo run summary in `examples/demo_out/run_summary.json`;
- ten end-to-end tests, verified both locally and in the Docker image.

The author names, affiliations, corresponding e-mail, ORCID iDs, funding, and
author contributions are **placeholders**. Replace them before submission. The
repository currently has no `LICENSE` file, so the manuscript says “publicly
available” rather than “open source”; adding a licence is an author decision.

The manuscript deliberately does **not** claim a head-to-head benchmark against
VEP, Haplosaurus, VCF2Prot, or ORFan. If you add such a comparison, use a shared
truth set and update the Results, Discussion, and cover letter together.
