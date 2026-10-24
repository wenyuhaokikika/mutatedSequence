# bioRxiv submission checklist

**Status: draft only. Do not submit until the items marked BLOCKER are resolved.**

## BLOCKER

- [ ] **Replace every author placeholder.** Add full names, affiliations, corresponding author email, and ORCID iDs. bioRxiv will not accept placeholders.
- [ ] **Add a code licence.** The repository currently has no `LICENSE` file; the manuscript therefore says “publicly available”, not “open source”. Add a licence before submission if you want to grant reuse rights. This is an authorship/legal decision, not something to infer from the code.
- [ ] **Confirm author approval and order.** All authors must approve the preprint and the author list.
- [ ] **Check the target category.** For this manuscript, “Bioinformatics” is the natural bioRxiv category.

## Scientific checks

- [ ] Confirm the chromosome 21 numbers against the current repository: 1,050/1,050 CDS exact; 952/1,050 protein exact; 91 `partial_3` and 7 `partial_5_3`.
- [ ] Confirm the genome-wide counts: 121,766 CDS records; 99,812 terminal stops; 21,954 without a terminal stop; 18,559 `mRNA_end_NF`; 758 out-of-frame stop triplets.
- [ ] Re-run `python -m unittest discover -s tests -q` and record the result in the submission notes.
- [ ] Decide whether to include a head-to-head comparison with VEP/Haplosaurus/VCF2Prot/ORFan. The current manuscript explicitly says this was not done; do not add claims that are not supported by a shared truth set.
- [ ] Verify that the demo numbers in the Results (40 records, 35 applied, 2 REF mismatches, 2 unresolved, 1 overlap, 7 symbolic, 2 stop-lost, 4 stop-gained) match `examples/demo_out/run_summary.json`.
- [ ] Check every citation against the source. In particular, VCF2Prot is cited as a 2022 bioRxiv preprint, and ORFan is cited as a GitHub resource because no peer-reviewed paper was identified at the time of writing.

## Files to upload

- [ ] Manuscript PDF or DOCX (`paper/manuscript.pdf` or `paper/manuscript.docx`).
- [ ] Figures as separate files, 300 dpi or vector: `paper/figures/figure1_workflow.png`, `figure2_validation.png`, `figure3_stoploss.png` (PDF versions are also provided).
- [ ] Tables: currently embedded in the manuscript; bioRxiv also accepts separate table files if preferred.
- [ ] Cover letter (`paper/cover-letter.md`; export to PDF or paste into the submission form).
- [ ] Supplementary material, if you add the test VCF/FASTA/GTF or a runnable example.

## Declarations

- [ ] Funding statement (or “This research received no specific grant”).
- [ ] Competing interests statement (currently “none”).
- [ ] Author contributions (currently a placeholder).
- [ ] Acknowledgements (currently a placeholder).
- [ ] Data/code availability (currently GitHub and Docker Hub; consider archiving a tagged release on Zenodo and adding the DOI).
- [ ] Ethics statement: public reference data only; no human participants or animals.
- [ ] AI-assistance disclosure if required by your institution or target journal. If used, state that an AI writing assistant helped with drafting/language editing and that the authors verified all scientific content.

## Formatting and metadata

- [ ] Title, running title, keywords.
- [ ] Abstract within the bioRxiv limit (the current draft is under 250 words).
- [ ] Figure legends match the final figures.
- [ ] Page numbers and line numbers for review (optional for a preprint).
- [ ] Author ORCID links.
- [ ] Confirm that the preprint has not been submitted elsewhere.
- [ ] Confirm that all authors agree to the preprint being posted and indexed.

## After posting

- [ ] Add the bioRxiv DOI to the GitHub README and the blog post.
- [ ] Link the preprint from the repository’s “Citation” section.
- [ ] Archive the exact code version used for the preprint (tag `v2.0.0-preprint` or similar) and add the DOI.
- [ ] Update the manuscript if the tool changes materially before journal submission.
