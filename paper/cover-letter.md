# Cover letter — bioRxiv preprint

[Date]

To the bioRxiv editorial team,

We are submitting the manuscript **“mutatedSequence: a genome-first, allele-aware tool for reconstructing full-length mutated transcripts, CDS, and proteins from VCF data”** for consideration as a bioRxiv preprint in the **Bioinformatics** category.

The manuscript describes an open, dependency-free tool that converts a single-sample VCF plus a genome FASTA and GTF/GFF3 annotation into complete reference and mutated cDNA, CDS, and protein sequences. The motivation is practical: variant callers and annotators describe changes in genomic coordinates, whereas structure prediction, neoantigen discovery, and proteogenomics need full-length sequences. The tool makes the conversion auditable and handles cases that commonly produce silent errors, including negative-strand genes, exon-junction variants, 5′-partial models, split stop codons, multiallelic genotypes, overlapping variants, common symbolic structural variants, and stop-codon loss with 3′ readthrough.

The manuscript reports reference reconstruction against Ensembl release 109 chromosome 21: all 1,050 CDS nucleotide sequences matched exactly, and 952 of 1,050 annotated proteins matched exactly, with the remaining differences confined to partial 3′ models. Ten end-to-end scenario tests cover the variant classes above. We also include a genome-wide annotation audit explaining why a terminal stop cannot be assumed for every transcript. We have been careful to describe the limits of the validation: this is not a population-scale accuracy estimate, and we did not perform a head-to-head benchmark against VEP, Haplosaurus, VCF2Prot, or ORFan. Those limitations are stated explicitly in the Discussion.

The source code, test data, bilingual documentation, and Docker instructions are publicly available at https://github.com/wenyuhaokikika/mutatedSequence. An environment-only Docker image is available as `wenyuhaobuemoji/mutatedsequence`. The work uses only public reference data; no human participants or animal experiments were involved.

This manuscript has not been submitted for peer review elsewhere, and all authors have approved the preprint. The authors declare no competing interests. [Add any funding or conflict-of-interest statements required by your institution.]

We suggest the preprint because the tool and its audit outputs may be immediately useful to groups converting VCFs into protein sequences for structure prediction and proteogenomics, and because we welcome community feedback before journal submission.

Sincerely,

[Corresponding author name]  
[Affiliation]  
[Email address]  
[ORCID]
