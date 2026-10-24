#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""DEPRECATED legacy concatenation/translation script.

The current implementation reconstructs transcripts and proteins from the
mutated genome, handles negative-strand alleles and variants spanning exon
junctions, and writes a per-variant report.  Use:

    python -m mutatedseq run ...

The original code is kept in git history.  This stub exits immediately so that
an old command cannot silently produce incorrect proteins.
"""

import sys

sys.exit(
    "DEPRECATED: use 'python -m mutatedseq run ...' (see README.md) instead of "
    "mutatedProteins.py."
)
