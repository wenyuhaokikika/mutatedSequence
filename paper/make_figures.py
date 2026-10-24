#!/usr/bin/env python3
"""Generate the manuscript figures (PNG + PDF, 300 dpi)."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

OUT = Path(__file__).resolve().parent / "figures"
OUT.mkdir(parents=True, exist_ok=True)

INK = "#22333B"
MUTED = "#5E6A71"
TEAL = "#2A9D8F"
CORAL = "#D95D4F"
GOLD = "#D8A23A"
SAND = "#F6F1E7"
PALE_TEAL = "#E3F2EF"
PALE_CORAL = "#FBE7E2"
PALE_GOLD = "#FBF0D8"
GRID = "#D8D2C4"

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.edgecolor": MUTED,
        "axes.labelcolor": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "text.color": INK,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def box(ax, x, y, w, h, title, body="", fc="white", ec=GRID, lw=1.0,
        title_size=10.5, body_size=8.2, align="left"):
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.012,rounding_size=0.035",
        linewidth=lw, edgecolor=ec, facecolor=fc, zorder=2,
    )
    ax.add_patch(patch)
    tx = x + 0.022 if align == "left" else x + w / 2
    ha = "left" if align == "left" else "center"
    if body:
        ax.text(tx, y + h * 0.66, title, ha=ha, va="center",
                fontsize=title_size, fontweight="bold", color=INK, zorder=3)
        ax.text(tx, y + h * 0.30, body, ha=ha, va="center",
                fontsize=body_size, color=MUTED, linespacing=1.28, zorder=3)
    else:
        ax.text(tx, y + h / 2, title, ha=ha, va="center",
                fontsize=title_size, fontweight="bold", color=INK, zorder=3)


def arrow(ax, x1, y1, x2, y2, color=TEAL, style="-|>", lw=1.6, ls="-"):
    ax.add_patch(
        FancyArrowPatch(
            (x1, y1), (x2, y2), arrowstyle=style, mutation_scale=12,
            linewidth=lw, color=color, linestyle=ls, zorder=1,
            shrinkA=0, shrinkB=0,
        )
    )


def step(ax, x, y, n):
    ax.scatter([x], [y], s=150, color=CORAL, zorder=4, linewidths=0)
    ax.text(x, y, str(n), ha="center", va="center", color="white",
            fontsize=8.5, fontweight="bold", zorder=5)


def figure1():
    fig, ax = plt.subplots(figsize=(12, 8.6))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    main_x, main_w = 0.045, 0.62
    cx = main_x + main_w / 2

    # Step 1: inputs
    box(ax, main_x, 0.845, main_w, 0.115, "", fc=SAND, ec=GRID)
    ax.text(main_x + 0.022, 0.935, "Inputs", fontsize=10.5, fontweight="bold")
    small = [
        ("Genome FASTA", "GRCh38 / GRCm39"),
        ("GTF or GFF3", "exons, CDS, stop"),
        ("VCF", "one sample, GT"),
    ]
    sw = (main_w - 0.044 - 0.03) / 3
    for i, (t, s) in enumerate(small):
        x = main_x + 0.022 + i * (sw + 0.015)
        box(ax, x, 0.855, sw, 0.060, t, s, fc="white", title_size=8.8, body_size=7.0)
    step(ax, main_x - 0.005, 0.965, 1)

    # Step 2: reference model
    box(
        ax, main_x, 0.690, main_w, 0.105,
        "Reference transcript model + validation",
        "exons, CDS, stop_codon, UTRs; optional cDNA/CDS/pep FASTA\n"
        "5′ partial models receive an N pad; stop is added to the CDS",
        fc=PALE_TEAL, ec=TEAL,
    )
    arrow(ax, cx, 0.845, cx, 0.797)
    step(ax, main_x - 0.005, 0.742, 2)

    # Step 3: variant layer
    box(
        ax, main_x, 0.525, main_w, 0.120,
        "Variant selection and genomic edits",
        "REF check on the forward genome; GT policy (alt / genotype / hap1 / hap2)\n"
        "explicit SNV, MNV and indel edits; symbolic DEL, INS, DUP and INV",
        fc=PALE_GOLD, ec=GOLD,
    )
    arrow(ax, cx, 0.690, cx, 0.647)
    step(ax, main_x - 0.005, 0.585, 3)

    # Step 4: blocks
    box(
        ax, main_x, 0.365, main_w, 0.115,
        "Genome-first reference / ALT blocks",
        "edits are applied once per chromosome in genomic coordinates;\n"
        "conflicting overlaps are skipped and reported, never silently merged",
        fc=PALE_TEAL, ec=TEAL,
    )
    arrow(ax, cx, 0.525, cx, 0.482)
    step(ax, main_x - 0.005, 0.422, 4)

    # Step 5: extraction
    box(
        ax, main_x, 0.205, main_w, 0.115,
        "Transcript extraction by interval intersection",
        "exons are cut from the edited genome and concatenated 5′→3′;\n"
        "negative-strand transcripts are reverse-complemented as a whole",
        fc="white", ec=GRID,
    )
    arrow(ax, cx, 0.365, cx, 0.322)
    step(ax, main_x - 0.005, 0.262, 5)

    # Step 6: sequences
    box(
        ax, main_x, 0.045, main_w, 0.115,
        "Sequences and audit outputs",
        "mutated cDNA · CDS · annotated protein · 3′-extended protein\n"
        "TSV + FASTA + normalized VCF with ref_match, stop_lost and overlap flags",
        fc=PALE_CORAL, ec=CORAL,
    )
    arrow(ax, cx, 0.205, cx, 0.162)
    step(ax, main_x - 0.005, 0.102, 6)

    # Right panel
    px, pw = 0.705, 0.275
    box(ax, px, 0.515, pw, 0.445, "", fc="white", ec=GRID)
    ax.text(px + 0.022, 0.925, "Handled explicitly", fontsize=10.5,
            fontweight="bold", color=INK)
    bullets = [
        "negative-strand alleles",
        "variants spanning exon junctions",
        "5′ partial models and split stop codons",
        "multiallelic and phased genotypes",
        "symbolic DEL / INS / DUP / INV",
        "overlapping edits with a fixed policy",
        "stop-codon loss and 3′ readthrough",
    ]
    for i, text in enumerate(bullets):
        y = 0.865 - i * 0.047
        ax.add_patch(plt.Circle((px + 0.026, y), 0.0055, color=TEAL, zorder=4))
        ax.text(px + 0.040, y, text, va="center", fontsize=8.0, color=INK)

    box(ax, px, 0.045, pw, 0.435, "", fc=SAND, ec=GRID)
    ax.text(px + 0.022, 0.445, "Reported, not hidden", fontsize=10.5,
            fontweight="bold", color=INK)
    bullets2 = [
        "REF mismatch (warning or --strict)",
        "unresolved or unsupported allele",
        "overlap_skipped",
        "partial 5′ / partial 3′ model",
        "no terminal stop (stop_lost)",
        "reference CDS / protein mismatch",
        "per-variant and aggregate effects",
    ]
    for i, text in enumerate(bullets2):
        y = 0.390 - i * 0.047
        ax.add_patch(plt.Circle((px + 0.026, y), 0.0055, color=CORAL, zorder=4))
        ax.text(px + 0.040, y, text, va="center", fontsize=8.0, color=INK)

    fig.savefig(OUT / "figure1_workflow.png", dpi=300, bbox_inches="tight",
                facecolor="white")
    fig.savefig(OUT / "figure1_workflow.pdf", bbox_inches="tight",
                facecolor="white")
    plt.close(fig)


def figure2():
    fig = plt.figure(figsize=(12, 4.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.15, 1.15], wspace=0.42)

    # Panel A
    ax = fig.add_subplot(gs[0, 0])
    labels = ["CDS nucleotide\nsequence", "Annotated\nprotein"]
    values = [100.0, 952 / 1050 * 100]
    colors = [TEAL, CORAL]
    bars = ax.barh(labels, values, color=colors, height=0.48, alpha=0.9)
    for bar, val, count in zip(bars, values, ["1,050 / 1,050", "952 / 1,050"]):
        ax.text(val + 1.8, bar.get_y() + bar.get_height() / 2,
                f"{count}\n({val:.1f}%)", va="center", fontsize=8.4, color=INK)
    ax.set_xlim(0, 122)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Exact match to Ensembl 109 (%)")
    ax.set_title("chr21 reference reconstruction", fontsize=10.5,
                 fontweight="bold", loc="left")
    ax.grid(axis="x", color=GRID, linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.text(0, -0.30, "98 protein differences are confined to partial 3′ models;\n"
                      "the underlying CDS nucleotides are identical.",
            transform=ax.transAxes, fontsize=7.4, color=MUTED, va="top")

    # Panel B
    ax = fig.add_subplot(gs[0, 1])
    diffs = ["+3", "0", "+5", "+4", "+1", "+2"]
    counts = [831, 136, 40, 29, 8, 6]
    bars = ax.bar(diffs, counts, color=GOLD, width=0.62, alpha=0.95)
    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2, count + 18,
                str(count), ha="center", fontsize=8.2, color=INK)
    ax.set_ylim(0, 930)
    ax.set_xlabel("Ensembl CDS length − raw GTF CDS length (bp)")
    ax.set_ylabel("Transcripts")
    ax.set_title("Why the CDS differs by 0–5 bp", fontsize=10.5,
                 fontweight="bold", loc="left")
    ax.grid(axis="y", color=GRID, linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.text(0, -0.30, "+3 bp is the stop codon added to the CDS; +1/+2 bp are\n"
                      "5′ partial-model N pads. Both are reconstructed explicitly.",
            transform=ax.transAxes, fontsize=7.4, color=MUTED, va="top")

    # Panel C
    ax = fig.add_subplot(gs[0, 2])
    total = 121766
    stop = 99812
    nostop = total - stop
    ax.barh([""], [stop], color=TEAL, height=0.42, alpha=0.9,
            label="terminal TAA/TAG/TGA")
    ax.barh([""], [nostop], left=[stop], color=CORAL, height=0.42, alpha=0.9,
            label="no terminal stop")
    ax.text(stop / 2, 0, f"{stop:,}\n({stop/total*100:.1f}%)",
            ha="center", va="center", color="white", fontsize=8.6,
            fontweight="bold")
    ax.text(stop + nostop / 2, 0, f"{nostop:,}\n({nostop/total*100:.1f}%)",
            ha="center", va="center", color="white", fontsize=8.6,
            fontweight="bold")
    ax.set_xlim(0, total)
    ax.set_yticks([])
    ax.set_xlabel("Ensembl 109 CDS records")
    ax.set_title("A stop cannot be assumed", fontsize=10.5,
                 fontweight="bold", loc="left")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.32),
              ncol=2, frameon=False, fontsize=7.8)
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)

    fig.text(0.01, -0.10,
             "A, reference reconstruction of all 1,050 chr21 CDS transcripts.  B, distribution of CDS length differences.  "
             "C, terminal-stop status across 121,766 Ensembl 109 CDS records; 18,559 of the 21,954 no-stop records carry mRNA_end_NF.",
             fontsize=7.2, color=MUTED, ha="left")

    fig.savefig(OUT / "figure2_validation.png", dpi=300, bbox_inches="tight",
                facecolor="white")
    fig.savefig(OUT / "figure2_validation.pdf", bbox_inches="tight",
                facecolor="white")
    plt.close(fig)


def figure3():
    fig, ax = plt.subplots(figsize=(11, 5.0))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.02, 0.94, "Reference", fontsize=11, fontweight="bold", color=INK)
    ax.text(0.02, 0.54, "Stop-loss variant", fontsize=11, fontweight="bold", color=CORAL)

    # wild type mRNA
    ax.plot([0.18, 0.93], [0.82, 0.82], color=GRID, lw=8, solid_capstyle="round")
    ax.add_patch(Rectangle((0.20, 0.785), 0.47, 0.07, facecolor=PALE_TEAL,
                           edgecolor=TEAL, lw=1.0))
    ax.text(0.435, 0.82, "CDS", ha="center", va="center", fontsize=9, color=INK)
    ax.add_patch(Rectangle((0.67, 0.785), 0.075, 0.07, facecolor=CORAL,
                           edgecolor=CORAL, lw=1.0))
    ax.text(0.7075, 0.82, "TAA", ha="center", va="center", fontsize=8, color="white")
    ax.text(0.78, 0.82, "3′ UTR", ha="left", va="center", fontsize=8.5, color=MUTED)
    ax.text(0.70, 0.90, "ribosome stops", ha="center", fontsize=8.5, color=TEAL)
    arrow(ax, 0.7075, 0.90, 0.7075, 0.86, color=TEAL, lw=1.2)

    # mutant mRNA
    ax.plot([0.18, 0.93], [0.40, 0.40], color=GRID, lw=8, solid_capstyle="round")
    ax.add_patch(Rectangle((0.20, 0.365), 0.47, 0.07, facecolor=PALE_TEAL,
                           edgecolor=TEAL, lw=1.0))
    ax.text(0.435, 0.40, "CDS", ha="center", va="center", fontsize=9, color=INK)
    ax.add_patch(Rectangle((0.67, 0.365), 0.075, 0.07, facecolor=PALE_CORAL,
                           edgecolor=CORAL, lw=1.0))
    ax.text(0.7075, 0.40, "CAA", ha="center", va="center", fontsize=8, color=CORAL)
    ax.add_patch(Rectangle((0.745, 0.365), 0.185, 0.07, facecolor=PALE_GOLD,
                           edgecolor=GOLD, lw=1.0, linestyle="--"))
    ax.text(0.8375, 0.40, "3′ UTR continues", ha="center", va="center",
            fontsize=7.8, color=INK)

    # outcomes
    box(ax, 0.18, 0.12, 0.345, 0.16,
        "next in-frame stop found",
        "extended protein ends there\nprotein_stop_ext_alt = true",
        fc=PALE_TEAL, ec=TEAL, title_size=8.8, body_size=7.6)
    box(ax, 0.555, 0.12, 0.375, 0.16,
        "no stop in available sequence",
        "translation ends at the extension boundary\nstop_lost = true; no fake stop is added",
        fc=PALE_CORAL, ec=CORAL, title_size=8.8, body_size=7.6)
    arrow(ax, 0.75, 0.365, 0.40, 0.28, color=TEAL, lw=1.4)
    arrow(ax, 0.80, 0.365, 0.74, 0.28, color=CORAL, lw=1.4)

    ax.text(0.18, 0.035,
            "Default extension: spliced 3′ UTR, ≤ 5 kb.  --extend-mode utr+genomic or genomic extends further; "
            "--extend-mode none stops at the annotated CDS end.",
            fontsize=7.6, color=MUTED)

    fig.savefig(OUT / "figure3_stoploss.png", dpi=300, bbox_inches="tight",
                facecolor="white")
    fig.savefig(OUT / "figure3_stoploss.pdf", bbox_inches="tight",
                facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    figure1()
    figure2()
    figure3()
    print("wrote", *(p.name for p in sorted(OUT.iterdir())))
