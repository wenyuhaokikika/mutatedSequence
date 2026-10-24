"""VCF parsing and genotype/allele selection.

This module intentionally implements only the parts of VCF needed for a
single-consensus sequence.  Genotype selection is explicit:

* ``alt``       apply the requested ALT allele (default: first ALT)
* ``genotype``  apply the first non-reference allele in GT
* ``hom``       apply only homozygous-alt genotypes
* ``hap1``      apply the allele on the first haplotype (GT order)
* ``hap2``      apply the allele on the second haplotype (GT order)

For unphased genotypes, haplotype assignment follows GT order; the result is a
convention, not a phased biological haplotype.  Use phased VCFs for real
haplotype-specific sequences.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple


def open_maybe_gzip(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


@dataclass
class VCFRecord:
    chrom: str
    pos: int
    id: str
    ref: str
    alts: List[str]
    qual: Optional[float]
    filter: str
    info: Dict[str, str]
    format_keys: List[str]
    samples: Dict[str, Dict[str, str]]
    line_no: int = 0
    raw_line: str = ""

    @property
    def alt_string(self) -> str:
        return ",".join(self.alts)

    def info_value(self, key: str) -> Optional[str]:
        value = self.info.get(key)
        if value is None or value is True:
            return None
        return str(value)

    def info_int(self, key: str) -> Optional[int]:
        value = self.info_value(key)
        if value is None:
            return None
        # Multi-allelic SV records can carry comma-separated END/SVLEN values;
        # use the first one for the selected (first) ALT.
        value = value.split(",")[0]
        try:
            return int(float(value))
        except ValueError:
            return None

    def info_float(self, key: str) -> Optional[float]:
        value = self.info_value(key)
        if value is None:
            return None
        try:
            return float(value)
        except ValueError:
            return None

    def gt(self, sample: Optional[str]) -> Tuple[List[Optional[int]], bool, str]:
        """Return ``(allele_indices, phased, gt_string)`` for a sample."""
        if not self.samples:
            return [], False, ""
        if sample is None:
            sample = next(iter(self.samples))
        data = self.samples.get(sample)
        if data is None or "GT" not in self.format_keys:
            return [], False, ""
        gt_text = data.get("GT", "")
        if not gt_text or gt_text in (".", "./.", ".|."):
            return [], gt_text.find("|") >= 0, gt_text
        phased = "|" in gt_text
        parts = re.split(r"[\/|]", gt_text)
        alleles: List[Optional[int]] = []
        for part in parts:
            if part == "." or part == "":
                alleles.append(None)
            else:
                try:
                    alleles.append(int(part))
                except ValueError:
                    alleles.append(None)
        return alleles, phased, gt_text


@dataclass
class ChosenAllele:
    index: int
    alt: str
    gt: str
    phased: bool
    reason: str = ""
    skip: bool = False


def parse_vcf(
    path: str,
    sample: Optional[str] = None,
    include_filtered: bool = False,
    filter_values: Optional[Sequence[str]] = None,
    min_qual: Optional[float] = None,
    chrom: Optional[str] = None,
) -> Tuple[List[VCFRecord], List[str], Optional[str]]:
    """Parse a VCF.

    Returns ``(records, sample_names, selected_sample)``.  ``selected_sample``
    is the first sample when ``sample`` is not given.
    """
    records: List[VCFRecord] = []
    sample_names: List[str] = []
    selected_sample = sample
    allowed_filters = set(filter_values or [])
    with open_maybe_gzip(path) as handle:
        for line_no, line in enumerate(handle, 1):
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n\r").split("\t")
                sample_names = fields[9:]
                if selected_sample is None and sample_names:
                    selected_sample = sample_names[0]
                continue
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) < 8:
                continue
            if chrom is not None and fields[0] != chrom:
                continue
            try:
                pos = int(fields[1])
            except ValueError:
                continue
            ref = fields[3].upper()
            alts = [a.strip().upper() for a in fields[4].split(",") if a.strip()]
            qual: Optional[float] = None
            if fields[5] not in (".", ""):
                try:
                    qual = float(fields[5])
                except ValueError:
                    qual = None
            filt = fields[6]
            if min_qual is not None and (qual is None or qual < min_qual):
                continue
            if not include_filtered:
                filt_parts = [f for f in filt.split(";") if f]
                if allowed_filters:
                    if not any(
                        f in allowed_filters or f == "." for f in filt_parts
                    ):
                        continue
                elif "PASS" not in filt_parts and filt != ".":
                    continue
            info: Dict[str, str] = {}
            if fields[7] and fields[7] != ".":
                for item in fields[7].split(";"):
                    if not item:
                        continue
                    if "=" in item:
                        key, value = item.split("=", 1)
                        info[key] = value
                    else:
                        info[item] = True  # type: ignore[assignment]
            fmt_keys: List[str] = []
            samples: Dict[str, Dict[str, str]] = {}
            if len(fields) >= 10:
                fmt_keys = fields[8].split(":")
                for name, text in zip(sample_names, fields[9:]):
                    values = text.split(":")
                    samples[name] = {
                        key: (values[i] if i < len(values) else "")
                        for i, key in enumerate(fmt_keys)
                    }
            records.append(
                VCFRecord(
                    chrom=fields[0],
                    pos=pos,
                    id=fields[2] if fields[2] != "." else f".",
                    ref=ref,
                    alts=alts,
                    qual=qual,
                    filter=filt,
                    info=info,
                    format_keys=fmt_keys,
                    samples=samples,
                    line_no=line_no,
                    raw_line=line.rstrip("\n\r"),
                )
            )
    return records, sample_names, selected_sample


def choose_allele(
    record: VCFRecord,
    sample: Optional[str],
    mode: str = "alt",
    alt_index: int = 1,
    allow_missing_gt: bool = False,
    include_hom_ref: bool = False,
) -> ChosenAllele:
    """Select one ALT allele from a record according to genotype/allele mode."""
    gt_alleles, phased, gt_text = record.gt(sample)
    if not record.alts:
        return ChosenAllele(0, "", gt_text, phased, "no_alt", skip=True)

    missing_gt = not gt_alleles or any(a is None for a in gt_alleles)
    hom_ref = bool(gt_alleles) and all(a == 0 for a in gt_alleles)

    def alt_at(index: int) -> Optional[str]:
        if 1 <= index <= len(record.alts):
            return record.alts[index - 1]
        return None

    if mode == "alt":
        if hom_ref and not include_hom_ref:
            return ChosenAllele(0, "", gt_text, phased, "hom_ref", skip=True)
        if missing_gt and not allow_missing_gt and not gt_alleles:
            # No GT at all is common in somatic/site-only VCFs; still apply the
            # requested ALT.  A missing genotype (./.) is treated as unknown.
            pass
        # --alt-index is for multiallelic records.  For a biallelic record the
        # only ALT is used, even if the user passed --alt-index 2.
        effective_index = 1 if len(record.alts) == 1 else alt_index
        alt = alt_at(effective_index)
        if alt is None:
            return ChosenAllele(
                effective_index, "", gt_text, phased,
                "alt_index_out_of_range", skip=True,
            )
        return ChosenAllele(effective_index, alt, gt_text, phased, "alt")

    if mode == "genotype":
        if missing_gt:
            if not allow_missing_gt:
                return ChosenAllele(0, "", gt_text, phased, "missing_gt", skip=True)
            alt = alt_at(alt_index)
            return ChosenAllele(alt_index, alt or "", gt_text, phased, "missing_gt_forced", skip=alt is None)
        non_ref = [a for a in gt_alleles if a not in (0, None)]
        if not non_ref:
            return ChosenAllele(0, "", gt_text, phased, "hom_ref", skip=True)
        idx = non_ref[0]
        alt = alt_at(idx)
        if alt is None:
            return ChosenAllele(idx, "", gt_text, phased, "gt_alt_out_of_range", skip=True)
        return ChosenAllele(idx, alt, gt_text, phased, "genotype")

    if mode == "hom":
        if not gt_alleles or any(a is None for a in gt_alleles):
            return ChosenAllele(0, "", gt_text, phased, "missing_gt", skip=True)
        if all(a == gt_alleles[0] for a in gt_alleles) and gt_alleles[0] not in (0, None):
            idx = gt_alleles[0]
            alt = alt_at(idx)
            if alt is None:
                return ChosenAllele(idx, "", gt_text, phased, "gt_alt_out_of_range", skip=True)
            return ChosenAllele(idx, alt, gt_text, phased, "hom")
        return ChosenAllele(0, "", gt_text, phased, "not_hom", skip=True)

    if mode in ("hap1", "hap2"):
        hap = 0 if mode == "hap1" else 1
        if len(gt_alleles) < 2:
            if not gt_alleles or gt_alleles[0] in (None, 0):
                return ChosenAllele(0, "", gt_text, phased, "missing_or_ref_hap", skip=True)
            idx = gt_alleles[0]
        else:
            idx = gt_alleles[hap]
        if idx is None:
            return ChosenAllele(0, "", gt_text, phased, "missing_allele", skip=True)
        if idx == 0:
            return ChosenAllele(0, "", gt_text, phased, "ref_hap", skip=True)
        alt = alt_at(idx)
        if alt is None:
            return ChosenAllele(idx, "", gt_text, phased, "gt_alt_out_of_range", skip=True)
        return ChosenAllele(idx, alt, gt_text, phased, mode)

    raise ValueError(
        f"unknown allele mode {mode!r}; use alt, genotype, hom, hap1 or hap2"
    )
