"""FASTA input with optional .fai random access.

A human genome should normally be supplied as an uncompressed FASTA plus a
samtools .fai index.  Gzipped FASTA is accepted for small test files and is
loaded into memory; for a 3 GB genome, decompress it once first.
"""

from __future__ import annotations

import gzip
import io
import os
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Set, Tuple


def open_text(path: str) -> io.TextIOBase:
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="ascii", errors="replace")
    return open(path, "rt", encoding="ascii", errors="replace")


def read_fasta(
    path: str,
    wanted: Optional[Set[str]] = None,
    uppercase: bool = True,
) -> Dict[str, str]:
    """Read a (possibly gzipped) FASTA into a dict.

    ``wanted`` limits memory use when only a subset of transcript IDs is needed.
    """
    out: Dict[str, str] = {}
    name: Optional[str] = None
    parts: List[str] = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None and (wanted is None or name in wanted):
                    seq = "".join(parts)
                    out[name] = seq.upper() if uppercase else seq
                name = line[1:].split()[0]
                parts = []
            else:
                if name is not None:
                    parts.append(line.strip())
        if name is not None and (wanted is None or name in wanted):
            seq = "".join(parts)
            out[name] = seq.upper() if uppercase else seq
    return out


def iter_fasta(path: str) -> Iterator[Tuple[str, str]]:
    """Yield ``(header_first_field, sequence)`` without loading everything."""
    name: Optional[str] = None
    parts: List[str] = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()
                name = line[1:].split()[0]
                parts = []
            else:
                if name is not None:
                    parts.append(line.strip())
        if name is not None:
            yield name, "".join(parts).upper()


def iter_fasta_records(path: str) -> Iterator[Tuple[str, str]]:
    """Yield ``(full_header, sequence)`` without loading everything."""
    header: Optional[str] = None
    parts: List[str] = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts).upper()
                header = line[1:].strip()
                parts = []
            else:
                if header is not None:
                    parts.append(line.strip())
        if header is not None:
            yield header, "".join(parts).upper()


def write_fasta(path: str, records: Iterable[Tuple[str, str]], width: int = 60) -> None:
    with open(path, "wt", encoding="ascii", newline="\n") as handle:
        for name, seq in records:
            handle.write(">" + name + "\n")
            for i in range(0, len(seq), width):
                handle.write(seq[i : i + width] + "\n")


@dataclass(frozen=True)
class FaiEntry:
    name: str
    length: int
    offset: int
    line_bases: int
    line_width: int


class FastaReader:
    """Random-access reader for plain FASTA, or in-memory reader for gzip."""

    def __init__(self, path: str, force_memory: bool = False):
        self.path = os.path.abspath(path)
        self._memory: Optional[Dict[str, str]] = None
        self._index: Dict[str, FaiEntry] = {}
        self._plain_path = self._resolve_plain(path)
        if force_memory or self._plain_path is None:
            sys.stderr.write(
                "[mutatedseq] WARNING: gzipped FASTA will be loaded into memory. "
                "For a human genome, decompress it and run 'samtools faidx' for "
                "low-memory random access.\n"
            )
            self._memory = read_fasta(path)
        else:
            self._ensure_fai(self._plain_path, path)
            self._index = read_fai(self._plain_path + ".fai")

    @staticmethod
    def _resolve_plain(path: str) -> Optional[str]:
        if not path.endswith(".gz"):
            return path
        # If a decompressed sibling exists, prefer it.
        for candidate in (path[:-3], path[:-3] + ".fa", path[:-3] + ".fasta"):
            if os.path.exists(candidate) and os.path.getsize(candidate) > 0:
                return candidate
        return None

    @staticmethod
    def _ensure_fai(plain_path: str, original_path: str) -> None:
        fai = plain_path + ".fai"
        if os.path.exists(fai):
            return
        sys.stderr.write(f"[mutatedseq] building FASTA index: {fai}\n")
        build_fai(plain_path, fai)

    @property
    def is_memory(self) -> bool:
        return self._memory is not None

    def chroms(self) -> List[str]:
        if self._memory is not None:
            return list(self._memory)
        return list(self._index)

    def has_chrom(self, chrom: str) -> bool:
        return self.resolve(chrom) is not None

    def resolve(self, chrom: str) -> Optional[str]:
        """Resolve common chromosome aliases (1/chr1, M/MT/chrM, ...)."""
        names = set(self._memory) if self._memory is not None else set(self._index)
        return resolve_chrom(chrom, names)

    def length(self, chrom: str) -> int:
        name = self.resolve(chrom)
        if name is None:
            raise KeyError(f"chromosome {chrom!r} not present in FASTA")
        if self._memory is not None:
            return len(self._memory[name])
        return self._index[name].length

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """Fetch ``[start, end)``, 0-based.  Out-of-range coordinates are clamped.

        Ensembl's toplevel genome is soft-masked; all sequence returned by this
        class is upper-cased so that VCF REF comparisons are case-insensitive in
        practice.
        """
        name = self.resolve(chrom)
        if name is None:
            raise KeyError(f"chromosome {chrom!r} not present in FASTA")
        if start < 0:
            start = 0
        if self._memory is not None:
            seq = self._memory[name]
            end = min(end, len(seq))
            if start >= end:
                return ""
            return seq[start:end].upper()
        entry = self._index[name]
        if end > entry.length:
            end = entry.length
        if start >= end:
            return ""
        if entry.line_bases <= 0:
            return ""
        pieces: List[bytes] = []
        got = 0
        need = end - start
        with open(self._plain_path, "rb") as handle:  # type: ignore[arg-type]
            while got < need:
                pos = (
                    entry.offset
                    + ((start + got) // entry.line_bases) * entry.line_width
                    + ((start + got) % entry.line_bases)
                )
                handle.seek(pos)
                # Read a little extra to cover newlines in the requested span.
                read_n = need - got + 2 * ((need - got) // entry.line_bases + 2)
                data = handle.read(read_n)
                if not data:
                    break
                data = data.replace(b"\r", b"").replace(b"\n", b"")
                take = data[: need - got]
                if not take:
                    break
                pieces.append(take)
                got += len(take)
        return b"".join(pieces).decode("ascii", errors="replace").upper()

    def sequence(self, chrom: str) -> str:
        """Return a whole chromosome sequence (use with care on large genomes)."""
        name = self.resolve(chrom)
        if name is None:
            raise KeyError(f"chromosome {chrom!r} not present in FASTA")
        if self._memory is not None:
            return self._memory[name]
        return self.fetch(name, 0, self._index[name].length)


def read_fai(path: str) -> Dict[str, FaiEntry]:
    out: Dict[str, FaiEntry] = {}
    with open(path, "rt", encoding="ascii") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            name, length, offset, line_bases, line_width = fields[:5]
            out[name] = FaiEntry(
                name, int(length), int(offset), int(line_bases), int(line_width)
            )
    return out


def build_fai(fasta_path: str, fai_path: Optional[str] = None) -> Dict[str, FaiEntry]:
    """Build a samtools-compatible .fai index for a plain FASTA."""
    if fai_path is None:
        fai_path = fasta_path + ".fai"
    entries: Dict[str, FaiEntry] = {}
    name: Optional[str] = None
    length = 0
    offset = 0
    line_bases = 0
    line_width = 0
    with open(fasta_path, "rb") as handle:
        while True:
            pos = handle.tell()
            line = handle.readline()
            if not line:
                break
            if line.startswith(b">"):
                if name is not None:
                    entries[name] = FaiEntry(name, length, offset, line_bases, line_width)
                name = line[1:].split()[0].decode("ascii", errors="replace")
                length = 0
                offset = handle.tell()
                line_bases = 0
                line_width = 0
            else:
                if not line.strip():
                    # Non-standard blank lines: skip but keep a consistent width.
                    continue
                if line_bases == 0:
                    line_bases = len(line.rstrip(b"\r\n"))
                    line_width = len(line)
                length += len(line.rstrip(b"\r\n"))
        if name is not None:
            entries[name] = FaiEntry(name, length, offset, line_bases, line_width)
    with open(fai_path, "wt", encoding="ascii", newline="\n") as out:
        for entry in entries.values():
            out.write(
                f"{entry.name}\t{entry.length}\t{entry.offset}\t"
                f"{entry.line_bases}\t{entry.line_width}\n"
            )
    return entries


def resolve_chrom(chrom: str, available: Set[str]) -> Optional[str]:
    """Resolve 1/chr1, MT/chrM/M, 23/X etc. against available names."""
    if chrom in available:
        return chrom
    candidates = [chrom]
    if chrom.startswith("chr"):
        candidates.append(chrom[3:])
    else:
        candidates.append("chr" + chrom)
    upper = chrom.upper()
    if upper in ("M", "MT", "CHRM", "CHRMT"):
        candidates.extend(["M", "MT", "chrM", "chrMT"])
    if chrom == "23":
        candidates.extend(["X", "chrX"])
    if chrom == "24":
        candidates.extend(["Y", "chrY"])
    if chrom == "25":
        candidates.extend(["MT", "M", "chrM", "chrMT"])
    for candidate in candidates:
        if candidate in available:
            return candidate
    # Case-insensitive last resort.
    lower_map = {name.lower(): name for name in available}
    return lower_map.get(chrom.lower())


def decompress_gzip(src: str, dst: str, chunk_size: int = 1 << 20) -> None:
    with gzip.open(src, "rb") as inp, open(dst, "wb") as out:
        while True:
            chunk = inp.read(chunk_size)
            if not chunk:
                break
            out.write(chunk)
