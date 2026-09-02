"""Streaming coordinate-window prefiltering for large TSV GWAS sources."""
from __future__ import annotations

from collections.abc import Iterable, Iterator, Sequence


def resolve_delimiter(value: str) -> str:
    """Resolve a human-readable delimiter name used in regional inputs."""
    return "\t" if value == "tab" else value


def filter_rows_for_regions(
    lines: Iterable[str], regions: Sequence[tuple[str, int, int]]
) -> Iterator[str]:
    iterator = iter(lines)
    header = next(iterator)
    columns = header.rstrip("\n").split("\t")
    chrom_index = columns.index("chromosome")
    position_index = columns.index("base_pair_location")
    yield header
    for line in iterator:
        fields = line.rstrip("\n").split("\t")
        if len(fields) <= max(chrom_index, position_index):
            continue
        chrom, position = fields[chrom_index], int(fields[position_index])
        if any(chrom == region_chrom and start <= position <= end for region_chrom, start, end in regions):
            yield line
