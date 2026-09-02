"""Streaming exact-gene eQTL prefiltering for unindexed Catalogue files."""
from __future__ import annotations

from collections.abc import Iterable, Iterator


def filter_rows_for_gene(lines: Iterable[str], gene_id: str) -> Iterator[str]:
    iterator = iter(lines)
    header = next(iterator)
    columns = header.rstrip("\n").split("\t")
    try:
        gene_index = columns.index("gene_id")
    except ValueError as error:
        raise ValueError("input header lacks gene_id") from error
    yield header
    for line in iterator:
        fields = line.rstrip("\n").split("\t")
        if len(fields) > gene_index and fields[gene_index] == gene_id:
            yield line
