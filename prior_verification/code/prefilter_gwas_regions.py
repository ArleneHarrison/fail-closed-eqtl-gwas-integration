"""Atomically extract one or more coordinate windows from a TSV GWAS file."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
from tempfile import NamedTemporaryFile

from hcsmr.reanalysis.gwas_prefilter import filter_rows_for_regions


def _region(value: str) -> tuple[str, int, int]:
    chrom, bounds = value.split(":", maxsplit=1)
    start, end = bounds.split("-", maxsplit=1)
    return chrom.removeprefix("chr"), int(start), int(end)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True)
    parser.add_argument("--region", action="append", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with Path(args.source).open(encoding="utf-8") as source, NamedTemporaryFile(
        "w", encoding="utf-8", delete=False, dir=out.parent, suffix=".tmp"
    ) as temporary:
        for line in filter_rows_for_regions(source, [_region(item) for item in args.region]):
            temporary.write(line)
        temporary_path = temporary.name
    os.replace(temporary_path, out)


if __name__ == "__main__":
    main()
