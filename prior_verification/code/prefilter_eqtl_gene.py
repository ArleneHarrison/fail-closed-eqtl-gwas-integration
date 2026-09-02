"""Atomically create a gene-specific local eQTL table from a gzip source."""
from __future__ import annotations

import argparse
import gzip
import os
from pathlib import Path
from tempfile import NamedTemporaryFile

from hcsmr.reanalysis.eqtl_prefilter import filter_rows_for_gene


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True)
    parser.add_argument("--gene-id", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.source, "rt", encoding="utf-8") as source, NamedTemporaryFile(
        "w", encoding="utf-8", delete=False, dir=out.parent, suffix=".tmp"
    ) as temporary:
        for line in filter_rows_for_gene(source, args.gene_id):
            temporary.write(line)
        temporary_path = temporary.name
    os.replace(temporary_path, out)


if __name__ == "__main__":
    main()
