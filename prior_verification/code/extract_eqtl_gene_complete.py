"""Stream a complete eQTL source and retain every row for one molecular trait.

The source archive remains compressed and immutable.  This utility does not
filter by a hand-selected genomic interval or silently deduplicate rows: the
source-provided cis associations for the requested trait are written verbatim,
with a compact audit record describing the scan.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True)
    parser.add_argument("--gene-id", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--audit-out", required=True)
    args = parser.parse_args()

    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    audit_path = Path(args.audit_out)
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    input_rows = 0
    retained_rows = 0
    malformed_rows = 0
    positions: list[int] = []
    trait_counts: Counter[str] = Counter()

    with gzip.open(args.source, "rt", encoding="utf-8", newline="") as source:
        reader = csv.reader(source, delimiter="\t")
        header = next(reader)
        required = {"gene_id", "molecular_trait_id", "chromosome", "position", "variant"}
        missing = required.difference(header)
        if missing:
            raise ValueError(f"source header missing required columns: {sorted(missing)}")
        width = len(header)
        gene_index = header.index("gene_id")
        trait_index = header.index("molecular_trait_id")
        position_index = header.index("position")
        with output.open("w", encoding="utf-8", newline="") as destination:
            writer = csv.writer(destination, delimiter="\t")
            writer.writerow(header)
            for row in reader:
                input_rows += 1
                if len(row) != width:
                    malformed_rows += 1
                    continue
                if row[gene_index] != args.gene_id:
                    continue
                retained_rows += 1
                trait_counts[row[trait_index]] += 1
                try:
                    positions.append(int(row[position_index]))
                except ValueError:
                    malformed_rows += 1
                writer.writerow(row)

    audit = {
        "source": args.source,
        "requested_gene_id": args.gene_id,
        "header": header,
        "source_rows_scanned_excluding_header": input_rows,
        "retained_rows": retained_rows,
        "malformed_width_or_position_rows": malformed_rows,
        "molecular_trait_row_counts": dict(trait_counts),
        "source_position_min": min(positions) if positions else None,
        "source_position_max": max(positions) if positions else None,
        "deduplication": "none; all source rows for requested gene retained",
        "status": "complete_stream_scan",
    }
    audit_path.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(audit, sort_keys=True))


if __name__ == "__main__":
    main()
