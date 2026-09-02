"""Extract every requested-gene row from an indexed, pre-specified eQTL region."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import subprocess
from collections import Counter
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tabix", required=True)
    parser.add_argument("--source", required=True)
    parser.add_argument("--region", required=True)
    parser.add_argument("--gene-id", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--audit-out", required=True)
    args = parser.parse_args()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    audit_out = Path(args.audit_out)
    audit_out.parent.mkdir(parents=True, exist_ok=True)
    # This archive's unprefixed header is not returned by tabix.  Read exactly
    # that header from the immutable gzip source, then stream indexed data rows.
    with gzip.open(args.source, "rt", encoding="utf-8", newline="") as source_handle:
        header = next(csv.reader(source_handle, delimiter="\t"))
    command = [args.tabix, args.source, args.region]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    assert process.stdout is not None
    reader = csv.reader(process.stdout, delimiter="\t")
    required = {"gene_id", "molecular_trait_id", "chromosome", "position", "variant"}
    missing = required.difference(header)
    if missing:
        process.kill()
        raise ValueError(f"tabix header missing required columns: {sorted(missing)}")
    gene_index = header.index("gene_id")
    trait_index = header.index("molecular_trait_id")
    width = len(header)
    region_rows = retained_rows = malformed_rows = 0
    trait_counts: Counter[str] = Counter()
    positions: list[int] = []
    with out.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.writer(destination, delimiter="\t")
        writer.writerow(header)
        for row in reader:
            region_rows += 1
            if len(row) != width:
                malformed_rows += 1
                continue
            if row[gene_index] != args.gene_id:
                continue
            retained_rows += 1
            trait_counts[row[trait_index]] += 1
            positions.append(int(row[header.index("position")]))
            writer.writerow(row)
    stderr = process.stderr.read() if process.stderr is not None else ""
    exit_code = process.wait()
    if exit_code:
        raise RuntimeError(f"tabix exit={exit_code}: {stderr.strip()}")
    audit = {
        "command": command,
        "requested_gene_id": args.gene_id,
        "source_region": args.region,
        "header": header,
        "region_rows_returned": region_rows,
        "retained_rows": retained_rows,
        "malformed_width_rows": malformed_rows,
        "molecular_trait_row_counts": dict(trait_counts),
        "source_position_min": min(positions) if positions else None,
        "source_position_max": max(positions) if positions else None,
        "deduplication": "none; all requested-gene rows returned by the complete pre-specified tabix region retained",
        "status": "complete_indexed_region_extract",
    }
    audit_out.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(audit, sort_keys=True))


if __name__ == "__main__":
    main()
