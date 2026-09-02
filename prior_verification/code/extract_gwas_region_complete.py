"""Stream an immutable TSV GWAS and retain every row in a locked region."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True)
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--start", required=True, type=int)
    parser.add_argument("--end", required=True, type=int)
    parser.add_argument("--out", required=True)
    parser.add_argument("--audit-out", required=True)
    args = parser.parse_args()
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    audit_path = Path(args.audit_out)
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    scanned = retained = malformed = 0
    positions: list[int] = []
    with Path(args.source).open(encoding="utf-8", newline="") as source:
        reader = csv.reader(source, delimiter="\t")
        header = next(reader)
        required = {"chromosome", "base_pair_location", "markername", "effect_allele", "other_allele", "beta", "standard_error", "n"}
        missing = required.difference(header)
        if missing:
            raise ValueError(f"GWAS header missing required columns: {sorted(missing)}")
        chrom_index = header.index("chromosome")
        position_index = header.index("base_pair_location")
        width = len(header)
        with output.open("w", encoding="utf-8", newline="") as destination:
            writer = csv.writer(destination, delimiter="\t")
            writer.writerow(header)
            for row in reader:
                scanned += 1
                if len(row) != width:
                    malformed += 1
                    continue
                try:
                    position = int(row[position_index])
                except ValueError:
                    malformed += 1
                    continue
                if str(row[chrom_index]).removeprefix("chr") != str(args.chromosome).removeprefix("chr"):
                    continue
                if not args.start <= position <= args.end:
                    continue
                retained += 1
                positions.append(position)
                writer.writerow(row)
    audit = {
        "source": args.source,
        "locked_region": f"{args.chromosome}:{args.start}-{args.end}",
        "header": header,
        "source_rows_scanned_excluding_header": scanned,
        "retained_rows": retained,
        "malformed_width_or_position_rows": malformed,
        "retained_position_min": min(positions) if positions else None,
        "retained_position_max": max(positions) if positions else None,
        "deduplication": "none; all raw GWAS rows in locked region retained",
        "status": "complete_stream_scan",
    }
    audit_path.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(audit, sort_keys=True))


if __name__ == "__main__":
    main()
