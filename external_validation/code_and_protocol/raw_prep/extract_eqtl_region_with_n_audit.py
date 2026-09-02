#!/usr/bin/env python3
"""Extract a GRCh38 eQTL region and derive per-row effective N from `an/2`.

The extractor refuses the region if any `an` value is missing, non-positive, or
odd. It preserves `an` and writes the explicit source rule into the audit.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import subprocess
from pathlib import Path


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--eqtl", required=True)
    p.add_argument("--region-grch38", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--audit-output", required=True)
    args = p.parse_args()

    eqtl = Path(args.eqtl)
    with gzip.open(eqtl, "rt", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\r\n").split("\t")
    if "an" not in header:
        raise SystemExit("E_EQTL_AN_MISSING")

    result = subprocess.run(
        ["tabix", str(eqtl), args.region_grch38],
        check=True,
        capture_output=True,
        text=True,
    )
    rows = []
    invalid = []
    for line_no, line in enumerate(result.stdout.splitlines(), start=1):
        values = line.split("\t")
        if len(values) != len(header):
            invalid.append({"regional_line": line_no, "error": "E_EQTL_FIELD_COUNT"})
            continue
        row = dict(zip(header, values))
        try:
            an = int(row["an"])
        except Exception:
            invalid.append({"regional_line": line_no, "error": "E_EQTL_AN_NOT_INTEGER", "an": row.get("an")})
            continue
        if an <= 0 or an % 2:
            invalid.append({"regional_line": line_no, "error": "E_EQTL_AN_NOT_POSITIVE_EVEN", "an": an})
            continue
        row["effective_n"] = str(an // 2)
        row["effective_n_source"] = "eQTL_Catalogue_an_divided_by_2"
        rows.append(row)

    output_header = header + ["effective_n", "effective_n_source"]
    with Path(args.output).open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=output_header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    n_values = sorted({int(r["effective_n"]) for r in rows})
    audit = {
        "eqtl_path": str(eqtl),
        "region_grch38": args.region_grch38,
        "n_source_field": "an",
        "n_derivation_rule": "effective_n = an / 2",
        "validation_rule": "an must be a positive even integer for every retained row",
        "rows_returned_by_tabix": len(result.stdout.splitlines()),
        "rows_retained": len(rows),
        "invalid_rows": invalid,
        "effective_n_unique": n_values,
        "status": "PASS" if not invalid and rows else ("E_NO_ROWS" if not rows and not invalid else "FAIL"),
    }
    with Path(args.audit_output).open("w", encoding="utf-8") as fh:
        json.dump(audit, fh, indent=2)
        fh.write("\n")
    if audit["status"] != "PASS":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
