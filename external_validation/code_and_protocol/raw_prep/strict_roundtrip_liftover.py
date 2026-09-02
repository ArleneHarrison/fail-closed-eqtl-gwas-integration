#!/usr/bin/env python3
"""Strict point/window liftOver with exact one-to-one round-trip auditing.

Input TSV columns:
  unit_id, chromosome, start_1based, end_1based

The output includes mapped 1-based inclusive boundaries and rejects ambiguous,
cross-chromosome, reversed, or non-exact round-trip mappings. No approximate
same-coordinate window fallback is permitted.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from pyliftover import LiftOver


def canonical_chrom(value: str) -> str:
    value = value.strip()
    if not value.startswith("chr"):
        value = "chr" + value
    return value


def unique_primary(lo: LiftOver, chrom: str, pos1: int):
    hits = lo.convert_coordinate(chrom, pos1 - 1)
    primary = [h for h in hits if h[0] in {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}]
    if len(primary) != 1:
        return None, len(primary)
    h = primary[0]
    return (h[0], int(h[1]) + 1, h[2], h[3]), 1


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--forward-chain", required=True)
    parser.add_argument("--reverse-chain", required=True)
    parser.add_argument("--source-build", required=True)
    parser.add_argument("--target-build", required=True)
    args = parser.parse_args()

    forward = LiftOver(args.forward_chain)
    reverse = LiftOver(args.reverse_chain)
    output_fields = [
        "unit_id", "source_build", "target_build", "source_chromosome",
        "source_start_1based", "source_end_1based", "target_chromosome",
        "target_start_1based", "target_end_1based", "start_forward_hits",
        "end_forward_hits", "start_reverse_hits", "end_reverse_hits",
        "start_roundtrip_exact", "end_roundtrip_exact", "status", "error_code",
    ]
    records = []
    with Path(args.input).open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            unit = row["unit_id"]
            chrom = canonical_chrom(row["chromosome"])
            start = int(row["start_1based"])
            end = int(row["end_1based"])
            rec = {k: "" for k in output_fields}
            rec.update(
                unit_id=unit, source_build=args.source_build, target_build=args.target_build,
                source_chromosome=chrom, source_start_1based=start, source_end_1based=end,
                status="INELIGIBLE", error_code="",
            )
            if start < 1 or end < start:
                rec["error_code"] = "E_SOURCE_WINDOW_INVALID"
                records.append(rec)
                continue
            mapped_start, n_start = unique_primary(forward, chrom, start)
            mapped_end, n_end = unique_primary(forward, chrom, end)
            rec["start_forward_hits"], rec["end_forward_hits"] = n_start, n_end
            if mapped_start is None or mapped_end is None:
                rec["error_code"] = "E_LIFTOVER_NOT_ONE_TO_ONE"
                records.append(rec)
                continue
            tchr_s, tstart, strand_s, _ = mapped_start
            tchr_e, tend, strand_e, _ = mapped_end
            rec.update(target_chromosome=tchr_s, target_start_1based=tstart, target_end_1based=tend)
            if tchr_s != tchr_e:
                rec["error_code"] = "E_LIFTOVER_CROSS_CHROMOSOME"
                records.append(rec)
                continue
            if strand_s != strand_e or strand_s != "+" or tend < tstart:
                rec["error_code"] = "E_LIFTOVER_INTERVAL_ORIENTATION"
                records.append(rec)
                continue
            back_start, nb_start = unique_primary(reverse, tchr_s, tstart)
            back_end, nb_end = unique_primary(reverse, tchr_e, tend)
            rec["start_reverse_hits"], rec["end_reverse_hits"] = nb_start, nb_end
            exact_start = back_start is not None and back_start[0] == chrom and back_start[1] == start
            exact_end = back_end is not None and back_end[0] == chrom and back_end[1] == end
            rec["start_roundtrip_exact"] = str(exact_start).upper()
            rec["end_roundtrip_exact"] = str(exact_end).upper()
            if not exact_start or not exact_end:
                rec["error_code"] = "E_LIFTOVER_ROUNDTRIP_MISMATCH"
                records.append(rec)
                continue
            rec["status"] = "PASS"
            records.append(rec)

    with Path(args.output).open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=output_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


if __name__ == "__main__":
    main()
