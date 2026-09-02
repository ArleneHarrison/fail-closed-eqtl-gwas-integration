"""Lift complete eQTL source rows from GRCh38 to GRCh37 with no silent drops."""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

from pyliftover import LiftOver


_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def transform(allele: str, strand: str) -> str:
    return allele.translate(_COMPLEMENT) if strand == "-" else allele


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--eqtl", required=True)
    parser.add_argument("--chain", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--summary-out", required=True)
    args = parser.parse_args()
    converter = LiftOver(args.chain)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    summary_out = Path(args.summary_out)
    summary_out.parent.mkdir(parents=True, exist_ok=True)
    status_counts: Counter[str] = Counter()
    target_positions: list[int] = []

    with Path(args.eqtl).open(encoding="utf-8", newline="") as source, out.open(
        "w", encoding="utf-8", newline=""
    ) as destination:
        reader = csv.DictReader(source, delimiter="\t")
        required = {"variant", "chromosome", "position", "ref", "alt"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"eQTL input lacks required columns: {sorted(missing)}")
        fieldnames = list(reader.fieldnames or []) + [
            "liftover_status", "liftover_mapping_count", "target_chromosome",
            "target_position", "target_ref", "target_alt", "target_variant", "liftover_strand",
            "liftover_query_chromosome",
        ]
        writer = csv.DictWriter(destination, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in reader:
            source_chrom = str(row["chromosome"]).removeprefix("chr")
            try:
                source_pos = int(row["position"])
            except ValueError:
                row.update(liftover_status="invalid_source_position", liftover_mapping_count="0")
                writer.writerow(row)
                status_counts["invalid_source_position"] += 1
                continue
            # UCSC chains typically use `chr11`, whereas some liftover APIs
            # accept bare `11`.  Record the convention that actually maps.
            mappings = converter.convert_coordinate(f"chr{source_chrom}", source_pos - 1) or []
            query_chrom = f"chr{source_chrom}"
            if not mappings:
                mappings = converter.convert_coordinate(source_chrom, source_pos - 1) or []
                query_chrom = source_chrom
            row["liftover_query_chromosome"] = query_chrom
            row["liftover_mapping_count"] = str(len(mappings))
            if len(mappings) != 1:
                row["liftover_status"] = "unmapped" if not mappings else "ambiguous_mapping"
                writer.writerow(row)
                status_counts[row["liftover_status"]] += 1
                continue
            target_chrom, target_zero, strand, _score = mappings[0]
            target_pos = target_zero + 1
            target_ref = transform(str(row["ref"]), strand)
            target_alt = transform(str(row["alt"]), strand)
            row.update(
                liftover_status="mapped",
                target_chromosome=str(target_chrom).removeprefix("chr"),
                target_position=str(target_pos),
                target_ref=target_ref,
                target_alt=target_alt,
                target_variant=f"chr{str(target_chrom).removeprefix('chr')}_{target_pos}_{target_ref}_{target_alt}",
                liftover_strand=strand,
            )
            writer.writerow(row)
            status_counts["mapped"] += 1
            target_positions.append(target_pos)

    with summary_out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "value"])
        writer.writeheader()
        for status, count in sorted(status_counts.items()):
            writer.writerow({"metric": f"rows_{status}", "value": count})
        writer.writerow({"metric": "mapped_target_position_min", "value": min(target_positions) if target_positions else ""})
        writer.writerow({"metric": "mapped_target_position_max", "value": max(target_positions) if target_positions else ""})
        writer.writerow({"metric": "allele_handling", "value": "negative-strand mappings are complemented; all unmapped and ambiguous rows retained"})


if __name__ == "__main__":
    main()
