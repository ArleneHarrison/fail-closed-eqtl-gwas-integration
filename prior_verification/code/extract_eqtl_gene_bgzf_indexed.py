"""Extract a fixed eQTL region with the project-local Tabix/BGZF fallback."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path

from hcsmr.reanalysis.tabix_bgzf import fetch_region


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True)
    parser.add_argument("--index", required=True)
    parser.add_argument("--region", required=True, help="one-based inclusive, e.g. 11:2000000-2100000")
    parser.add_argument("--gene-id", required=True, help="gene ID, or __AUTO_LEXICOGRAPHIC_FIRST__ for the locked coordinate-control rule")
    parser.add_argument("--out", required=True)
    parser.add_argument("--audit-out", required=True)
    args = parser.parse_args()
    chrom, span = args.region.split(":", 1)
    start, end = map(int, span.split("-", 1))
    with gzip.open(args.source, "rt", encoding="utf-8", newline="") as source:
        header = next(csv.reader(source, delimiter="\t"))
    chunk_rows = [line.split("\t") for line in fetch_region(args.source, args.index, chrom, start, end)]
    required = {"gene_id", "molecular_trait_id", "chromosome", "position", "variant"}
    missing = required.difference(header)
    if missing:
        raise ValueError(f"header missing required columns: {sorted(missing)}")
    width, gene_index, trait_index, position_index, chromosome_index = (
        len(header), header.index("gene_id"), header.index("molecular_trait_id"), header.index("position"), header.index("chromosome")
    )
    malformed = sum(len(row) != width for row in chunk_rows)
    rows = []
    for row in chunk_rows:
        if len(row) != width:
            continue
        try:
            position = int(row[position_index])
        except ValueError:
            malformed += 1
            continue
        if str(row[chromosome_index]).removeprefix("chr") == chrom.removeprefix("chr") and start <= position <= end:
            rows.append(row)
    gene_counts = Counter(row[gene_index] for row in rows if row[gene_index])
    if args.gene_id == "__AUTO_LEXICOGRAPHIC_FIRST__":
        eligible_genes = sorted(gene for gene, count in gene_counts.items() if count >= 2)
        if not eligible_genes:
            raise ValueError("no gene has at least two complete-width records in the locked region")
        selected_gene = eligible_genes[0]
        selection_rule = "lexicographically_first_nonempty_gene_id_with_at_least_two_complete_width_rows"
    else:
        selected_gene = args.gene_id
        selection_rule = "explicit_gene_id"
    retained = [row for row in rows if row[gene_index] == selected_gene]
    output, audit = Path(args.out), Path(args.audit_out)
    output.parent.mkdir(parents=True, exist_ok=True)
    audit.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerows(retained)
    payload = {
        "access_method": "project_local_tbi_bgzf_reader", "source": args.source, "index": args.index,
        "source_region": args.region, "requested_gene_id": args.gene_id, "selected_gene_id": selected_gene, "selection_rule": selection_rule,
        "chunk_rows_returned": len(chunk_rows), "region_rows_returned": len(rows), "complete_width_gene_row_counts": dict(sorted(gene_counts.items())),
        "retained_rows": len(retained), "malformed_width_rows": malformed,
        "molecular_trait_row_counts": dict(Counter(row[trait_index] for row in retained)),
        "source_position_min": min((int(row[position_index]) for row in retained), default=None),
        "source_position_max": max((int(row[position_index]) for row in retained), default=None),
        "deduplication": "none; every requested-gene row in the complete fixed region is retained",
        "status": "complete_indexed_region_extract",
    }
    audit.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, sort_keys=True))


if __name__ == "__main__":
    main()
