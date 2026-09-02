"""Summarize archived multi-resource workflow statuses without biological interpretation."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--status", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    status_path = Path(args.status)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(status_path)
    parsed = frame["region_id"].str.split("|", n=2, expand=True)
    frame[["gwas", "gene", "context"]] = parsed
    frame["resource_family"] = frame["context"].str.split("_", n=1).str[0]
    summary = (
        frame.groupby(["gwas", "resource_family", "status"], dropna=False)
        .size()
        .reset_index(name="cases")
        .sort_values(["gwas", "resource_family", "status"])
    )
    summary.to_csv(out_dir / "cross_resource_status_summary.tsv", sep="\t", index=False)
    frame[["region_id", "gwas", "gene", "context", "resource_family", "status", "ineligibility_reason"]].to_csv(
        out_dir / "cross_resource_status_cases.tsv", sep="\t", index=False
    )
    digest = hashlib.sha256(status_path.read_bytes()).hexdigest()
    manifest = {
        "source_sha256": digest,
        "cases": int(len(frame)),
        "gwas_outcomes": sorted(frame["gwas"].unique().tolist()),
        "resource_families": sorted(frame["resource_family"].unique().tolist()),
        "contexts": sorted(frame["context"].unique().tolist()),
        "status_counts": {str(key): int(value) for key, value in frame["status"].value_counts().items()},
        "claim_boundary": (
            "Retrospective portability/status audit only. Posterior values are deliberately not summarized, "
            "and no association, colocalization, tissue, or biological claim is made."
        ),
    }
    (out_dir / "cross_resource_status_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, sort_keys=True))


if __name__ == "__main__":
    main()
