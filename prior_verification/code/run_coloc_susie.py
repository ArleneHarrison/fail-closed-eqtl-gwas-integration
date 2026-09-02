"""Validate and dispatch formal coloc-SuSiE regional analyses.

The script intentionally writes `ineligible` records for regions that lack
ancestry-matched LD instead of producing a PIP-based substitute. A production
coloc execution is allowed only through a version-pinned R environment that
provides `coloc.susie` and compatible SuSiE credible sets.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from hcsmr.reanalysis.coloc_contract import validate_coloc_inputs


OUTPUT_COLUMNS = [
    "region_id",
    "status",
    "ineligibility_reason",
    "pp_h0",
    "pp_h1",
    "pp_h2",
    "pp_h3",
    "pp_h4",
    "credible_set_ids",
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, help="JSON list of regional input records")
    parser.add_argument("--out", default=None, help="CSV output path")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()

    regions = json.loads(Path(args.manifest).read_text(encoding="utf-8"))
    if not isinstance(regions, list):
        raise ValueError("colocalization manifest must contain a JSON list")
    rows = []
    for region in regions:
        result = validate_coloc_inputs(region)
        rows.append(
            {
                "region_id": region.get("region_id", "unknown"),
                "status": result.status,
                "ineligibility_reason": result.ineligibility_reason,
                "pp_h0": None,
                "pp_h1": None,
                "pp_h2": None,
                "pp_h3": None,
                "pp_h4": None,
                "credible_set_ids": "",
            }
        )
    frame = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(args.out, index=False)
    print(frame.to_string(index=False))
    if not args.validate_only and (frame["status"] == "ready").any():
        raise RuntimeError(
            "Validated regions require the version-pinned R coloc.susie runner; "
            "this Python command does not fabricate colocalization posteriors."
        )


if __name__ == "__main__":
    main()
