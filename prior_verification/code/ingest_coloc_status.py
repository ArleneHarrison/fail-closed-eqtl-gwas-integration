"""Merge a completed regional coloc status into the candidate status registry."""
from __future__ import annotations

import argparse

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry", required=True)
    parser.add_argument("--region-id", required=True)
    parser.add_argument("--status-file", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    registry = pd.read_csv(args.registry)
    status = pd.read_csv(args.status_file).iloc[0]
    mask = registry["region_id"] == args.region_id
    if mask.sum() != 1:
        raise ValueError(f"expected exactly one registry row for {args.region_id}")
    registry.loc[mask, "status"] = status["status"]
    registry.loc[mask, "ineligibility_reason"] = status["reason"]
    registry.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
