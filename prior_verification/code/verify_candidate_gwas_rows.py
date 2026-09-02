"""Stream a raw CAD GWAS file and audit every listed candidate row."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from hcsmr.reanalysis.raw_gwas_verification import (
    RAW_COLUMNS,
    candidate_query_keys,
    filter_raw_rows_for_candidates,
    verify_candidate_rows,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--raw-gwas", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--chunksize", type=int, default=500_000)
    args = parser.parse_args()
    candidates = pd.read_csv(args.candidates)
    candidates = candidates.loc[candidates["candidate_id"].astype(str).str.startswith("CAD|")].copy()
    if candidates.empty:
        raise ValueError("no CAD candidates found")
    wanted_markers, wanted_coordinates = candidate_query_keys(candidates)
    matched_chunks = []
    for chunk in pd.read_csv(args.raw_gwas, sep="\t", usecols=list(RAW_COLUMNS), chunksize=args.chunksize):
        matched = filter_raw_rows_for_candidates(chunk, wanted_markers, wanted_coordinates)
        if not matched.empty:
            matched_chunks.append(matched)
    raw_matches = (
        pd.concat(matched_chunks, ignore_index=True)
        if matched_chunks
        else pd.DataFrame(columns=list(RAW_COLUMNS))
    )
    result = verify_candidate_rows(candidates, raw_matches)
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, index=False)
    print(result["match_status"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
