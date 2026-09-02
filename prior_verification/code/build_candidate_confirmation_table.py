"""Build the locked-evidence gate table used by the revised manuscript."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from hcsmr.reanalysis.candidate_confirmation import build_confirmation_table
from hcsmr.reanalysis.context_registry import load_context_registry


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--diagnostics", required=True)
    parser.add_argument("--contexts", required=True)
    parser.add_argument("--coloc", required=True)
    parser.add_argument("--coordinates", default=None)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    result = build_confirmation_table(
        pd.read_csv(args.diagnostics),
        load_context_registry(args.contexts),
        pd.read_csv(args.coloc),
        pd.read_csv(args.coordinates) if args.coordinates else None,
    )
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(out, index=False)


if __name__ == "__main__":
    main()
