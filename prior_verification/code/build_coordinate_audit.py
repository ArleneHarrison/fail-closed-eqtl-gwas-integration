"""Create a source-build audit for every inherited exploratory candidate."""
from __future__ import annotations

import argparse

import pandas as pd

from hcsmr.reanalysis.coordinate_contract import (
    SourceCoordinateRecord,
    build_coordinate_audit,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--sources", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    source_frame = pd.read_csv(args.sources)
    records = [
        SourceCoordinateRecord(row.source, row.genome_build, row.role)
        for row in source_frame.itertuples(index=False)
    ]
    candidates = pd.read_csv(args.candidates)
    build_coordinate_audit(candidates, records).to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
