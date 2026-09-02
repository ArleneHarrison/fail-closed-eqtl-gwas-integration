"""Small, auditable helpers for ancestry-matched LD reference construction."""
from __future__ import annotations

import re

import numpy as np


_REGION = re.compile(r"^chr?([0-9XYM]+):(\d+)-(\d+)$", re.IGNORECASE)


def parse_region(value: str) -> tuple[str, int, int]:
    """Parse a one-based, inclusive genomic region and reject invalid bounds."""
    match = _REGION.match(value.strip())
    if not match:
        raise ValueError(f"invalid region: {value}")
    chrom, start, end = match.groups()
    start_int, end_int = int(start), int(end)
    if start_int < 1 or end_int < start_int:
        raise ValueError(f"invalid region bounds: {value}")
    return chrom, start_int, end_int


def build_variant_key(chrom: str, position: int, ref: str, alt: str) -> str:
    return f"chr{chrom}_{position}_{ref.upper()}_{alt.upper()}"


def select_requested_variants(
    observed: list[str], requested: set[str]
) -> tuple[list[str], list[str]]:
    """Keep observed variants in requested order and report unmatched requests."""
    observed_set = set(observed)
    selected = [variant for variant in observed if variant in requested]
    missing = sorted(requested.difference(observed_set))
    return selected, missing


def compute_ld(dosages: np.ndarray) -> np.ndarray:
    """Compute pairwise-complete Pearson LD from variants-by-samples dosages."""
    if dosages.ndim != 2 or dosages.shape[0] < 2:
        raise ValueError("dosages must contain at least two variants")
    n_variants = dosages.shape[0]
    ld = np.eye(n_variants, dtype=float)
    for left in range(n_variants):
        for right in range(left + 1, n_variants):
            mask = np.isfinite(dosages[left]) & np.isfinite(dosages[right])
            if mask.sum() < 3:
                correlation = np.nan
            else:
                x, y = dosages[left, mask], dosages[right, mask]
                if np.std(x) == 0 or np.std(y) == 0:
                    correlation = np.nan
                else:
                    correlation = float(np.corrcoef(x, y)[0, 1])
            ld[left, right] = ld[right, left] = correlation
    return ld
