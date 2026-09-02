"""Small, explicit allele transformations used after coordinate liftover."""
from __future__ import annotations


_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def complement_allele(allele: str) -> str:
    normalized = allele.upper()
    if any(base not in "ACGT" for base in normalized):
        raise ValueError(f"non-DNA allele cannot be complemented: {allele}")
    return normalized.translate(_COMPLEMENT)


def transformed_alleles(ref: str, alt: str, strand: str) -> tuple[str, str]:
    if strand == "+":
        return ref.upper(), alt.upper()
    if strand == "-":
        return complement_allele(ref), complement_allele(alt)
    raise ValueError(f"unexpected liftover strand: {strand}")
