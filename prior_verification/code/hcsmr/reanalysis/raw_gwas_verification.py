"""Verify candidate-level GWAS values against an immutable raw summary file."""
from __future__ import annotations

import re

import pandas as pd


RAW_COLUMNS = (
    "chromosome",
    "base_pair_location",
    "effect_allele",
    "other_allele",
    "beta",
    "standard_error",
    "p_value",
    "markername",
)


def _expected_marker(variant: str) -> str:
    match = re.fullmatch(r"chr([^_]+)_([0-9]+)_([^_]+)_([^_]+)", str(variant))
    if match is None:
        raise ValueError(f"invalid candidate variant identifier: {variant!r}")
    chrom, position, ref, alt = match.groups()
    return f"{chrom}:{position}_{ref}_{alt}"


def _coordinate_key(variant: str) -> tuple[str, int]:
    marker = _expected_marker(variant)
    coordinate = marker.split("_", 1)[0]
    chrom, position = coordinate.split(":", 1)
    return chrom, int(position)


def _same_number(left: object, right: object, tolerance: float = 1e-12) -> bool:
    try:
        return bool(abs(float(left) - float(right)) <= tolerance)
    except (TypeError, ValueError):
        return False


def candidate_query_keys(candidates: pd.DataFrame) -> tuple[set[str], set[tuple[str, int]]]:
    """Derive exact-marker and coordinate keys without assuming allele agreement."""
    if "variant" not in candidates.columns:
        raise ValueError("candidate table is missing: ['variant']")
    markers = {_expected_marker(value) for value in candidates["variant"]}
    coordinates = {_coordinate_key(value) for value in candidates["variant"]}
    return markers, coordinates


def filter_raw_rows_for_candidates(
    raw: pd.DataFrame, markers: set[str], coordinates: set[tuple[str, int]]
) -> pd.DataFrame:
    """Keep exact markers and alternative alleles at candidate coordinates."""
    required = {"chromosome", "base_pair_location", "markername"}
    missing = required.difference(raw.columns)
    if missing:
        raise ValueError(f"raw GWAS table is missing: {sorted(missing)}")
    coordinate_keys = list(
        zip(
            raw["chromosome"].astype(str),
            pd.to_numeric(raw["base_pair_location"], errors="raise").astype(int),
        )
    )
    keep = raw["markername"].astype(str).isin(markers) | pd.Series(
        [key in coordinates for key in coordinate_keys], index=raw.index
    )
    return raw.loc[keep].copy()


def _allele_relation(candidate_ref: str, candidate_alt: str, raw_effect: object, raw_other: object) -> str:
    effect, other = str(raw_effect).upper(), str(raw_other).upper()
    if effect == candidate_alt and other == candidate_ref:
        return "effect_matches_candidate_alt"
    if effect == candidate_ref and other == candidate_alt:
        return "effect_matches_candidate_ref"
    return "allele_set_or_representation_mismatch"


def _resolution_for_relation(relation: str) -> str:
    if relation == "effect_matches_candidate_alt":
        return "none"
    if relation == "effect_matches_candidate_ref":
        return "validate_build_strand_and_indel_normalization_before_effect_flip"
    return "validate_build_rsid_strand_and_indel_normalization"


def verify_candidate_rows(candidates: pd.DataFrame, raw: pd.DataFrame) -> pd.DataFrame:
    """Return one row per candidate, retaining missing and discordant matches."""
    required_candidate_columns = {
        "candidate_id", "variant", "gwas_beta", "gwas_se", "gwas_p"
    }
    required_raw_columns = set(RAW_COLUMNS)
    missing_candidate = required_candidate_columns.difference(candidates.columns)
    missing_raw = required_raw_columns.difference(raw.columns)
    if missing_candidate:
        raise ValueError(f"candidate table is missing: {sorted(missing_candidate)}")
    if missing_raw:
        raise ValueError(f"raw GWAS table is missing: {sorted(missing_raw)}")

    raw_by_marker = {
        str(marker): group.copy()
        for marker, group in raw.groupby(raw["markername"].astype(str), sort=False)
    }
    raw_by_coordinate = {
        (str(chrom), int(position)): group.copy()
        for (chrom, position), group in raw.groupby(
            [raw["chromosome"].astype(str), pd.to_numeric(raw["base_pair_location"], errors="raise")],
            sort=False,
        )
    }
    rows: list[dict[str, object]] = []
    for candidate in candidates.itertuples(index=False):
        candidate_dict = candidate._asdict()
        expected_marker = _expected_marker(candidate_dict["variant"])
        _, _, candidate_ref, candidate_alt = re.fullmatch(
            r"chr([^_]+)_([0-9]+)_([^_]+)_([^_]+)", str(candidate_dict["variant"])
        ).groups()
        candidate_ref, candidate_alt = candidate_ref.upper(), candidate_alt.upper()
        expected_chrom, expected_position = _coordinate_key(candidate_dict["variant"])
        matching_rows = raw_by_marker.get(expected_marker)
        if matching_rows is None:
            coordinate_rows = raw_by_coordinate.get((expected_chrom, expected_position))
            if coordinate_rows is None:
                rows.append(
                    {
                        "candidate_id": candidate_dict["candidate_id"],
                        "variant": candidate_dict["variant"],
                        "expected_markername": expected_marker,
                        "raw_markername": "",
                        "raw_row_count_at_coordinate": 0,
                        "match_status": "raw_row_not_found",
                        "effect_allele_match": False,
                        "other_allele_match": False,
                        "beta_match": False,
                        "se_match": False,
                        "p_value_match": False,
                        "raw_effect_alleles": "",
                        "raw_other_alleles": "",
                        "candidate_ref": candidate_ref,
                        "candidate_alt": candidate_alt,
                        "allele_relation": "raw_row_not_found",
                        "required_resolution": "validate_build_and_rsid_mapping",
                        "raw_beta": float("nan"),
                        "raw_se": float("nan"),
                        "raw_p_value": float("nan"),
                    }
                )
                continue
            coordinate_source = coordinate_rows.iloc[0]
            coordinate_relation = (
                _allele_relation(
                    candidate_ref,
                    candidate_alt,
                    coordinate_source["effect_allele"],
                    coordinate_source["other_allele"],
                )
                if len(coordinate_rows) == 1
                else "multiple_raw_rows_at_coordinate"
            )
            coordinate_resolution = (
                _resolution_for_relation(coordinate_relation)
                if len(coordinate_rows) == 1
                else "validate_build_rsid_strand_and_indel_normalization"
            )
            rows.append(
                {
                    "candidate_id": candidate_dict["candidate_id"],
                    "variant": candidate_dict["variant"],
                    "expected_markername": expected_marker,
                    "raw_markername": ";".join(coordinate_rows["markername"].astype(str)),
                    "raw_row_count_at_coordinate": int(len(coordinate_rows)),
                    "match_status": "coordinate_found_marker_allele_mismatch",
                    "effect_allele_match": False,
                    "other_allele_match": False,
                    "beta_match": False,
                    "se_match": False,
                    "p_value_match": False,
                    "raw_effect_alleles": ";".join(coordinate_rows["effect_allele"].astype(str)),
                    "raw_other_alleles": ";".join(coordinate_rows["other_allele"].astype(str)),
                    "candidate_ref": candidate_ref,
                    "candidate_alt": candidate_alt,
                    "allele_relation": coordinate_relation,
                    "required_resolution": coordinate_resolution,
                    "raw_beta": float(coordinate_source["beta"]),
                    "raw_se": float(coordinate_source["standard_error"]),
                    "raw_p_value": float(coordinate_source["p_value"]),
                }
            )
            continue
        source = matching_rows.iloc[0]
        relation = _allele_relation(candidate_ref, candidate_alt, source["effect_allele"], source["other_allele"])
        effect_match = relation == "effect_matches_candidate_alt"
        other_match = effect_match
        beta_match = _same_number(candidate_dict["gwas_beta"], source["beta"])
        se_match = _same_number(candidate_dict["gwas_se"], source["standard_error"])
        p_value_match = _same_number(candidate_dict["gwas_p"], source["p_value"])
        status = (
            "verified"
            if all((effect_match, other_match, beta_match, se_match, p_value_match))
            else "raw_row_value_or_allele_mismatch"
        )
        rows.append(
            {
                "candidate_id": candidate_dict["candidate_id"],
                "variant": candidate_dict["variant"],
                "expected_markername": expected_marker,
                "raw_markername": str(source["markername"]),
                "raw_row_count_at_coordinate": int(len(raw_by_coordinate[(expected_chrom, expected_position)])),
                "match_status": status,
                "effect_allele_match": effect_match,
                "other_allele_match": other_match,
                "beta_match": beta_match,
                "se_match": se_match,
                "p_value_match": p_value_match,
                "raw_effect_alleles": str(source["effect_allele"]),
                "raw_other_alleles": str(source["other_allele"]),
                "candidate_ref": candidate_ref,
                "candidate_alt": candidate_alt,
                "allele_relation": relation,
                "required_resolution": _resolution_for_relation(relation),
                "raw_beta": float(source["beta"]),
                "raw_se": float(source["standard_error"]),
                "raw_p_value": float(source["p_value"]),
            }
        )
    return pd.DataFrame(rows)


def summarize_candidate_qc(qc_rows: pd.DataFrame) -> pd.DataFrame:
    """Close a candidate-level input gate only when every row is verified."""
    required = {"candidate_id", "match_status"}
    missing = required.difference(qc_rows.columns)
    if missing:
        raise ValueError(f"QC table is missing: {sorted(missing)}")
    rows: list[dict[str, object]] = []
    for candidate_id, group in qc_rows.groupby("candidate_id", sort=False):
        statuses = group["match_status"].astype(str).tolist()
        failed = sorted({status for status in statuses if status != "verified"})
        rows.append(
            {
                "candidate_id": candidate_id,
                "raw_gwas_qc_rows": int(len(group)),
                "raw_gwas_qc_verified_rows": int(sum(status == "verified" for status in statuses)),
                "raw_gwas_qc_pass": not failed,
                "raw_gwas_qc_status": "pass" if not failed else "blocked_" + ";".join(failed),
            }
        )
    return pd.DataFrame(rows)
