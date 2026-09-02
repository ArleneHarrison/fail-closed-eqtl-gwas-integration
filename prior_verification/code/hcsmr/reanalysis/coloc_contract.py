"""Input validation for formal multi-signal colocalization jobs."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping


@dataclass(frozen=True)
class ColocResult:
    status: str
    ineligibility_reason: str = ""


def validate_coloc_inputs(region: Mapping[str, object]) -> ColocResult:
    """Validate that a region has an LD reference before invoking coloc-SuSiE.

    This only validates inputs. Posterior probabilities are written later by a
    compatible, version-pinned `coloc.susie` execution and are never inferred
    from PIP or GWAS p-value proxies.
    """
    if not region.get("ld_path"):
        return ColocResult("ineligible", "missing_ld_reference")
    if not region.get("ld_ancestry"):
        return ColocResult("ineligible", "missing_ld_ancestry")
    if not Path(str(region["ld_path"])).is_file():
        return ColocResult("ineligible", "ld_file_not_found")
    if not region.get("eqtl_path") or not region.get("gwas_path"):
        return ColocResult("ineligible", "missing_summary_statistics")
    return ColocResult("ready")
