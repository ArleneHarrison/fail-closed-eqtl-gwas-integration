"""Prespecified simulation grid and uncertainty summaries for hcsMR reanalysis."""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product


@dataclass(frozen=True)
class SimulationScenarioV2:
    n_instruments: int
    pleiotropy: str
    ld: str
    sample_overlap: float


def required_scenarios() -> list[SimulationScenarioV2]:
    """Return the locked stress-test grid required before primary claims.

    The actual runner may add true effect size, eQTL covariance, and prior
    sensitivity subgrids, but it must not remove sparse-IV or adverse
    pleiotropy conditions from this protocol.
    """
    return [
        SimulationScenarioV2(
            n_instruments=n_instruments,
            pleiotropy=pleiotropy,
            ld=ld,
            sample_overlap=sample_overlap,
        )
        for n_instruments, pleiotropy, ld, sample_overlap in product(
            (3, 4, 5, 10, 20),
            ("balanced", "directional", "correlated", "inside_violation"),
            ("independent", "residual_ld"),
            (0.0, 0.25),
        )
    ]


def binomial_interval(
    successes: int, trials: int, z: float = 1.959963984540054
) -> tuple[float, float]:
    """Return a two-sided Wilson confidence interval for a binomial rate."""
    if trials <= 0 or successes < 0 or successes > trials:
        raise ValueError("successes must satisfy 0 <= successes <= trials and trials > 0")
    proportion = successes / trials
    denominator = 1.0 + z * z / trials
    centre = (proportion + z * z / (2.0 * trials)) / denominator
    half_width = (
        z
        * (
            proportion * (1.0 - proportion) / trials
            + z * z / (4.0 * trials * trials)
        )
        ** 0.5
        / denominator
    )
    lower = 0.0 if successes == 0 else max(0.0, centre - half_width)
    upper = 1.0 if successes == trials else min(1.0, centre + half_width)
    return lower, upper
