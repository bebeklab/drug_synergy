"""Equal-weight pairwise QUBO scoring and exact ranking."""

from __future__ import annotations

import pandas as pd


PRIMARY_COSTS = (
    "gi50_percentile_cost",
    "pathway_percentile_cost",
    "mechanism_percentile_cost",
)


def score_equal_weight_qubo(
    frame: pd.DataFrame, cost_columns: tuple[str, ...] = PRIMARY_COSTS
) -> pd.DataFrame:
    """Assign equal-weight energy and exact rank to every feasible drug pair."""
    result = frame.copy()
    result["qubo_energy"] = result[list(cost_columns)].mean(axis=1)
    result["qubo_rank"] = result["qubo_energy"].rank(
        method="min", ascending=True
    ).astype(int)
    return result
