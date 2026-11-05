"""Construction of the three pair-level QUBO inputs."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .io import cancer_cells, canonical_key, order_pairs
from .paths import INPUTS


def gi50_profiles(cancer: str, standardize: bool = True) -> dict[str, np.ndarray]:
    """Return median-imputed cancer-panel GI50 vectors for every available drug."""
    matrix = pd.read_csv(INPUTS / "nci60_gi50_matrix.csv")
    cells = cancer_cells(list(matrix.columns[2:]), cancer)
    values = matrix[cells].apply(pd.to_numeric, errors="coerce")
    values = values.fillna(values.median(axis=0))
    if standardize:
        scale = values.std(axis=0, ddof=0).replace(0, 1)
        values = (values - values.mean(axis=0)) / scale
    return {
        canonical_key(drug): row.to_numpy(float)
        for drug, (_, row) in zip(matrix["almanac_drug"], values.iterrows())
    }


def add_gi50_activity_cost(frame: pd.DataFrame, cancer: str) -> pd.DataFrame:
    """Add negative mean joint single-drug activity; lower values are preferred."""
    profiles = gi50_profiles(cancer, standardize=True)
    result = frame.copy()
    result["gi50_activity_cost"] = [
        -0.5
        * (
            float(profiles[canonical_key(row.drug1)].mean())
            + float(profiles[canonical_key(row.drug2)].mean())
        )
        for row in result.itertuples(index=False)
    ]
    return result


def load_primary_pair_inputs() -> pd.DataFrame:
    """Merge outcomes, pathway-module costs, and mechanism-redundancy costs."""
    from .io import load_outcomes

    outcomes = load_outcomes()
    pathways = order_pairs(pd.read_csv(INPUTS / "pathway_module_pair_costs.csv"))
    mechanisms = order_pairs(pd.read_csv(INPUTS / "mechanism_redundancy_pair_costs.csv"))
    merged = outcomes.merge(
        pathways[["cancer", "drug1", "drug2", "structured_pathway_cost"]],
        on=["cancer", "drug1", "drug2"],
        how="left",
        validate="one_to_one",
    ).merge(
        mechanisms[["drug1", "drug2", "mechanism_redundancy"]],
        on=["drug1", "drug2"],
        how="left",
        validate="many_to_one",
    )
    required = ["structured_pathway_cost", "mechanism_redundancy"]
    if merged[required].isna().any().any():
        raise ValueError(f"Missing required inputs:\n{merged[required].isna().sum()}")
    return merged


def percentile_rank_costs(frame: pd.DataFrame) -> pd.DataFrame:
    """Convert the three raw costs to comparable within-cancer percentiles."""
    result = frame.copy()
    mappings = {
        "gi50_activity_cost": "gi50_percentile_cost",
        "structured_pathway_cost": "pathway_percentile_cost",
        "mechanism_redundancy": "mechanism_percentile_cost",
    }
    for raw, percentile in mappings.items():
        result[percentile] = result[raw].rank(method="average", pct=True)
    return result
