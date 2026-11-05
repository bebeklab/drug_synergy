"""Input loading and deterministic drug-pair normalization."""

from __future__ import annotations

import json

import pandas as pd

from .paths import CONFIG, INPUTS


CANCER_PREFIX = {
    "Brain": "CNS:",
    "Breast": "BR:",
    "Colon": "CO:",
    "Kidney": "RE:",
    "Leukemia": "LE:",
    "Lung": "LC:",
    "Melanoma": "ME:",
    "Ovarian": "OV:",
}


def load_config() -> dict:
    """Load the fixed analysis parameters recorded for the manuscript."""
    with CONFIG.open(encoding="utf-8") as handle:
        return json.load(handle)


def canonical_key(value: object) -> str:
    """Normalize a drug name for matching without changing its display name."""
    return "".join(character for character in str(value).lower() if character.isalnum())


def order_pairs(frame: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with drug names in deterministic alphabetical order."""
    result = frame.copy()
    swap = result["drug1"].astype(str) > result["drug2"].astype(str)
    left = result.loc[swap, "drug1"].copy()
    result.loc[swap, "drug1"] = result.loc[swap, "drug2"].to_numpy()
    result.loc[swap, "drug2"] = left.to_numpy()
    return result


def cancer_cells(columns: list[str], cancer: str) -> list[str]:
    """Return NCI-60 cell-line columns assigned to one manuscript cancer group."""
    cells = [column for column in columns if column.startswith(CANCER_PREFIX[cancer])]
    # Preserve the panel assignment used in the paper analysis.
    if cancer == "Breast":
        cells.append("ME:MDA-MB-435")
    elif cancer == "Melanoma" and "ME:MDA-MB-435" in cells:
        cells.remove("ME:MDA-MB-435")
    return cells


def load_outcomes() -> pd.DataFrame:
    """Load the final cancer-pair cohort and continuous consensus outcomes."""
    frame = pd.read_csv(INPUTS / "cancer_pair_consensus_outcomes.csv")
    frame = frame.rename(columns={"consensus_percentile": "continuous"})
    return order_pairs(frame)
