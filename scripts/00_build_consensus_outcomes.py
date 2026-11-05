#!/usr/bin/env python3
"""Build cancer-specific consensus-synergy outcomes from NCI-ALMANAC records."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "data" / "raw" / "drugcomb" / "summary_df_ALMANAC.csv"
DEFAULT_OUTPUT = (
    ROOT / "results" / "generated" / "consensus_outcomes"
    / "cancer_pair_consensus_outcomes.csv"
)
SYNERGY_COLUMNS = [
    "synergy_bliss",
    "synergy_loewe",
    "synergy_zip",
    "synergy_hsa",
]
TISSUE_NAMES = {
    "brain": "Brain",
    "breast": "Breast",
    "large_intestine": "Colon",
    "kidney": "Kidney",
    "haematopoietic_and_lymphoid": "Leukemia",
    "lung": "Lung",
    "skin": "Melanoma",
    "ovary": "Ovarian",
}


def normalized_name(value: object) -> str:
    """Normalize spacing and case for drug-name ordering and tissue matching."""
    return " ".join(str(value).lower().split())


def ordered_pair(left: object, right: object) -> tuple[str, str]:
    """Return one order-invariant representation of a drug pair."""
    first, second = str(left), str(right)
    if normalized_name(first) <= normalized_name(second):
        return first, second
    return second, first


def build_consensus_outcomes(data: pd.DataFrame) -> pd.DataFrame:
    """Aggregate four synergy definitions into cancer-specific pair outcomes."""
    required = {
        "drug_row", "drug_col", "cell_line_name", "tissue_name",
        *SYNERGY_COLUMNS,
    }
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"ALMANAC table is missing columns: {sorted(missing)}")

    data = data.dropna(
        subset=["drug_row", "drug_col", "cell_line_name", "tissue_name"]
    ).copy()
    for metric in SYNERGY_COLUMNS:
        data[metric] = pd.to_numeric(data[metric], errors="coerce")

    tissue_keys = data["tissue_name"].map(
        lambda value: normalized_name(value).replace(" ", "_")
    )
    data["cancer_type"] = tissue_keys.map(TISSUE_NAMES)
    data = data.loc[data["cancer_type"].notna()].copy()

    percentile_columns = []
    for metric in SYNERGY_COLUMNS:
        column = f"{metric}_percentile"
        data[column] = data.groupby("cancer_type")[metric].rank(
            method="average", pct=True
        )
        percentile_columns.append(column)

    data["n_metrics"] = data[percentile_columns].notna().sum(axis=1)
    data["screen_consensus"] = data[percentile_columns].mean(axis=1, skipna=True)
    data = data.loc[data["n_metrics"] >= 3].copy()

    pairs = data.apply(
        lambda row: ordered_pair(row["drug_row"], row["drug_col"]), axis=1
    )
    data[["pair_a", "pair_b"]] = pd.DataFrame(pairs.tolist(), index=data.index)

    # Collapse repeated screens within each cell line before averaging across
    # cell lines, so cell lines with more records do not receive more weight.
    by_cell = (
        data.groupby(
            ["cancer_type", "pair_a", "pair_b", "cell_line_name"],
            as_index=False,
        )
        .agg(
            cell_line_consensus=("screen_consensus", "median"),
            n_screens=("screen_consensus", "size"),
        )
    )
    scores = (
        by_cell.groupby(["cancer_type", "pair_a", "pair_b"], as_index=False)
        .agg(
            consensus_score=("cell_line_consensus", "mean"),
            n_cell_lines=("cell_line_name", "nunique"),
            n_screens=("n_screens", "sum"),
        )
    )
    scores["percentile_within_cancer"] = scores.groupby("cancer_type")[
        "consensus_score"
    ].rank(method="average", pct=True)
    scores["top1_flag"] = (
        scores["percentile_within_cancer"] >= 0.99
    ).astype(int)
    return scores


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    columns = [
        "drug_row", "drug_col", "cell_line_name", "tissue_name",
        *SYNERGY_COLUMNS,
    ]
    records = pd.read_csv(args.input, usecols=columns, na_values=["\\N"])
    scores = build_consensus_outcomes(records)
    analysis_table = scores.rename(columns={
        "cancer_type": "cancer",
        "pair_a": "drug1",
        "pair_b": "drug2",
        "percentile_within_cancer": "consensus_percentile",
    })[["cancer", "drug1", "drug2", "consensus_percentile"]]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    analysis_table.to_csv(args.output, index=False)
    print(
        f"Wrote {len(analysis_table):,} cancer-specific pair outcomes "
        f"to {args.output}"
    )


if __name__ == "__main__":
    main()
