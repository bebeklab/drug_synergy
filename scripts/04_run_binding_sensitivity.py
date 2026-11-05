#!/usr/bin/env python3
"""Compare primary, binding-inclusive, and PFA-binding-only rankings."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from drug_synergy.metrics import classification_metrics  # noqa: E402
from drug_synergy.paths import INPUTS, RESULTS  # noqa: E402


def pair_value(matrix: pd.DataFrame, left: str, right: str) -> float:
    return float(matrix.loc[left, right])


def main() -> None:
    source = RESULTS / "primary_qubo" / "ranked_pairs.csv"
    if not source.exists():
        raise FileNotFoundError("Run 01_run_primary_qubo.py first")
    output = RESULTS / "binding_sensitivity"
    output.mkdir(parents=True, exist_ok=True)
    pairs = pd.read_csv(source)
    similarity = pd.read_csv(INPUTS / "pfa2_binding_similarity.csv", index_col=0)
    mask = pd.read_csv(INPUTS / "pfa2_binding_pair_mask.csv", index_col=0)

    known = set(similarity.index) & set(similarity.columns)
    pairs = pairs[pairs["drug1"].isin(known) & pairs["drug2"].isin(known)].copy()
    pairs["binding_supported"] = [
        pair_value(mask, left, right) == 1
        for left, right in zip(pairs["drug1"], pairs["drug2"])
    ]
    pairs = pairs[pairs["binding_supported"]].copy()
    pairs["binding_similarity_cost"] = [
        pair_value(similarity, left, right)
        for left, right in zip(pairs["drug1"], pairs["drug2"])
    ]

    ranked_groups = []
    metric_rows = []
    for cancer, raw in pairs.groupby("cancer", sort=True):
        frame = raw.reset_index(drop=True).copy()
        raw_to_ranked = {
            "gi50_activity_cost": "gi50_p",
            "structured_pathway_cost": "pathway_p",
            "mechanism_redundancy": "mechanism_p",
            "binding_similarity_cost": "binding_p",
        }
        for raw_name, ranked_name in raw_to_ranked.items():
            frame[ranked_name] = frame[raw_name].rank(method="average", pct=True)

        frame["energy_primary_same_cohort"] = frame[
            ["gi50_p", "pathway_p", "mechanism_p"]
        ].mean(axis=1)
        frame["energy_binding_inclusive"] = frame[
            ["gi50_p", "pathway_p", "mechanism_p", "binding_p"]
        ].mean(axis=1)
        labels = (frame["continuous"].to_numpy(float) >= 0.99).astype(int)
        score_columns = {
            "Primary QUBO": "energy_primary_same_cohort",
            "Binding-inclusive QUBO": "energy_binding_inclusive",
            "PFA binding only": "binding_p",
        }
        for model, energy_column in score_columns.items():
            metrics = classification_metrics(
                labels, -frame[energy_column].to_numpy(float)
            )
            metric_rows.append({
                "cancer": cancer,
                "configuration": model,
                "n_pairs": len(frame),
                "n_positive": int(labels.sum()),
                **metrics,
            })
        frame["top1_label"] = labels
        ranked_groups.append(frame)

    metrics = pd.DataFrame(metric_rows)
    macro = metrics.groupby("configuration", sort=False).agg(
        auroc_macro_mean=("auroc", "mean"),
        auprc_macro_mean=("auprc", "mean"),
        prevalence_macro_mean=("prevalence", "mean"),
    ).reset_index()
    pd.concat(ranked_groups, ignore_index=True).to_csv(
        output / "ranked_pairs.csv", index=False
    )
    metrics.to_csv(output / "cancer_specific_metrics.csv", index=False)
    macro.to_csv(output / "macro_metrics.csv", index=False)
    print(macro.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
