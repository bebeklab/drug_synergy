#!/usr/bin/env python3
"""Evaluate early retrieval and leave-one-input-out QUBO ablations."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from drug_synergy.metrics import (  # noqa: E402
    classification_metrics,
    retrieval_metrics,
)
from drug_synergy.paths import RESULTS  # noqa: E402


CONFIGURATIONS = {
    "Complete model": (
        "gi50_percentile_cost",
        "pathway_percentile_cost",
        "mechanism_percentile_cost",
    ),
    "Without pathway": ("gi50_percentile_cost", "mechanism_percentile_cost"),
    "Without GI50": ("pathway_percentile_cost", "mechanism_percentile_cost"),
    "Without mechanism": ("gi50_percentile_cost", "pathway_percentile_cost"),
}


def main() -> None:
    source = RESULTS / "primary_qubo" / "ranked_pairs.csv"
    if not source.exists():
        raise FileNotFoundError("Run 01_run_primary_qubo.py first")
    output = RESULTS / "evaluation"
    output.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(source)

    retrieval_rows = []
    for cancer, frame in data.groupby("cancer", sort=True):
        labels = frame["top1_label"].to_numpy(int)
        scores = -frame["qubo_energy"].to_numpy(float)
        row = {"cancer": cancer}
        for fraction in (0.01, 0.05, 0.10, 0.20):
            precision, recall, enrichment = retrieval_metrics(labels, scores, fraction)
            suffix = int(100 * fraction)
            row[f"precision_at_{suffix}_percent"] = precision
            row[f"recall_at_{suffix}_percent"] = recall
            row[f"ef_at_{suffix}_percent"] = enrichment
        retrieval_rows.append(row)

    ablation_rows = []
    for name, columns in CONFIGURATIONS.items():
        cancer_rows = []
        for _, frame in data.groupby("cancer", sort=True):
            labels = frame["top1_label"].to_numpy(int)
            scores = -frame[list(columns)].mean(axis=1).to_numpy(float)
            metrics = classification_metrics(labels, scores)
            _, recall, enrichment = retrieval_metrics(labels, scores, 0.10)
            cancer_rows.append((metrics["auroc"], metrics["auprc"], recall, enrichment))
        averages = np.mean(cancer_rows, axis=0)
        ablation_rows.append({
            "configuration": name,
            "auroc": averages[0],
            "auprc": averages[1],
            "recall_at_10_percent": averages[2],
            "ef_at_10_percent": averages[3],
        })

    pd.DataFrame(retrieval_rows).to_csv(output / "early_retrieval.csv", index=False)
    pd.DataFrame(ablation_rows).to_csv(output / "input_ablation.csv", index=False)
    print(pd.DataFrame(ablation_rows).round(4).to_string(index=False))


if __name__ == "__main__":
    main()
