#!/usr/bin/env python3
"""Generate the three-input label-free QUBO ranking and metrics."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from drug_synergy.features import (  # noqa: E402
    add_gi50_activity_cost,
    load_primary_pair_inputs,
    percentile_rank_costs,
)
from drug_synergy.metrics import classification_metrics  # noqa: E402
from drug_synergy.paths import RESULTS  # noqa: E402
from drug_synergy.qubo import score_equal_weight_qubo  # noqa: E402


def main() -> None:
    output = RESULTS / "primary_qubo"
    output.mkdir(parents=True, exist_ok=True)
    data = load_primary_pair_inputs()
    ranked_groups = []
    metric_rows = []

    for cancer, raw in data.groupby("cancer", sort=True):
        frame = add_gi50_activity_cost(raw.reset_index(drop=True), cancer)
        frame = percentile_rank_costs(frame)
        frame = score_equal_weight_qubo(frame)
        frame["top1_label"] = (frame["continuous"].astype(float) >= 0.99).astype(int)
        labels = frame["top1_label"].to_numpy(int)
        scores = -frame["qubo_energy"].to_numpy(float)
        metrics = classification_metrics(labels, scores)
        metric_rows.append({
            "cancer": cancer,
            "n_drugs": len(set(frame["drug1"]) | set(frame["drug2"])),
            "n_pairs": len(frame),
            "n_positive": int(labels.sum()),
            **metrics,
        })
        ranked_groups.append(frame)

    ranked = pd.concat(ranked_groups, ignore_index=True)
    metrics = pd.DataFrame(metric_rows)
    macro = pd.DataFrame([{
        "cancers": metrics["cancer"].nunique(),
        "unique_evaluable_drugs": len(set(data["drug1"]) | set(data["drug2"])),
        "auroc_macro_mean": metrics["auroc"].mean(),
        "auprc_macro_mean": metrics["auprc"].mean(),
        "prevalence_macro_mean": metrics["prevalence"].mean(),
    }])
    ranked.to_csv(output / "ranked_pairs.csv", index=False)
    metrics.to_csv(output / "cancer_specific_metrics.csv", index=False)
    macro.to_csv(output / "macro_metrics.csv", index=False)
    print(metrics.round(4).to_string(index=False))
    print("\nMacro means")
    print(macro.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
