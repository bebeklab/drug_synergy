#!/usr/bin/env python3
"""Run five-fold pair-disjoint RF and neural-network regression comparisons."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import StratifiedKFold
from sklearn.neural_network import MLPRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from drug_synergy.io import load_config  # noqa: E402
from drug_synergy.metrics import classification_metrics  # noqa: E402
from drug_synergy.paths import RESULTS  # noqa: E402
from drug_synergy.qubo import PRIMARY_COSTS  # noqa: E402


def main() -> None:
    source = RESULTS / "primary_qubo" / "ranked_pairs.csv"
    if not source.exists():
        raise FileNotFoundError("Run 01_run_primary_qubo.py first")
    output = RESULTS / "ml_comparison"
    output.mkdir(parents=True, exist_ok=True)
    config = load_config()
    data = pd.read_csv(source)
    predictions = []

    for cancer_index, (cancer, raw) in enumerate(data.groupby("cancer", sort=True)):
        frame = raw.reset_index(drop=True)
        features = frame[list(PRIMARY_COSTS)].to_numpy(float)
        continuous = frame["continuous"].to_numpy(float)
        labels = frame["top1_label"].to_numpy(int)
        folds = StratifiedKFold(
            n_splits=config["ml_folds"],
            shuffle=True,
            random_state=config["ml_seed"] + cancer_index,
        )

        for fold, (train, test) in enumerate(folds.split(features, labels), start=1):
            forest = RandomForestRegressor(
                n_estimators=config["random_forest_trees"],
                criterion="squared_error",
                max_features=None,
                min_samples_leaf=1,
                n_jobs=-1,
                random_state=config["ml_seed"] + 100 * cancer_index + fold,
            )
            forest.fit(features[train], continuous[train])

            neural = make_pipeline(
                StandardScaler(),
                MLPRegressor(
                    hidden_layer_sizes=tuple(config["neural_hidden_layers"]),
                    activation="relu",
                    solver="sgd",
                    learning_rate_init=1e-3,
                    momentum=0.5,
                    alpha=1e-4,
                    batch_size=64,
                    max_iter=1000,
                    early_stopping=True,
                    validation_fraction=0.2,
                    n_iter_no_change=30,
                    random_state=config["ml_seed"] + 1000 + 100 * cancer_index + fold,
                ),
            )
            neural.fit(features[train], continuous[train])

            for model, scores in (
                ("Random-forest regression", forest.predict(features[test])),
                ("Neural-network regression", neural.predict(features[test])),
            ):
                for local_index, row_index in enumerate(test):
                    predictions.append({
                        "cancer": cancer,
                        "fold": fold,
                        "model": model,
                        "drug1": frame.iloc[row_index]["drug1"],
                        "drug2": frame.iloc[row_index]["drug2"],
                        "top1_label": int(labels[row_index]),
                        "continuous_outcome": float(continuous[row_index]),
                        "score": float(scores[local_index]),
                    })

    predictions = pd.DataFrame(predictions)
    metric_rows = []
    for (cancer, model), group in predictions.groupby(["cancer", "model"], sort=True):
        labels = group["top1_label"].to_numpy(int)
        scores = group["score"].to_numpy(float)
        metric_rows.append({
            "cancer": cancer,
            "model": model,
            "n_pairs": len(group),
            "n_positive": int(labels.sum()),
            **classification_metrics(labels, scores),
        })

    for cancer, frame in data.groupby("cancer", sort=True):
        labels = frame["top1_label"].to_numpy(int)
        metric_rows.append({
            "cancer": cancer,
            "model": "Equal-weight QUBO",
            "n_pairs": len(frame),
            "n_positive": int(labels.sum()),
            **classification_metrics(labels, -frame["qubo_energy"].to_numpy(float)),
        })

    metrics = pd.DataFrame(metric_rows).sort_values(["cancer", "model"])
    macro = metrics.groupby("model", sort=False).agg(
        auroc_mean=("auroc", "mean"),
        auprc_mean=("auprc", "mean"),
        prevalence_mean=("prevalence", "mean"),
    ).reset_index()
    predictions.to_csv(output / "out_of_fold_predictions.csv", index=False)
    metrics.to_csv(output / "cancer_specific_metrics.csv", index=False)
    macro.to_csv(output / "macro_metrics.csv", index=False)
    print(macro.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
