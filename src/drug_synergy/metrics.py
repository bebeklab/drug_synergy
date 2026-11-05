"""Evaluation metrics used in the manuscript."""

from __future__ import annotations

import numpy as np
from sklearn.metrics import auc, precision_recall_curve, roc_auc_score


def trapezoidal_auprc(labels: np.ndarray, scores: np.ndarray) -> float:
    """Calculate area under the precision-recall curve by trapezoidal integration."""
    precision, recall, _ = precision_recall_curve(labels, scores)
    return float(auc(recall[::-1], precision[::-1]))


def classification_metrics(labels: np.ndarray, scores: np.ndarray) -> dict[str, float]:
    """Return AUROC, trapezoidal AUPRC, and positive prevalence."""
    return {
        "prevalence": float(labels.mean()),
        "auroc": float(roc_auc_score(labels, scores)),
        "auprc": trapezoidal_auprc(labels, scores),
    }


def retrieval_metrics(
    labels: np.ndarray, scores: np.ndarray, fraction: float
) -> tuple[float, float, float]:
    """Return precision, recall, and enrichment in the top-ranked fraction."""
    count = max(1, int(np.ceil(len(labels) * fraction)))
    selected = np.argsort(-scores)[:count]
    true_positives = int(labels[selected].sum())
    precision = true_positives / count
    recall = true_positives / int(labels.sum())
    enrichment = precision / float(labels.mean())
    return precision, recall, enrichment
