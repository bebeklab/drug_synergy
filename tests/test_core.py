"""Small deterministic tests for the central scoring calculations."""

import numpy as np
import pandas as pd

from drug_synergy.mechanism import mechanism_pair_score
from drug_synergy.metrics import retrieval_metrics
from drug_synergy.qubo import score_equal_weight_qubo


def test_equal_weight_energy_is_mean_of_percentile_costs():
    frame = pd.DataFrame({
        "gi50_percentile_cost": [0.1, 0.9],
        "pathway_percentile_cost": [0.2, 0.6],
        "mechanism_percentile_cost": [0.3, 0.3],
    })
    scored = score_equal_weight_qubo(frame)
    np.testing.assert_allclose(scored["qubo_energy"], [0.2, 0.6])
    assert scored["qubo_rank"].tolist() == [1, 2]


def test_retrieval_enrichment_uses_prevalence_baseline():
    labels = np.array([1, 0, 0, 0])
    scores = np.array([1.0, 0.8, 0.2, 0.1])
    precision, recall, enrichment = retrieval_metrics(labels, scores, 0.25)
    assert precision == 1.0
    assert recall == 1.0
    assert enrichment == 4.0


def test_mechanism_score_uses_class_when_a_target_set_is_absent():
    left = pd.Series({
        "almanac_drug": "Drug A",
        "primary_targets": "Target 1",
        "mechanism_class": "DNA crosslinker",
    })
    right = pd.Series({
        "almanac_drug": "Drug B",
        "primary_targets": np.nan,
        "mechanism_class": "DNA crosslinker",
    })
    score = mechanism_pair_score(left, right)
    assert np.isnan(score["target_jaccard"])
    assert score["shared_mechanism_class"] == 1.0
    assert score["mechanism_redundancy"] == 1.0
