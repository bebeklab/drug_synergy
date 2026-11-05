"""Pairwise redundancy calculations for reviewed ALMANAC drug annotations.

This module calculates scores from a local annotation table. It does not query
DrugComb, ChEMBL, the NCI Drug Dictionary, or PubChem at run time.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd


def annotation_set(value: object) -> set[str]:
    """Split a semicolon-delimited annotation field into cleaned terms."""
    if pd.isna(value):
        return set()
    return {term.strip() for term in str(value).split(";") if term.strip()}


def mechanism_pair_score(left: pd.Series, right: pd.Series) -> dict[str, object]:
    """Calculate target overlap and mechanism redundancy for one drug pair.

    When both drugs have protein-target annotations, redundancy is the mean
    of target-set Jaccard overlap and the shared-class indicator. If either
    target set is absent, the score uses the shared-class indicator alone.
    """
    left_targets = annotation_set(left["primary_targets"])
    right_targets = annotation_set(right["primary_targets"])
    left_classes = annotation_set(left["mechanism_class"])
    right_classes = annotation_set(right["mechanism_class"])

    shared_targets = left_targets & right_targets
    shared_classes = left_classes & right_classes
    target_overlap_available = bool(left_targets and right_targets)
    annotation_complete = bool(left_classes and right_classes)

    if target_overlap_available:
        target_jaccard = len(shared_targets) / len(left_targets | right_targets)
    else:
        target_jaccard = np.nan

    shared_class = float(bool(shared_classes)) if annotation_complete else np.nan
    if not annotation_complete:
        redundancy = np.nan
    elif target_overlap_available:
        redundancy = 0.5 * target_jaccard + 0.5 * shared_class
    else:
        redundancy = shared_class

    return {
        "drug1": left["almanac_drug"],
        "drug2": right["almanac_drug"],
        "target_jaccard": target_jaccard,
        "shared_mechanism_class": shared_class,
        "mechanism_redundancy": redundancy,
        "annotation_complete": annotation_complete,
        "target_overlap_available": target_overlap_available,
        "shared_targets": "; ".join(sorted(shared_targets)),
        "shared_classes": "; ".join(sorted(shared_classes)),
    }


def build_mechanism_pair_table(annotations: pd.DataFrame) -> pd.DataFrame:
    """Calculate mechanism redundancy for every unordered annotated drug pair."""
    required = {"almanac_drug", "primary_targets", "mechanism_class"}
    missing = required - set(annotations.columns)
    if missing:
        raise ValueError(f"Annotation table is missing columns: {sorted(missing)}")
    if annotations["almanac_drug"].duplicated().any():
        duplicates = annotations.loc[
            annotations["almanac_drug"].duplicated(), "almanac_drug"
        ].tolist()
        raise ValueError(f"Duplicate ALMANAC drug names: {duplicates}")

    rows = [
        mechanism_pair_score(left, right)
        for (_, left), (_, right) in combinations(annotations.iterrows(), 2)
    ]
    return pd.DataFrame(rows)
