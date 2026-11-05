#!/usr/bin/env python3
"""Build pairwise mechanism-redundancy costs from reviewed drug annotations.

The ``almanac_drug`` column identifies the DrugComb/NCI-ALMANAC drugs used in
the combination screen. Target, action-type, and mechanism-class annotations
come primarily from ChEMBL 35, with reviewed additions from the NCI Drug
Dictionary and PubChem. The source used for each drug is retained in the
``annotation_source`` column of the input table.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from drug_synergy.mechanism import build_mechanism_pair_table  # noqa: E402


DEFAULT_ANNOTATIONS = ROOT / "data" / "annotations" / "mechanism_annotations.csv"
DEFAULT_OUTPUT = (
    ROOT / "results" / "generated" / "mechanism_rebuild"
    / "mechanism_redundancy_pair_costs.csv"
)
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotations", type=Path, default=DEFAULT_ANNOTATIONS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    annotations = pd.read_csv(args.annotations)
    pair_scores = build_mechanism_pair_table(annotations)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    pair_scores.to_csv(args.output, index=False)
    print(f"Wrote {len(pair_scores):,} pair scores to {args.output}")


if __name__ == "__main__":
    main()
