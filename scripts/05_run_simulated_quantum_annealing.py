#!/usr/bin/env python3
"""Test the three-input QUBO with OpenJij SQA."""

from __future__ import annotations

import argparse
import itertools
import sys
import time
from pathlib import Path

import numpy as np
import openjij as oj
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from drug_synergy.io import load_config  # noqa: E402
from drug_synergy.paths import RESULTS  # noqa: E402


CANCER_ORDER = {
    "Brain": 0,
    "Breast": 1,
    "Colon": 2,
    "Kidney": 3,
    "Leukemia": 4,
    "Lung": 5,
    "Melanoma": 6,
    "Ovarian": 7,
}


def build_problem(
    frame: pd.DataFrame,
    cardinality_penalty: float,
    forbidden_energy: float,
) -> tuple[list[str], dict[tuple[int, int], float], dict[tuple[int, int], float]]:
    """Construct the pairwise QUBO with a two-drug cardinality penalty."""
    drugs = sorted(set(frame["drug1"]) | set(frame["drug2"]))
    index = {drug: position for position, drug in enumerate(drugs)}
    energies = {
        tuple(sorted((index[row.drug1], index[row.drug2]))): float(row.qubo_energy)
        for row in frame.itertuples(index=False)
    }
    qubo = {(i, i): -3.0 * cardinality_penalty for i in range(len(drugs))}
    for i, j in itertools.combinations(range(len(drugs)), 2):
        objective = energies.get((i, j), forbidden_energy)
        qubo[(i, j)] = objective + 2.0 * cardinality_penalty
    return drugs, energies, qubo


def decode_samples(response, energies):
    """Collect valid two-drug samples and their occurrence counts."""
    found, counts, valid_reads = {}, {}, 0
    for record in response.data(["sample", "energy", "num_occurrences"]):
        chosen = tuple(sorted(i for i, value in record.sample.items() if value == 1))
        if len(chosen) != 2 or chosen not in energies:
            continue
        valid_reads += int(record.num_occurrences)
        found[chosen] = energies[chosen]
        counts[chosen] = counts.get(chosen, 0) + int(record.num_occurrences)
    return found, counts, valid_reads


def run(selected_cancers: list[str] | None = None) -> None:
    config = load_config()
    source = RESULTS / "primary_qubo" / "ranked_pairs.csv"
    if not source.exists():
        raise FileNotFoundError("Run 01_run_primary_qubo.py first")
    data = pd.read_csv(source)
    if selected_cancers:
        data = data[data["cancer"].isin(selected_cancers)]
    output = RESULTS / "sqa"
    output.mkdir(parents=True, exist_ok=True)
    run_rows, summary_rows, candidate_rows = [], [], []

    for cancer, frame in data.groupby("cancer", sort=True):
        frame = frame.sort_values("qubo_rank").reset_index(drop=True)
        drugs, energies, qubo = build_problem(
            frame,
            config["sqa_cardinality_penalty"],
            config["sqa_forbidden_pair_energy"],
        )
        exact_top = [
            tuple(sorted((drugs.index(row.drug1), drugs.index(row.drug2))))
            for row in frame.head(5).itertuples(index=False)
        ]
        exact_best = exact_top[0]
        all_found, all_counts = {}, {}

        for repeat in range(1, config["sqa_repeats"] + 1):
            sampler = oj.SQASampler()
            started = time.perf_counter()
            found, counts, valid_reads = {}, {}, 0
            first_seed = (
                config["sqa_seed"]
                + 1_000_000 * CANCER_ORDER[cancer]
                + 10_000 * repeat
            )
            for read_index in range(config["sqa_reads_per_repeat"]):
                response = sampler.sample_qubo(
                    qubo,
                    num_reads=1,
                    num_sweeps=config["sqa_sweeps"],
                    trotter=config["sqa_trotter"],
                    beta=config["sqa_beta"],
                    gamma=config["sqa_gamma"],
                    seed=first_seed + read_index,
                )
                one_found, one_counts, one_valid = decode_samples(response, energies)
                found.update(one_found)
                valid_reads += one_valid
                for pair, count in one_counts.items():
                    counts[pair] = counts.get(pair, 0) + count
            all_found.update(found)
            for pair, count in counts.items():
                all_counts[pair] = all_counts.get(pair, 0) + count
            ranked = sorted(found, key=found.get)
            best = ranked[0] if ranked else None
            run_rows.append({
                "cancer": cancer,
                "repeat": repeat,
                "seed": first_seed,
                "valid_read_fraction": valid_reads / config["sqa_reads_per_repeat"],
                "exact_optimum_recovered": int(exact_best in found),
                "best_pair_matches_exact": int(best == exact_best),
                "best_energy_gap": np.nan if best is None else found[best] - energies[exact_best],
                "top5_overlap": len(set(ranked[:5]) & set(exact_top)),
                "runtime_seconds": time.perf_counter() - started,
            })

        aggregate_top = sorted(all_found, key=all_found.get)[:5]
        total_reads = config["sqa_repeats"] * config["sqa_reads_per_repeat"]
        for rank, pair in enumerate(aggregate_top, start=1):
            candidate_rows.append({
                "cancer": cancer,
                "sampled_energy_rank": rank,
                "drug1": drugs[pair[0]],
                "drug2": drugs[pair[1]],
                "qubo_energy": energies[pair],
                "energy_gap_from_exact": energies[pair] - energies[exact_best],
                "sample_count": all_counts.get(pair, 0),
                "sample_frequency_all_reads": all_counts.get(pair, 0) / total_reads,
                "is_exact_optimum": int(pair == exact_best),
                "is_in_exact_top5": int(pair in exact_top),
            })
        cancer_runs = [row for row in run_rows if row["cancer"] == cancer]
        summary_rows.append({
            "cancer": cancer,
            "exact_best_pair": f"{drugs[exact_best[0]]} + {drugs[exact_best[1]]}",
            "exact_optimum_recovery_frequency": np.mean(
                [row["exact_optimum_recovered"] for row in cancer_runs]
            ),
            "exact_optimum_read_frequency": all_counts.get(exact_best, 0) / total_reads,
            "aggregate_top5_overlap": len(set(aggregate_top) & set(exact_top)),
            "mean_energy_gap": np.nanmean([row["best_energy_gap"] for row in cancer_runs]),
            "mean_runtime_seconds": np.mean([row["runtime_seconds"] for row in cancer_runs]),
        })
        print(f"Completed {cancer}", flush=True)

    pd.DataFrame(run_rows).to_csv(output / "sqa_independent_runs.csv", index=False)
    pd.DataFrame(summary_rows).to_csv(output / "sqa_summary_by_cancer.csv", index=False)
    pd.DataFrame(candidate_rows).to_csv(output / "sqa_low_energy_candidates.csv", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cancer", nargs="*", default=None)
    arguments = parser.parse_args()
    run(arguments.cancer)
