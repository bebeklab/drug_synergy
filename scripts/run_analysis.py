#!/usr/bin/env python3
"""Run the main analyses from the prepared input tables."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = (
    "01_run_primary_qubo.py",
    "02_run_retrieval_and_ablation.py",
    "03_run_supervised_comparison.py",
    "04_run_binding_sensitivity.py",
)


def main() -> None:
    for script in SCRIPTS:
        print(f"\n==> {script}", flush=True)
        subprocess.run([sys.executable, str(ROOT / "scripts" / script)], check=True)
    print(
        "\nCore analysis complete. Run 05_run_simulated_quantum_annealing.py "
        "separately because it is computationally intensive."
    )


if __name__ == "__main__":
    main()
