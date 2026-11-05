# Analysis map

Each script has one role in the analysis workflow.

| Analysis | Script | Main output |
|---|---|---|
| Consensus-synergy outcomes | `scripts/00_build_consensus_outcomes.py` | `results/generated/consensus_outcomes/` |
| Mechanism-redundancy scores | `scripts/00_build_mechanism_costs.py` | `results/generated/mechanism_rebuild/` |
| Pathway-module scores | `scripts/build_pathway_modules.R` | `results/generated/pathway_rebuild/` |
| PFA binding representation | `scripts/00_fit_pfa_binding.py` | `results/generated/pfa_binding/` |
| Three-input QUBO | `scripts/01_run_primary_qubo.py` | `results/generated/primary_qubo/` |
| Retrieval and ablation | `scripts/02_run_retrieval_and_ablation.py` | `results/generated/evaluation/` |
| RF and neural-network comparison | `scripts/03_run_supervised_comparison.py` | `results/generated/ml_comparison/` |
| Binding sensitivity | `scripts/04_run_binding_sensitivity.py` | `results/generated/binding_sensitivity/` |
| Simulated quantum annealing | `scripts/05_run_simulated_quantum_annealing.py` | `results/generated/sqa/` |

Shared data loading, GI50 processing, QUBO scoring, and metric calculations are
in `src/drug_synergy/`. `scripts/run_analysis.py` runs the primary QUBO,
retrieval and ablation, supervised comparison, and binding sensitivity in
order. PFA fitting and SQA are separate because they take longer.
