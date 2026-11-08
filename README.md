# Quantum-Inspired Optimization of Cancer-Specific Drug Combinations Without Combination-Response Training Data

The analysis ranks cancer-specific drug pairs from three single-drug inputs:

- GI50 activity in the relevant NCI-60 cancer panel;
- pathway-response modules; and
- mechanism-of-action redundancy.

Each input is converted to a within-cancer percentile cost and receives equal
weight in the QUBO. Combination-response data are used to evaluate the ranking
and train the supervised comparison models. 

## Folder guide

```text
config/             Seeds, thresholds, and model settings
data/inputs/        Tables read by the analysis scripts
data/annotations/   Reviewed drug names, targets, and mechanisms
docs/               Data sources and analysis map
scripts/            Analysis steps in execution order
src/drug_synergy/   Shared loading, scoring, and metric functions
tests/              Checks for the central calculations
```

Run the analysis through `scripts/`. Shared loading, scoring, and evaluation
functions are in `src/`. The scripts create `results/generated/` locally when
needed. That output directory is not stored in GitHub.

## Install

Python 3.11 or newer is recommended.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[sqa,test]'
```

## Prepare the analysis inputs

The files in `data/inputs/` are inputs to the main QUBO and evaluation scripts;
they are not raw database downloads. Several are created by the preparation
scripts described below:

```text
included data/inputs -> main analysis -> local outputs
source data -> preparation scripts -> rebuilt tables for inspection
```

Prepared copies are included so the main analysis runs without downloading
large databases or repeating manual drug-identifier review. Use the included
files to reproduce the reported analysis, then continue to **Run the main
analysis** below.

The following sections explain how the consensus, pathway, mechanism, and
binding tables were created. Run them if you want to inspect those preparation
steps or create new copies from the source data. Rebuilt tables are written
under `results/generated/` for inspection and do not overwrite the included
files in `data/inputs/`.

### Create the consensus-synergy outcomes

Download [`summary_df_ALMANAC.csv`](https://zenodo.org/api/records/23986100/files/summary_df_ALMANAC.csv/content)
from the [project dataset](https://doi.org/10.5281/zenodo.23986099). The file
contains the NCI-ALMANAC combination records used to define the observed
synergy outcome. Place it at
`data/raw/drugcomb/summary_df_ALMANAC.csv`, then run:

```bash
python scripts/00_build_consensus_outcomes.py
```

The script ranks Bliss, Loewe, ZIP, and HSA values separately within each
cancer. For a screen with at least three values, it averages their percentile
ranks. It then takes the median across repeated screens in the same cell line
and the mean across cell lines. The resulting score describes the observed
synergy of each drug pair within one cancer group. The script saves these
scores to
`results/generated/consensus_outcomes/cancer_pair_consensus_outcomes.csv` in
the same four-column format read by the main analysis.

### Create the pathway-module input

The pathway-module cost is created from the ALMANAC summary records,
GSE116450 expression matrix, and the KEGG collection downloaded through
MSigDB. Install the R packages once:

```bash
Rscript scripts/install_r_dependencies.R
```

Place the source files in the locations listed in
`docs/RAW_DATA_DOWNLOADS.md`. Download the public inputs with:

```bash
python scripts/download_pathway_sources.py
```

MSigDB requires login, so place
`c2.cp.kegg.v2023.1.Hs.symbols.gmt` in `data/raw/msigdb/` manually. Then run:

```bash
Rscript scripts/build_pathway_modules.R
```

The script uses untreated (`0nM`), 24-hour samples from GSE116450 to compare
gene expression between sensitive and resistant cell lines. It performs KEGG
pathway enrichment, groups drugs into Walktrap modules, calculates pathway
costs for each drug pair, and saves the results to
`results/generated/pathway_rebuild/`. Intermediate gene-level tables are not
saved.

Pathway enrichment uses 10,000 permutations and may take several minutes.
Detailed execution settings and notes about rebuilding this input are provided
in `docs/RAW_DATA_DOWNLOADS.md`.

### Create the mechanism-redundancy input

The DrugComb/NCI-ALMANAC drug names and their reviewed target and
mechanism-class annotations are stored in
`data/annotations/mechanism_annotations.csv`. ChEMBL 35 supplied most
annotations; the NCI Drug Dictionary and PubChem were used for reviewed
additions. The source for each drug is recorded in the table. Compute the
pairwise target overlap and mechanism-redundancy scores with:

```bash
python scripts/00_build_mechanism_costs.py
```

For each drug pair, the script measures shared protein targets and checks
whether the drugs share a mechanism class. It combines those two values into
the mechanism-redundancy score described in the paper. The script saves the
new table in `results/generated/mechanism_rebuild/`. The script reads the local
reviewed annotations; it does not download new records or change that
annotation file.

### Create the binding representation

Binding was tested as an additional input, but it is not part of the primary
three-input QUBO. Fit the two-dimensional PFA representation with:

```bash
python scripts/00_fit_pfa_binding.py
```

The script starts with the harmonized IC50, Ki, and Kd table and saves the PFA
parameters and drug-pair binding similarities in
`results/generated/pfa_binding/`. The binding sensitivity analysis reads the
saved similarity matrix in `data/inputs/`, which is the matrix used for the
paper.

## Run the main analysis

After installing the software, run the analysis from the repository root:

```bash
python scripts/run_analysis.py
```

This command reads the prepared input tables and runs the three-input QUBO,
early-retrieval and ablation analyses, five-fold supervised comparisons, and
the binding sensitivity analysis. It creates `results/generated/` and its
subdirectories automatically.

To run the analysis steps separately:

```bash
python scripts/01_run_primary_qubo.py
python scripts/02_run_retrieval_and_ablation.py
python scripts/03_run_supervised_comparison.py
python scripts/04_run_binding_sensitivity.py
```

## Simulated quantum annealing

The SQA step is run separately because it takes longer:

```bash
python scripts/05_run_simulated_quantum_annealing.py
```

To check the setup with one cancer:

```bash
python scripts/05_run_simulated_quantum_annealing.py --cancer Brain
```

The paper used 10 batches per cancer, 500 reads per batch, 1,000 sweeps per
read, and seed `20261030`. These settings are in `config/analysis.json`.

## Tests

```bash
pytest
```

The tests check the equal-weight QUBO calculation, retrieval metrics, and the
mechanism-redundancy calculation.

## Data sources

Inputs come from NCI-ALMANAC/DrugComb, CellMiner, GSE116450, KEGG/MSigDB,
ChEMBL, BindingDB, the NCI Drug Dictionary, and PubChem. The processed analysis
files are included. Creating every input from its source requires separate
downloads and some manual identifier review. Source links and versions are in
`docs/RAW_DATA_DOWNLOADS.md`; `data/README.md` maps each input file to its
source.

The source ALMANAC table used by this project is archived on Zenodo:
[doi:10.5281/zenodo.23986099](https://doi.org/10.5281/zenodo.23986099).

## Use restrictions

This repository is provided solely for anonymous peer review. Reuse,
modification, redistribution, and publication are not permitted. See `LICENSE`
for the complete notice.
