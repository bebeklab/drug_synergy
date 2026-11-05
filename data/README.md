# Data files

The files in `inputs/` are the prepared tables read by the main analysis. They
are called inputs because of their role in that analysis, even when a separate
preparation script creates them from raw data. Large source downloads are not
stored in the GitHub repository.

The primary analysis reads four files: the observed cancer-pair outcomes, the
GI50 matrix, the pathway-module pair costs, and the mechanism-redundancy pair
costs. The additional binding files support the separate binding analysis.
Intermediate tables created by the preparation scripts are written to
`results/generated/` and are not stored here.

| File | What it contains | Source |
|---|---|---|
| `inputs/cancer_pair_consensus_outcomes.csv` | Cancer, drug pair, and observed consensus-synergy percentile for the 30,208 evaluated pairs. | Created from `summary_df_ALMANAC.csv` by `scripts/00_build_consensus_outcomes.py` |
| `inputs/nci60_gi50_matrix.csv` | Single-drug GI50 profiles across NCI-60 cell lines. | CellMiner/NCI-60 |
| `inputs/pathway_module_pair_costs.csv` | Cancer-specific pathway-module costs for each pair. | NCI-60 response, untreated 0-nM 24-hour GSE116450 expression, and KEGG/MSigDB |
| `inputs/mechanism_redundancy_pair_costs.csv` | Target overlap, shared mechanism class, and the resulting redundancy score for every drug pair. | Created from the reviewed annotations by `scripts/00_build_mechanism_costs.py` |
| `inputs/binding_assays_raw_harmonized.tsv` | Harmonized IC50, Ki, and Kd records. | BindingDB May 2026 and ChEMBL 35 |
| `inputs/pfa2_binding_similarity.csv` | Pairwise similarity between PFA binding profiles. | Derived from the binding records |
| `inputs/pfa2_binding_pair_mask.csv` | Whether both drugs have a binding profile. | Derived from the binding records |

Included preparation scripts create the consensus outcomes, pathway-module
costs, mechanism-redundancy costs, and PFA binding files. The GI50 matrix and
reviewed annotation tables are included because their preparation involved
source-specific identifier matching or manual review.

`annotations/mechanism_annotations.csv` links the DrugComb/NCI-ALMANAC drug
names to the reviewed annotations used for the mechanism calculation. Its
`annotation_source` column records whether each annotation came from ChEMBL 35
or a reviewed external source.

`annotations/almanac_cellminer_nsc_crosswalk.csv` records the NSC identifier
selected for each ALMANAC drug, the corresponding CellMiner name, the matching
method, data coverage, and CellMiner quality-control status. This table
documents how the GI50 matrix was assembled from the CellMiner workbook.

Creating these files from their original sources requires NCI-ALMANAC/DrugComb records, CellMiner data,
GSE116450 expression, MSigDB `c2.cp.kegg.v2023.1.Hs`, ChEMBL 35, and the May
2026 BindingDB release. The mechanism analysis also uses
`annotations/mechanism_manual_overrides.csv`.

Download links, exact releases, expected local filenames, and access notes are
in [`docs/RAW_DATA_DOWNLOADS.md`](../docs/RAW_DATA_DOWNLOADS.md). Some sources
require a login, a web-form download, or manual identifier review. The main
workflow starts from the input files listed above.
