# Source data

The main analysis reads the prepared tables in `data/inputs/`. This page
explains where each table came from and what must be downloaded to create it
again. A newer database release may contain different drugs, targets, genes,
or pathway definitions and may therefore produce different results.

Create a local staging area before downloading anything:

```bash
mkdir -p data/raw/{drugcomb,cellminer,geo,msigdb,chembl,bindingdb}
```

`data/raw/` is ignored by Git because several source files are large and some
have redistribution terms set by their providers.

## Source checklist

| Source | Release or file used | Official access | Suggested local location |
|---|---|---|---|
| DrugComb / NCI-ALMANAC | `summary_df_ALMANAC.csv`, the project copy of the ALMANAC summary subset | [DrugComb v1.4 archive](https://zenodo.org/records/11102665) for the broader public release | `data/raw/drugcomb/summary_df_ALMANAC.csv` |
| CellMiner / NCI-60 | Compound activity data, including `AVGDATA` where available | [CellMiner download page](https://discover.nci.nih.gov/cellminer/loadDownload.do) | `data/raw/cellminer/` |
| GEO | `GSE116450_series_matrix.txt.gz` | [GSE116450 record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116450) | `data/raw/geo/GSE116450_series_matrix.txt.gz` |
| Ensembl | Release 115 Affymetrix HG-U133 Plus 2 to HGNC mapping | [Ensembl 115 archive](https://sep2025.archive.ensembl.org/) | `data/raw/ensembl/ensembl115_affy_hgnc.tsv` |
| MSigDB | `c2.cp.kegg.v2023.1.Hs.symbols.gmt` | [MSigDB downloads](https://www.gsea-msigdb.org/gsea/downloads.jsp) | `data/raw/msigdb/c2.cp.kegg.v2023.1.Hs.symbols.gmt` |
| ChEMBL | ChEMBL 35 SQLite release | [ChEMBL 35 archive](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_35/) | `data/raw/chembl/chembl_35_sqlite/chembl_35.db` |
| BindingDB | May 2026 all-data TSV, `BindingDB_All_202605_tsv.zip` | [BindingDB downloads](https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp) and [archived releases](https://library.ucsd.edu/dc/collection/bb03870458) | `data/raw/bindingdb/BindingDB_All_202605_tsv.zip` |
| NCI Drug Dictionary | Individual drug records used to review incomplete mechanism annotations | [NCI dictionaries](https://www.cancer.gov/publications/dictionaries) | No bulk file required |
| PubChem | Individual compound records used to review names, structures, and mechanisms | [PubChem](https://pubchem.ncbi.nlm.nih.gov/) | No bulk file required |

## Download steps

### 1. DrugComb outcomes

The analysis started from `summary_df_ALMANAC.csv`. It contains the
NCI-ALMANAC records selected from DrugComb, including Bliss, Loewe, ZIP, and
HSA scores. The file is about 291 MB, so it belongs in the accompanying Zenodo
data deposit rather than in GitHub. After downloading it, save it as
`data/raw/drugcomb/summary_df_ALMANAC.csv`.

Run `python scripts/00_build_consensus_outcomes.py` to create one observed
consensus-synergy score for each cancer and drug pair. The public DrugComb v1.4
archive contains a broader summary table. Using that broader file requires
selecting the same ALMANAC records before running the script.

### 2. CellMiner GI50 data

Open the [CellMiner download page](https://discover.nci.nih.gov/cellminer/loadDownload.do),
select the compound-activity data for the NCI-60 panel, and save the downloaded
files under `data/raw/cellminer/`. Match the drug records with the NSC
identifiers in `data/annotations/almanac_cellminer_nsc_crosswalk.csv`. Use the
aggregate `AVGDATA` profile when available;
otherwise, take the median of the available screening measurements for each
drug--cell-line combination.

The matrix read by the scripts is
`data/inputs/nci60_gi50_matrix.csv`. The current repository does not
automate the CellMiner web-form download or the initial NSC identifier review.

### 3. GSE116450 expression

Download the public pathway files that do not require login with:

```bash
python scripts/download_pathway_sources.py
```

The pathway workflow selects the untreated control samples labeled `0nM` at
24 hours from the expression matrix. The download script also saves the fixed
Ensembl 115 BioMart mapping used to convert HG-U133 Plus 2 probe identifiers to
HGNC symbols. The pathway build then runs locally from these versioned files.

### 4. KEGG pathways from MSigDB

MSigDB requires a free account. Sign in on the
[MSigDB downloads page](https://www.gsea-msigdb.org/gsea/downloads.jsp), choose
release 2023.1 for human gene symbols, and download
`c2.cp.kegg.v2023.1.Hs.symbols.gmt`. Do not replace it with the current KEGG
collection when reproducing the reported results.

The pathway pair costs read by the analysis are stored in
`data/inputs/pathway_module_pair_costs.csv`.

Create the module assignments and pair costs from the ALMANAC response, GEO
expression, fixed Ensembl probe map, and GMT download with:

```bash
Rscript scripts/install_r_dependencies.R
Rscript scripts/build_pathway_modules.R
```

The script uses seed 42 and one worker by default. The seed used to create
the pathway file for the paper was not recorded. As a result, running the
script may change pathways close to the significance threshold and may
produce different modules. Use the pathway files in `data/inputs/` when
running the analysis with the inputs reported in the paper.

The pathway workflow uses ALMANAC single-drug relative-inhibition values to
define sensitive and resistant cell lines. CellMiner GI50 is a separate QUBO
input and is not substituted for relative inhibition in this pathway step.

### 5. ChEMBL 35

Download the SQLite archive from the
[ChEMBL 35 release directory](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_35/),
extract it under `data/raw/chembl/`, and confirm that the database is available
at:

```text
data/raw/chembl/chembl_35_sqlite/chembl_35.db
```

ChEMBL 35 supplied most mechanism annotations and part of the binding-assay
data. The drug-level annotation table read by the mechanism script is
`data/annotations/mechanism_annotations.csv`; manual additions and their source
URLs are recorded in `data/annotations/mechanism_manual_overrides.csv`.

Run `python scripts/00_build_mechanism_costs.py` to compute the score for every
drug pair. This script uses the included reviewed table. It does not repeat the
manual database review.

### 6. BindingDB May 2026

Obtain `BindingDB_All_202605_tsv.zip` from the BindingDB download page or its
long-term archive and place it in `data/raw/bindingdb/`. BindingDB publishes
regular updates, so a differently dated all-data file is not an exact
substitute for the May 2026 release.

BindingDB and ChEMBL records were converted to a common set of drug names,
target names, assay types, and negative-log concentration values in
`data/inputs/binding_assays_raw_harmonized.tsv`. Fit the two-dimensional
PFA representation with:

```bash
python scripts/00_fit_pfa_binding.py
```

## Outputs created by the included code

Use the included prepared inputs to run the QUBO, retrieval, ablation,
supervised comparison, binding sensitivity, and simulated annealing analyses.
Use the preparation scripts to create the consensus outcomes, pathway modules,
mechanism pair scores, and PFA binding representation from the input tables
listed above.

The remaining manual steps are the initial NSC identifier review used to make
the CellMiner GI50 matrix and the database review used to make the drug-level
mechanism annotation table. Both reviewed tables are included. Use them to run
the analyses without repeating those manual steps.

