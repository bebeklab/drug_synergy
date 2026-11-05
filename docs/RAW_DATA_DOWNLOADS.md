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

| Source | Release or file used | Official access |
|---|---|---|
| DrugComb / NCI-ALMANAC | `summary_df_ALMANAC.csv`, the project copy of the ALMANAC summary subset | [Project dataset](https://doi.org/10.5281/zenodo.23986099), [direct CSV download](https://zenodo.org/api/records/23986100/files/summary_df_ALMANAC.csv/content), and [original DrugComb v1.4 record](https://zenodo.org/records/11102665) |
| CellMiner / NCI-60 | Compound activity data, including `AVGDATA` where available | [DTP NCI-60 raw compound-activity download](https://discover.nci.nih.gov/cellminer/download/rawdataset/DTP_NCI60_RAW.zip) and [CellMiner dataset index](https://discover.nci.nih.gov/cellminer/loadDownload.do) |
| GEO | `GSE116450_series_matrix.txt.gz` | [Direct series-matrix download](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/GSE116450/matrix/GSE116450_series_matrix.txt.gz) and [GSE116450 record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116450) |
| Ensembl | Release 115 Affymetrix HG-U133 Plus 2 to HGNC mapping | [Ensembl 115 BioMart service](https://sep2025.archive.ensembl.org/biomart/martservice) |
| MSigDB | `c2.cp.kegg.v2023.1.Hs.symbols.gmt` | [Exact MSigDB 2023.1.Hs file (login required)](https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2023.1.Hs/c2.cp.kegg.v2023.1.Hs.symbols.gmt) |
| ChEMBL | ChEMBL 35 SQLite release | [Direct ChEMBL 35 SQLite archive](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_35/chembl_35_sqlite.tar.gz) and [release directory](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_35/) |
| BindingDB | May 2026 all-data TSV, `BindingDB_All_202605_tsv.zip` | The exact monthly file is no longer on the live download page; use the included harmonized table for this analysis. See [BindingDB downloads](https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp) and the [UCSD quarterly archive](https://library.ucsd.edu/dc/collection/bb03870458). |
| NCI Drug Dictionary | Individual drug records used to review incomplete mechanism annotations | [NCI Drug Dictionary](https://www.cancer.gov/publications/dictionaries/cancer-drug) |
| PubChem | Individual compound records used to review names, structures, and mechanisms | [PubChem](https://pubchem.ncbi.nlm.nih.gov/) |

## Download steps

### 1. DrugComb outcomes

The analysis started from [`summary_df_ALMANAC.csv`](https://zenodo.org/api/records/23986100/files/summary_df_ALMANAC.csv/content).
It contains the NCI-ALMANAC records selected from DrugComb, including Bliss,
Loewe, ZIP, and HSA scores. The file is about 291 MB and is archived in the
[project's Zenodo dataset](https://doi.org/10.5281/zenodo.23986099) rather than
in GitHub. After downloading it, save it as
`data/raw/drugcomb/summary_df_ALMANAC.csv`.

Run `python scripts/00_build_consensus_outcomes.py` to create one observed
consensus-synergy score for each cancer and drug pair. The public DrugComb v1.4
archive contains a broader summary table. Using that broader file requires
selecting the same ALMANAC records before running the script.

### 2. CellMiner GI50 data

Download the [DTP NCI-60 raw compound-activity archive](https://discover.nci.nih.gov/cellminer/download/rawdataset/DTP_NCI60_RAW.zip)
and save it under `data/raw/cellminer/`. The broader
[CellMiner dataset index](https://discover.nci.nih.gov/cellminer/loadDownload.do)
also lists the current raw and processed NCI-60 files. Match the drug records with the NSC
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

MSigDB requires a free account. Sign in and download the
[exact 2023.1.Hs KEGG file](https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2023.1.Hs/c2.cp.kegg.v2023.1.Hs.symbols.gmt),
`c2.cp.kegg.v2023.1.Hs.symbols.gmt`. The link redirects to the MSigDB login
page when no active session is present. Do not replace this file with the
current KEGG collection when reproducing the reported pathway input.

The pathway pair costs read by the analysis are stored in
`data/inputs/pathway_module_pair_costs.csv`.

Create the module assignments and pair costs from the ALMANAC response, GEO
expression, fixed Ensembl probe map, and GMT download with:

```bash
Rscript scripts/install_r_dependencies.R
Rscript scripts/build_pathway_modules.R
```

The script uses seed 42 and one worker by default so new runs are reproducible.
The seed used for the earlier pathway build was not retained, so rebuilding
the pathway input may produce slightly different results. The included
`data/inputs/pathway_module_pair_costs.csv` is the pathway input used for the
reported analysis.

The pathway workflow uses ALMANAC single-drug relative-inhibition values to
define sensitive and resistant cell lines. CellMiner GI50 is a separate QUBO
input and is not substituted for relative inhibition in this pathway step.

### 5. ChEMBL 35

Download the [ChEMBL 35 SQLite archive](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_35/chembl_35_sqlite.tar.gz)
from the [ChEMBL 35 release directory](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_35/),
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

The binding input was assembled from the May 2026 monthly file
`BindingDB_All_202605_tsv.zip`. BindingDB's live download page now serves newer
monthly files, and its UCSD archive contains quarterly snapshots rather than
this exact May file. Therefore, neither the current live download nor the
nearest quarterly snapshot is an exact replacement. The harmonized records
used by this repository are included in
`data/inputs/binding_assays_raw_harmonized.tsv`. Use the
[BindingDB download page](https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp)
for current data or the [UCSD BindingDB archive](https://library.ucsd.edu/dc/collection/bb03870458)
for stable quarterly releases.

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

