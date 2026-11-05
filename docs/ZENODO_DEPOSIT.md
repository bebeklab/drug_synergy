# Zenodo deposit

Use separate Zenodo records for the software and the large data file.

## Software record

1. Push this repository to GitHub.
2. In Zenodo, open the GitHub integration page and enable the repository.
3. Create a GitHub release, such as `v1.0.0`.
4. Zenodo will archive that release and assign a software DOI.

Do not add analysis results or large source data to the GitHub release.

## Data record

Create a Zenodo upload with resource type **Dataset** and upload these files:

| File | Purpose |
|---|---|
| `summary_df_ALMANAC.csv` | Exact NCI-ALMANAC/DrugComb summary table read by `scripts/00_build_consensus_outcomes.py` |
| `README.txt` | Describes the file, its source, the expected local path, and the command that reads it |
| `SHA256SUMS.txt` | Records the checksum used to confirm file integrity |

The checksum for the project copy of `summary_df_ALMANAC.csv` is:

```text
9fdb360a98d38ff40943bb709197be5cd862ce7582fbfb5d46863e9f7bc999ae  summary_df_ALMANAC.csv
```

After downloading the dataset, place the CSV at:

```text
data/raw/drugcomb/summary_df_ALMANAC.csv
```

Then run:

```bash
python scripts/00_build_consensus_outcomes.py
```

Do not upload the following to this data record:

- `results/` outputs;
- files already stored in `data/inputs/` or `data/annotations/` on GitHub;
- GSE116450, ChEMBL, BindingDB, or MSigDB downloads;
- Python environments, caches, or temporary files.

Use the official source links in `docs/RAW_DATA_DOWNLOADS.md` for the external
databases. Check the source terms before redistributing any third-party file.

## Suggested metadata

- **Title:** Source data for Quantum-Inspired Optimization of Cancer-Specific Drug Combinations Without Combination-Response Training Data
- **Resource type:** Dataset
- **Creators:** Use the manuscript author order and ORCID identifiers
- **Description:** State that the record supplies the exact NCI-ALMANAC summary table used to create the cancer-specific consensus-synergy outcomes
- **Related identifiers:** Add the software DOI, GitHub repository, and manuscript DOI when available
- **Version:** `1.0.0`

Reserve the data DOI before publishing if it must appear in the repository or
manuscript. Add the final DOI to the main README and the manuscript's Data and
Code Availability section.
