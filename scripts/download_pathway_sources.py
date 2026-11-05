#!/usr/bin/env python3
"""Download public pathway inputs that do not require an account."""

from __future__ import annotations

import argparse
import urllib.parse
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"

FILES = {
    RAW / "drugcomb" / "summary_table_v1.4.csv": (
        "https://zenodo.org/records/11102665/files/summary_table_v1.4.csv?download=1"
    ),
    RAW / "geo" / "GSE116450_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/GSE116450/"
        "matrix/GSE116450_series_matrix.txt.gz"
    ),
}

ENSEMBL_URL = "https://sep2025.archive.ensembl.org/biomart/martservice"
ENSEMBL_QUERY = (
    '<Query virtualSchemaName="default" formatter="TSV" header="1" '
    'uniqueRows="1" count="" datasetConfigVersion="0.6">'
    '<Dataset name="hsapiens_gene_ensembl" interface="default">'
    '<Attribute name="hgnc_symbol"/>'
    '<Attribute name="affy_hg_u133_plus_2"/>'
    "</Dataset></Query>"
)


def download(url: str, destination: Path, force: bool) -> None:
    if destination.exists() and not force:
        print(f"Exists: {destination}")
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".part")
    print(f"Downloading {url}")
    urllib.request.urlretrieve(url, temporary)
    temporary.replace(destination)
    print(f"Wrote: {destination}")


def download_ensembl_mapping(destination: Path, force: bool) -> None:
    if destination.exists() and not force:
        print(f"Exists: {destination}")
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = urllib.parse.urlencode({"query": ENSEMBL_QUERY}).encode("utf-8")
    request = urllib.request.Request(ENSEMBL_URL, data=payload, method="POST")
    print("Downloading fixed Ensembl 115 probe mapping")
    with urllib.request.urlopen(request, timeout=300) as response:
        content = response.read()
    if not content.startswith(b"HGNC symbol\tAFFY HG U133 Plus 2 probe"):
        raise RuntimeError("Unexpected response from the Ensembl BioMart archive")
    destination.write_bytes(content)
    print(f"Wrote: {destination}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true", help="replace existing files")
    args = parser.parse_args()
    for destination, url in FILES.items():
        download(url, destination, args.force)
    download_ensembl_mapping(
        RAW / "ensembl" / "ensembl115_affy_hgnc.tsv", args.force
    )
    print(
        "MSigDB requires login. Download "
        "c2.cp.kegg.v2023.1.Hs.symbols.gmt manually into data/raw/msigdb/."
    )


if __name__ == "__main__":
    main()
