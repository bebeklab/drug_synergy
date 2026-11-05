#!/usr/bin/env Rscript

cran_packages <- c("data.table", "dplyr", "tidyr", "igraph", "matrixStats")
bioconductor_packages <- c("limma", "fgsea")

missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_cran) > 0) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
missing_bioc <- bioconductor_packages[
  !vapply(bioconductor_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

