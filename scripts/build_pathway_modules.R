#!/usr/bin/env Rscript

# Build cancer-specific pathway-module costs from the source files. The fixed
# seed makes new runs repeatable. Permutation results near the significance
# threshold can differ from an earlier build whose random-number state was not
# retained. Only pair costs and module assignments are written to disk;
# gene-level and enrichment tables stay in memory.

required_packages <- c(
  "data.table", "dplyr", "tidyr", "limma", "fgsea", "igraph",
  "matrixStats"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    ". Run: Rscript scripts/install_r_dependencies.R"
  )
}

args <- commandArgs(trailingOnly = TRUE)
argument <- function(name, default) {
  prefix <- paste0("--", name, "=")
  match <- args[startsWith(args, prefix)]
  if (length(match) == 0) return(default)
  sub(prefix, "", match[[1]], fixed = TRUE)
}

script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]))
project <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

almanac_path <- argument(
  "almanac",
  file.path(project, "data/raw/drugcomb/summary_df_ALMANAC.csv")
)
expression_path <- argument(
  "expression",
  file.path(project, "data/raw/geo/GSE116450_series_matrix.txt.gz")
)
kegg_path <- argument(
  "kegg",
  file.path(project, "data/raw/msigdb/c2.cp.kegg.v2023.1.Hs.symbols.gmt")
)
probe_map_path <- argument(
  "probe-map",
  file.path(project, "data/raw/ensembl/ensembl115_affy_hgnc.tsv")
)
output_dir <- argument(
  "output-dir", file.path(project, "results/generated/pathway_rebuild")
)
seed <- as.integer(argument("seed", "42"))
nperm <- as.integer(argument("permutations", "10000"))
adjusted_p_cutoff <- as.numeric(argument("adjusted-p-cutoff", "0.05"))
workers <- as.integer(argument("workers", "1"))

for (path in c(almanac_path, expression_path, kegg_path, probe_map_path)) {
  if (!file.exists(path)) stop("Required input is missing: ", path)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

canonical <- function(value) {
  trimws(tolower(as.character(value)))
}

cancer_names <- c(
  brain = "Brain",
  breast = "Breast",
  haematopoietic_and_lymphoid = "Leukemia",
  kidney = "Kidney",
  large_intestine = "Colon",
  lung = "Lung",
  ovary = "Ovarian",
  skin = "Melanoma"
)

cell_name_map <- c(
  "COLO 205" = "COLO-205", "HCT116" = "HCT-116",
  "HL-60(TB)" = "HL-60", "HS 578T" = "HS-578T",
  "IGROV1" = "IGR-OV1", "LOX IMVI" = "LOX",
  "NCI/ADR-RES" = "NCI-ADR-RES", "NCIH23" = "NCI-H23",
  "OVCAR3" = "OVCAR-3", "RXF 393" = "RXF-393",
  "UACC62" = "UACC-62"
)

message("Reading ALMANAC single-drug responses")
required_almanac <- c(
  "drug_row", "drug_col", "cell_line_name", "tissue_name", "ri_row", "ri_col"
)
almanac <- data.table::fread(
  almanac_path,
  select = required_almanac,
  na.strings = c("", "NA", "\\N"),
  showProgress = FALSE,
  data.table = FALSE
)
missing_columns <- setdiff(required_almanac, names(almanac))
if (length(missing_columns) > 0) {
  stop("ALMANAC file is missing columns: ", paste(missing_columns, collapse = ", "))
}

response <- dplyr::bind_rows(
  dplyr::transmute(
    almanac, tissue_name, cell_line_name, drug = drug_row,
    relative_inhibition = as.numeric(ri_row)
  ),
  dplyr::transmute(
    almanac, tissue_name, cell_line_name, drug = drug_col,
    relative_inhibition = as.numeric(ri_col)
  )
) |>
  dplyr::mutate(
    cancer = unname(cancer_names[gsub(" ", "_", canonical(tissue_name))]),
    cell_line_name = dplyr::recode(cell_line_name, !!!cell_name_map)
  ) |>
  dplyr::filter(!is.na(cancer), !is.na(relative_inhibition)) |>
  dplyr::group_by(cancer, drug, cell_line_name) |>
  dplyr::summarise(
    relative_inhibition = mean(relative_inhibition, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::group_by(cancer, drug) |>
  dplyr::arrange(dplyr::desc(relative_inhibition), .by_group = TRUE) |>
  dplyr::mutate(
    response_rank = dplyr::row_number(),
    n_cell_lines = dplyr::n(),
    response_group = dplyr::case_when(
      n_cell_lines == 1 ~ NA_character_,
      n_cell_lines == 2 & response_rank == 1 ~ "Sensitive",
      n_cell_lines == 2 & response_rank == 2 ~ "Resistant",
      dplyr::ntile(dplyr::desc(relative_inhibition), 3) == 1 ~ "Sensitive",
      dplyr::ntile(dplyr::desc(relative_inhibition), 3) == 3 ~ "Resistant",
      TRUE ~ "Intermediate"
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::filter(response_group %in% c("Sensitive", "Resistant"))

message("Reading untreated 0-nM, 24-hour GSE116450 expression")
expression <- data.table::fread(
  expression_path,
  sep = "\t",
  skip = "ID_REF",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE,
  showProgress = FALSE
)
expression <- expression[expression[[1]] != "!series_matrix_table_end", , drop = FALSE]
# The matrix columns are GEO sample accessions. Replace them with the sample
# titles, which encode cell line, topotecan concentration, and exposure time.
metadata <- readLines(expression_path, warn = FALSE)
title_line <- metadata[grepl("^!Sample_title", metadata)]
if (length(title_line) != 1) stop("Could not identify one GEO !Sample_title line")
sample_titles <- strsplit(title_line, "\t", fixed = TRUE)[[1]][-1]
sample_titles <- gsub('^"|"$', "", sample_titles)
if (length(sample_titles) != ncol(expression) - 1) {
  stop("GEO sample-title count does not match the expression matrix")
}
names(expression)[-1] <- sample_titles
sample_names <- names(expression)[-1]
baseline_samples <- grepl("_topotecan_0nM_24h$", sample_names)
if (sum(baseline_samples) != 60) {
  stop("Expected 60 untreated 24-hour samples; found ", sum(baseline_samples))
}
baseline <- expression[, c(TRUE, baseline_samples), drop = FALSE]
cell_lines <- sub("_topotecan_0nM_24h$", "", names(baseline)[-1])
probe_ids <- as.character(baseline[[1]])
expression_matrix <- t(as.matrix(baseline[-1]))
storage.mode(expression_matrix) <- "double"
rownames(expression_matrix) <- cell_lines
colnames(expression_matrix) <- probe_ids

message("Mapping probes with the fixed Ensembl 115 BioMart export")
probe_map <- data.table::fread(
  probe_map_path,
  sep = "\t",
  header = TRUE,
  select = c("HGNC symbol", "AFFY HG U133 Plus 2 probe"),
  na.strings = "",
  data.table = FALSE,
  showProgress = FALSE
)
names(probe_map) <- c("gene", "probe")
probe_map <- probe_map |>
  dplyr::filter(
    !is.na(gene), !is.na(probe),
    gene != "", probe %in% colnames(expression_matrix)
  ) |>
  dplyr::distinct(gene, probe)

# Preserve one probe per gene by selecting the probe with the largest IQR
# across the 60 untreated samples, matching the manuscript workflow.
probe_iqr <- matrixStats::colIQRs(expression_matrix)
iqr_table <- data.frame(
  probe = colnames(expression_matrix), iqr = probe_iqr,
  stringsAsFactors = FALSE
)
best_probes <- probe_map |>
  dplyr::left_join(iqr_table, by = "probe") |>
  dplyr::group_by(gene) |>
  dplyr::slice_max(iqr, n = 1, with_ties = FALSE) |>
  dplyr::ungroup()
gene_expression <- t(expression_matrix[, best_probes$probe, drop = FALSE])
rownames(gene_expression) <- best_probes$gene

missing_expression <- setdiff(unique(response$cell_line_name), colnames(gene_expression))
if (length(missing_expression) > 0) {
  stop(
    "Response cell lines missing from GSE116450: ",
    paste(missing_expression, collapse = ", ")
  )
}

message("Calculating sensitive-versus-resistant expression contrasts")
# Preserve the ordering used by the original analysis. This matters for
# permutation-based enrichment because each group consumes the next segment of
# R's random-number stream.
response_groups <- split(
  response,
  list(response$cancer, response$drug),
  drop = TRUE
)
gene_rankings <- lapply(response_groups, function(group) {
  group <- group[!duplicated(group$cell_line_name), , drop = FALSE]
  if (length(unique(group$response_group)) < 2) return(NULL)
  design_group <- factor(
    group$response_group,
    levels = c("Resistant", "Sensitive")
  )
  design <- stats::model.matrix(~ 0 + design_group)
  colnames(design) <- c("Resistant", "Sensitive")
  fit <- limma::lmFit(
    gene_expression[, group$cell_line_name, drop = FALSE], design
  )
  fit <- limma::contrasts.fit(
    fit,
    limma::makeContrasts(Sensitive - Resistant, levels = design)
  )
  statistic <- limma::eBayes(fit)$t[, 1]
  sort(statistic, decreasing = TRUE)
})
gene_rankings <- gene_rankings[!vapply(gene_rankings, is.null, logical(1))]

kegg <- fgsea::gmtPathways(kegg_path)
message(
  "Running KEGG enrichment for ", length(gene_rankings),
  " cancer-drug groups with ", nperm, " permutations"
)
enrichment_one <- function(key) {
  result <- suppressWarnings(fgsea::fgsea(
    pathways = kegg,
    stats = gene_rankings[[key]],
    minSize = 10,
    maxSize = 500,
    nperm = nperm
  ))
  result$key <- key
  result
}
if (workers > 1 && .Platform$OS.type != "windows") {
  enrichment <- parallel::mclapply(
    names(gene_rankings), enrichment_one,
    mc.cores = workers, mc.set.seed = TRUE
  )
} else {
  enrichment <- lapply(names(gene_rankings), enrichment_one)
}
enrichment <- data.table::rbindlist(enrichment, fill = TRUE)
enrichment <- enrichment[!is.na(padj) & padj < adjusted_p_cutoff]
if (nrow(enrichment) == 0) stop("No pathways passed adjusted P < 0.05")

key_lookup <- response |>
  dplyr::distinct(cancer, drug) |>
  dplyr::mutate(key = interaction(cancer, drug, drop = TRUE)) |>
  dplyr::mutate(key = as.character(key))
enrichment <- dplyr::left_join(
  as.data.frame(enrichment), key_lookup, by = "key"
)

message("Building pathway networks and Walktrap modules")
pair_outputs <- list()
module_outputs <- list()
for (cancer_name in sort(unique(enrichment$cancer))) {
  cancer_enrichment <- enrichment[enrichment$cancer == cancer_name, , drop = FALSE]
  collapsed <- cancer_enrichment |>
    dplyr::group_by(drug, pathway) |>
    dplyr::summarise(NES = mean(NES, na.rm = TRUE), .groups = "drop")
  drugs <- sort(unique(collapsed$drug))
  pathways <- sort(unique(collapsed$pathway))
  nes <- tidyr::pivot_wider(
    tidyr::expand_grid(pathway = pathways, drug = drugs) |>
      dplyr::left_join(collapsed, by = c("pathway", "drug")) |>
      dplyr::mutate(NES = tidyr::replace_na(NES, 0)),
    id_cols = pathway,
    names_from = drug,
    values_from = NES,
    values_fill = 0
  ) |>
    dplyr::arrange(pathway)
  pathway_names <- nes$pathway
  nes <- as.matrix(nes[, -1, drop = FALSE])
  rownames(nes) <- pathway_names
  maximum <- apply(abs(nes), 2, max, na.rm = TRUE)
  maximum[maximum == 0] <- 1
  nes <- sweep(nes, 2, maximum, "/")

  signs <- sign(nes)
  shared_direction <- crossprod(signs > 0) + crossprod(signs < 0)
  diag(shared_direction) <- 0
  drugs <- colnames(shared_direction)
  edge_index <- which(upper.tri(shared_direction) & shared_direction > 0, arr.ind = TRUE)
  edges <- data.frame(
    from = drugs[edge_index[, 1]],
    to = drugs[edge_index[, 2]],
    weight = shared_direction[edge_index],
    stringsAsFactors = FALSE
  )
  graph <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = drugs)
  membership_all <- setNames(rep(NA_integer_, length(drugs)), drugs)
  connected <- igraph::V(graph)$name[igraph::degree(graph) > 0]
  if (length(connected) > 1) {
    connected_graph <- igraph::induced_subgraph(graph, connected)
    clusters <- igraph::cluster_walktrap(
      connected_graph,
      weights = igraph::E(connected_graph)$weight
    )
    membership_all[names(igraph::membership(clusters))] <- igraph::membership(clusters)
  }
  next_module <- if (all(is.na(membership_all))) 1L else max(membership_all, na.rm = TRUE) + 1L
  for (drug in names(membership_all)[is.na(membership_all)]) {
    membership_all[[drug]] <- next_module
    next_module <- next_module + 1L
  }

  strength <- igraph::strength(graph, vids = drugs, weights = igraph::E(graph)$weight)
  confidence <- rank(strength, ties.method = "average") / length(strength)
  confidence[strength == 0] <- 0
  names(confidence) <- drugs

  pair_index <- which(upper.tri(shared_direction), arr.ind = TRUE)
  drug1 <- drugs[pair_index[, 1]]
  drug2 <- drugs[pair_index[, 2]]
  cross_module <- membership_all[drug1] != membership_all[drug2]
  pair_confidence <- pmin(confidence[drug1], confidence[drug2])
  pair_outputs[[cancer_name]] <- data.frame(
    cancer = cancer_name,
    drug1 = drug1,
    drug2 = drug2,
    cross_module = as.integer(cross_module),
    pair_pathway_confidence = as.numeric(pair_confidence),
    structured_pathway_cost = as.numeric(!cross_module) + (1 - pair_confidence),
    stringsAsFactors = FALSE
  )
  module_outputs[[cancer_name]] <- data.frame(
    cancer = cancer_name,
    drug = drugs,
    pathway_module = as.integer(membership_all[drugs]),
    pathway_strength = as.numeric(strength[drugs]),
    pathway_confidence = as.numeric(confidence[drugs]),
    stringsAsFactors = FALSE
  )
}

pair_output <- do.call(rbind, pair_outputs)
module_output <- do.call(rbind, module_outputs)
utils::write.csv(
  pair_output,
  file.path(output_dir, "pathway_module_pair_costs.csv"),
  row.names = FALSE
)
utils::write.csv(
  module_output,
  file.path(output_dir, "pathway_drug_modules.csv"),
  row.names = FALSE
)
message("Wrote ", nrow(pair_output), " pair costs to ", output_dir)

