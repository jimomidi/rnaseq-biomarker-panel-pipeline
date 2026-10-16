# ============================================================
# Step 1 — Download and prepare RNA-seq expression data (cohort 1)
# ============================================================
#
# Retrieves gene-level RNA-seq quantification for cohort 1 from the
# GDC (Genomic Data Commons) via TCGAbiolinks:
#   - Tumor samples restricted to a specific subtype, identified by
#     keyword matching against the clinical annotation.
#   - All available solid-tissue-normal samples from the same project
#     (unfiltered by subtype).
#
# Outputs (written under COHORT1_RAW_DIR):
#   cohort1_tumor_<assay>.csv   - one CSV per expression assay (tumor)
#   cohort1_normal_<assay>.csv  - one CSV per expression assay (normal)
#   cohort1_tumor_clinical.csv  - clinical annotation for the selected tumor cases
#   cohort1_tumor_SE.rds        - full SummarizedExperiment object (tumor)
#   cohort1_normal_SE.rds       - full SummarizedExperiment object (normal)

source("config.R")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
for (pkg in c("TCGAbiolinks", "SummarizedExperiment", "dplyr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE)
}

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)

ensure_dir(COHORT1_RAW_DIR)

#' Identify the clinical column that best describes tumor site/diagnosis.
#'
#' @param clinical Clinical data frame returned by `GDCquery_clinic()`.
#' @param candidates Candidate column names, in order of preference.
#' @return The name of the first matching column, or `NA` if none exist.
select_site_column <- function(clinical, candidates) {
  match <- intersect(candidates, colnames(clinical))
  if (length(match) == 0) return(NA_character_)
  match[1]
}

#' Filter clinical records to cases whose site/diagnosis column matches
#' any of the given keywords (case-insensitive).
#'
#' @param clinical Clinical data frame.
#' @param site_col Column to search.
#' @param keywords Character vector of keywords (combined with OR).
#' @return Character vector of matching `submitter_id` values.
filter_case_ids_by_keyword <- function(clinical, site_col, keywords) {
  pattern <- paste(keywords, collapse = "|")
  clinical %>%
    filter(grepl(pattern, .data[[site_col]], ignore.case = TRUE)) %>%
    pull(submitter_id)
}

#' Flatten any list-valued columns to semicolon-separated strings so the
#' data frame can be written safely with `write.csv()`.
flatten_list_columns <- function(df) {
  as.data.frame(lapply(df, function(col) {
    if (is.list(col)) sapply(col, function(x) paste(unlist(x), collapse = ";")) else col
  }))
}

#' Download, prepare, and save all expression assays for a GDC query.
#'
#' @param query A `GDCquery()` object.
#' @param out_dir Output directory.
#' @param file_prefix Prefix used for the saved CSV/RDS files.
#' @return The prepared `SummarizedExperiment` object.
download_and_save_expression <- function(query, out_dir, file_prefix) {
  GDCdownload(query, method = "api", files.per.chunk = 10)
  se <- GDCprepare(query)

  message(sprintf("  %s: %d samples, %d genes", file_prefix, ncol(se), nrow(se)))

  for (assay_name in assayNames(se)) {
    mat <- assay(se, assay_name)
    rownames(mat) <- rowData(se)$gene_name
    out_path <- file.path(out_dir, paste0(file_prefix, "_", assay_name, ".csv"))
    write.csv(mat, file = out_path, quote = FALSE)
    message(sprintf("  [saved] %s (%d genes x %d samples)", basename(out_path), nrow(mat), ncol(mat)))
  }

  saveRDS(se, file = file.path(out_dir, paste0(file_prefix, "_SE.rds")))
  se
}

# ---- 1. Clinical annotation and tumor-subtype selection -----------------

message(">>> Retrieving clinical annotation for cohort 1 (", GDC_PROJECT_COHORT1, ") ...")
clinical <- GDCquery_clinic(project = GDC_PROJECT_COHORT1, type = "clinical")

site_col <- select_site_column(
  clinical,
  c("tissue_or_organ_of_origin", "primary_site", "primary_diagnosis")
)
if (is.na(site_col)) {
  stop("Could not find a tumor-site/diagnosis column in the clinical data. ",
       "Inspect colnames(clinical) and update select_site_column().")
}
message(">>> Using clinical column for subtype selection: ", site_col)

case_ids <- filter_case_ids_by_keyword(clinical, site_col, COHORT1_SUBTYPE_KEYWORDS)
message(">>> Cases matching the configured subtype keywords: ", length(case_ids))

if (length(case_ids) == 0) {
  message(">>> No cases matched. Available values in '", site_col, "':")
  print(table(clinical[[site_col]]))
  stop("No cases matched COHORT1_SUBTYPE_KEYWORDS in config.R; update the keyword list.")
}

# ---- 2. Tumor expression -------------------------------------------------

message(">>> Querying tumor RNA-seq data ...")
query_tumor <- GDCquery(
  project       = GDC_PROJECT_COHORT1,
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type   = "Primary Tumor",
  barcode       = case_ids
)
se_tumor <- download_and_save_expression(query_tumor, COHORT1_RAW_DIR, "cohort1_tumor")

case_clinical <- clinical %>%
  filter(submitter_id %in% case_ids) %>%
  flatten_list_columns()

write.csv(
  case_clinical,
  file = file.path(COHORT1_RAW_DIR, "cohort1_tumor_clinical.csv"),
  quote = FALSE, row.names = FALSE
)
message("  [saved] cohort1_tumor_clinical.csv")

# ---- 3. Normal expression (all normal samples in the project) -----------

message(">>> Querying normal-tissue RNA-seq data (all samples in the project) ...")
query_normal <- GDCquery(
  project       = GDC_PROJECT_COHORT1,
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type   = "Solid Tissue Normal"
)
se_normal <- download_and_save_expression(query_normal, COHORT1_RAW_DIR, "cohort1_normal")

# ---- 4. Summary -----------------------------------------------------------

message("\n========== Cohort 1 download summary ==========")
message("GDC project        : ", GDC_PROJECT_COHORT1)
message("Tumor samples       : ", ncol(se_tumor))
message("Normal samples      : ", ncol(se_normal))
message("Genes                : ", nrow(se_tumor))
message("Output directory     : ", COHORT1_RAW_DIR)
message("=================================================")
