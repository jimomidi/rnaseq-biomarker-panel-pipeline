# ============================================================
# Step 2 — Download and prepare RNA-seq expression data (cohort 2)
# ============================================================
#
# Cohort 2 does not have enough matched normal samples in its own GDC
# project, so an external control cohort (retrieved via recount3) is
# used as the normal/control reference instead. Because the two
# expression sets come from different sequencing studies, a batch
# effect is expected. This script merges the two sources and TMM-
# normalizes them, and runs a diagnostic (PCA + PERMANOVA) confirming
# that a batch effect is present. The correction itself
# (`limma::removeBatchEffect()`, with the tumor/normal group protected
# as a covariate) is applied consistently by every downstream script
# that compares tumor vs. normal for cohort 2 — see
# `load_cohort2_expression()` in R/utils_expression.R, used by
# R/04 (differential expression), R/06 (panel selection), R/07 (ROC),
# and R/08 (boxplots) — rather than baked into a single "corrected"
# matrix here. This keeps one batch-correction implementation, applied
# identically everywhere the tumor/normal batch confound matters,
# instead of producing a separately corrected matrix that only some
# downstream steps would use.
#
# Steps:
#   1. Download tumor RNA-seq data for cohort 2, restricted to a
#      specific histological subtype (GDC).
#   2. Download the external control expression set (recount3).
#   3. Merge the two count matrices on shared genes, TMM-normalize.
#   4. Run PCA and a PERMANOVA-based diagnostic to quantify how much
#      of the variance is explained by batch (study of origin) vs.
#      biological group (tumor/normal).
#
# Outputs (written under COHORT2_RAW_DIR / PROCESSED_DIR):
#   cohort2_tumor_<assay>.csv                  - GDC tumor expression assays
#   cohort2_tumor_clinical.csv                 - clinical annotation for selected cases
#   control_reference_counts.csv / _tpm.csv    - external control expression (recount3)
#   merged_metadata.csv                        - sample-level batch/group labels
#   merged_logcpm_before_correction.csv        - TMM log-CPM (batch not yet removed)
#   batch_pca_before_correction.png            - PCA colored by batch and group

source("config.R")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
for (pkg in c("TCGAbiolinks", "SummarizedExperiment", "recount3",
              "edgeR", "dplyr", "ggplot2", "vegan")) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE)
}

library(TCGAbiolinks)
library(SummarizedExperiment)
library(recount3)
library(edgeR)
library(dplyr)
library(ggplot2)

set.seed(RANDOM_SEED)
ensure_dir(COHORT2_RAW_DIR)
ensure_dir(PROCESSED_DIR)

#' Flatten any list-valued columns to semicolon-separated strings.
flatten_list_columns <- function(df) {
  as.data.frame(lapply(df, function(col) {
    if (is.list(col)) sapply(col, function(x) paste(unlist(x), collapse = ";")) else col
  }))
}

# ---- 1. Tumor expression (GDC, histological subtype filter) --------------

message(">>> Retrieving clinical annotation for cohort 2 (", GDC_PROJECT_COHORT2, ") ...")
clinical <- GDCquery_clinic(project = GDC_PROJECT_COHORT2, type = "clinical")

hist_col <- intersect(
  c("primary_diagnosis", "histological_type", "disease_type"),
  colnames(clinical)
)[1]
if (is.na(hist_col)) {
  stop("Could not find a histology column in the clinical data. ",
       "Inspect colnames(clinical) and update the candidate list.")
}
message(">>> Using clinical column for histology selection: ", hist_col)

case_ids <- clinical %>%
  filter(grepl(COHORT2_HISTOLOGY_KEYWORD, .data[[hist_col]], ignore.case = TRUE)) %>%
  pull(submitter_id)
message(">>> Cases matching the configured histology keyword: ", length(case_ids))

if (length(case_ids) == 0) {
  stop("No cases matched COHORT2_HISTOLOGY_KEYWORD in config.R; inspect table(clinical[[hist_col]]).")
}

message(">>> Querying tumor RNA-seq data ...")
query_tumor <- GDCquery(
  project       = GDC_PROJECT_COHORT2,
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type   = "Primary Tumor",
  barcode       = case_ids
)
GDCdownload(query_tumor, method = "api", files.per.chunk = 10)
se_tumor <- GDCprepare(query_tumor)
message(sprintf("  Tumor: %d samples, %d genes", ncol(se_tumor), nrow(se_tumor)))

for (assay_name in assayNames(se_tumor)) {
  mat <- assay(se_tumor, assay_name)
  rownames(mat) <- rowData(se_tumor)$gene_name
  out_path <- file.path(COHORT2_RAW_DIR, paste0("cohort2_tumor_", assay_name, ".csv"))
  write.csv(mat, file = out_path, quote = FALSE)
  message(sprintf("  [saved] %s (%d genes x %d samples)", basename(out_path), nrow(mat), ncol(mat)))
}

case_clinical <- clinical %>%
  filter(submitter_id %in% case_ids) %>%
  flatten_list_columns()
write.csv(
  case_clinical,
  file = file.path(COHORT2_RAW_DIR, "cohort2_tumor_clinical.csv"),
  quote = FALSE, row.names = FALSE
)
saveRDS(se_tumor, file = file.path(COHORT2_RAW_DIR, "cohort2_tumor_SE.rds"))

tumor_counts <- assay(se_tumor, "unstranded")
rownames(tumor_counts) <- rowData(se_tumor)$gene_name

# ---- 2. External control expression (recount3) ---------------------------

message(">>> Downloading external control expression data (recount3) ...")
rse_control <- create_rse_manual(
  project      = RECOUNT3_CONTROL_PROJECT,
  project_home = RECOUNT3_CONTROL_PROJECT_HOME,
  organism     = "human",
  annotation   = RECOUNT3_ANNOTATION,
  type         = "gene"
)
message(">>> Control samples retrieved: ", ncol(rse_control))

control_counts <- assay(rse_control, "raw_counts")
rownames(control_counts) <- rowData(rse_control)$gene_name

write.csv(
  control_counts,
  file = file.path(COHORT2_RAW_DIR, "control_reference_counts.csv"),
  quote = FALSE
)

gene_lengths_kb <- rowData(rse_control)$bp_length / 1000
rpk <- control_counts / gene_lengths_kb
control_tpm <- t(t(rpk) / (colSums(rpk) / 1e6))
rownames(control_tpm) <- rowData(rse_control)$gene_name

write.csv(
  control_tpm,
  file = file.path(COHORT2_RAW_DIR, "control_reference_tpm.csv"),
  quote = FALSE
)
message("  [saved] control_reference_counts.csv / control_reference_tpm.csv")

# ---- 3. Merge cohorts and TMM-normalize -----------------------------------

message("\n>>> Merging tumor and control expression matrices ...")

common_genes <- intersect(rownames(tumor_counts), rownames(control_counts))
message(">>> Shared genes: ", length(common_genes))

tumor_common   <- tumor_counts[common_genes, ]
control_common <- control_counts[common_genes, ]
merged_counts  <- cbind(tumor_common, control_common)

batch <- c(rep("cohort2_gdc", ncol(tumor_common)), rep("external_control", ncol(control_common)))
group <- c(rep("Tumor", ncol(tumor_common)), rep("Normal", ncol(control_common)))

metadata <- data.frame(sample = colnames(merged_counts), batch = batch, group = group)
write.csv(
  metadata,
  file = file.path(PROCESSED_DIR, "merged_metadata.csv"),
  quote = FALSE, row.names = FALSE
)

message(">>> TMM normalization ...")
dge <- calcNormFactors(DGEList(counts = merged_counts, group = group), method = "TMM")
logcpm_before <- cpm(dge, log = TRUE, prior.count = 1)
write.csv(
  logcpm_before,
  file = file.path(PROCESSED_DIR, "merged_logcpm_before_correction.csv"),
  quote = FALSE
)

message(sprintf(
  ">>> Merged matrix: %d genes x %d samples (%d tumor, %d normal)",
  nrow(logcpm_before), ncol(logcpm_before), ncol(tumor_common), ncol(control_common)
))

# ---- 4. Batch-effect diagnostic --------------------------------------------
#
# This only demonstrates that a batch effect is present (batch and
# biological group are perfectly confounded here: every tumor sample
# is cohort2_gdc, every normal sample is external_control). It does
# not correct the data — the correction is applied downstream, per
# analysis, via `load_cohort2_expression()`.

#' PCA scatter plot colored by biological group and shaped by batch.
plot_batch_pca <- function(mat, meta, title) {
  mat <- mat[apply(mat, 1, var) > 0, ]
  pca <- prcomp(t(mat), scale. = TRUE)
  var_explained <- round(summary(pca)$importance[2, 1:2] * 100, 1)

  df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Batch = meta$batch, Group = meta$group)

  ggplot(df, aes(x = PC1, y = PC2, color = Group, shape = Batch)) +
    geom_point(size = 3, alpha = 0.8) +
    labs(
      title = title,
      x = paste0("PC1 (", var_explained[1], "%)"),
      y = paste0("PC2 (", var_explained[2], "%)")
    ) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

p_before <- plot_batch_pca(logcpm_before, metadata, "PCA - Before Batch Correction")
ggsave(file.path(PROCESSED_DIR, "batch_pca_before_correction.png"), p_before, width = 7, height = 5, dpi = 150)
message("  [saved] batch_pca_before_correction.png")

#' Quantify the batch effect via PC correlation and PERMANOVA.
#'
#' Reports the correlation of PC1/PC2 with batch and group, plus a
#' PERMANOVA (via vegan::adonis2) partitioning variance between them.
check_batch_effect <- function(mat, meta, label) {
  if (!requireNamespace("vegan", quietly = TRUE)) {
    message("Package 'vegan' not installed; skipping PERMANOVA diagnostic for ", label)
    return(invisible(NULL))
  }
  library(vegan)

  mat_t <- t(mat[apply(mat, 1, var) > 0, ])
  pca <- prcomp(mat_t, scale. = TRUE)

  batch_num <- ifelse(meta$batch == "cohort2_gdc", 1, 0)
  group_num <- ifelse(meta$group == "Tumor", 1, 0)

  cat("\n======", label, "======\n")
  cat("cor(PC1, Batch):", round(cor(pca$x[, 1], batch_num), 3), "\n")
  cat("cor(PC2, Batch):", round(cor(pca$x[, 2], batch_num), 3), "\n")
  cat("cor(PC1, Group):", round(cor(pca$x[, 1], group_num), 3), "\n")
  cat("cor(PC2, Group):", round(cor(pca$x[, 2], group_num), 3), "\n")

  dist_mat <- dist(mat_t)
  perm <- adonis2(dist_mat ~ batch + group, data = meta, permutations = 99)
  cat("\nPERMANOVA:\n")
  print(perm)
}

check_batch_effect(logcpm_before, metadata, "Batch effect diagnostic (uncorrected data)")

message("\n========== Cohort 2 preparation summary ==========")
message("Tumor samples          : ", ncol(tumor_common))
message("Control samples        : ", ncol(control_common))
message("Genes                  : ", nrow(logcpm_before))
message("Processed data written : ", PROCESSED_DIR)
message("Note: batch/group are confounded by design (every tumor sample")
message("is 'cohort2_gdc', every normal sample is 'external_control').")
message("Batch correction is applied downstream at the point of analysis;")
message("see load_cohort2_expression() in R/utils_expression.R.")
message("====================================================")
