# ============================================================
# Shared helpers for loading cohort expression matrices
# ============================================================
#
# Differential expression (R/04), gene panel selection, ROC evaluation,
# and boxplot visualization (R/06-08) all compare tumor vs. normal for
# cohort 2 and therefore all need the same batch-corrected, low-
# expression-filtered log-CPM matrix. Centralizing that loading logic
# here (rather than repeating or re-deriving it per script) guarantees
# every analysis that makes this comparison sees identically prepared
# data, and keeps R/04, R/06, R/07, and R/08 free of repeated
# data-preparation code.
#
# Requires config.R to have been sourced already (for
# LOW_EXPRESSION_LOGCPM / LOW_EXPRESSION_FRACTION / PROCESSED_DIR /
# COHORT1_RAW_DIR) and the `edgeR` and `limma` packages to be loaded.

#' Load a single-assay expression CSV where the gene identifier is the
#' first (unnamed) column, deduplicating repeated gene names.
load_expression_matrix <- function(path) {
  m <- read.csv(path, row.names = NULL, check.names = FALSE)
  m <- m[!duplicated(m[, 1]), ]
  rownames(m) <- m[, 1]
  as.matrix(m[, -1])
}

#' Apply the pipeline's standard low-expression filter to a log-CPM matrix.
#'
#' A gene is retained if log-CPM > LOW_EXPRESSION_LOGCPM in at least
#' LOW_EXPRESSION_FRACTION of samples.
filter_low_expression <- function(logcpm_mat) {
  min_samples <- ceiling(LOW_EXPRESSION_FRACTION * ncol(logcpm_mat))
  keep <- rowSums(logcpm_mat > LOW_EXPRESSION_LOGCPM) >= min_samples
  logcpm_mat[keep, ]
}

#' Load cohort 1's low-expression-filtered, TMM log-CPM matrix and
#' tumor/normal labels, built from raw tumor and normal counts.
#'
#' @return List with `mat` (genes x samples) and `label`
#'   (character vector of "Tumor"/"Normal", one per column of `mat`).
load_cohort1_expression <- function() {
  tumor_mat  <- load_expression_matrix(file.path(COHORT1_RAW_DIR, "cohort1_tumor_unstranded.csv"))
  normal_mat <- load_expression_matrix(file.path(COHORT1_RAW_DIR, "cohort1_normal_unstranded.csv"))

  common_genes <- intersect(rownames(tumor_mat), rownames(normal_mat))
  raw_counts <- cbind(tumor_mat[common_genes, ], normal_mat[common_genes, ])
  label <- c(rep("Tumor", ncol(tumor_mat)), rep("Normal", ncol(normal_mat)))

  dge <- calcNormFactors(DGEList(raw_counts), method = "TMM")
  logcpm <- cpm(dge, log = TRUE, prior.count = 1)

  list(mat = filter_low_expression(logcpm), label = label)
}

#' Load cohort 2's low-expression-filtered, batch-corrected log-CPM
#' matrix and tumor/normal labels.
#'
#' Tumor samples and the external control samples come from different
#' studies (see R/02_download_prepare_cohort2.R), so this applies the
#' same `limma::removeBatchEffect()` correction used for differential
#' expression (R/04) here too, with the tumor/normal group protected
#' as a covariate. Every analysis that compares tumor vs. normal for
#' cohort 2 (DEG, gene panel ranking, ROC, boxplots) should use this
#' loader rather than reading the raw merged log-CPM file directly, so
#' that none of them are confounded by the uncorrected batch effect.
#'
#' @return List with `mat` (genes x samples, batch-corrected),
#'   `label` (character vector of "Tumor"/"Normal", one per column of
#'   `mat`), and `meta` (the metadata rows for those columns).
load_cohort2_expression <- function() {
  mat  <- as.matrix(read.csv(file.path(PROCESSED_DIR, "merged_logcpm_before_correction.csv"),
                              row.names = 1, check.names = FALSE))
  meta <- read.csv(file.path(PROCESSED_DIR, "merged_metadata.csv"))
  mat  <- mat[, meta$sample]

  mat_filtered <- filter_low_expression(mat)

  group <- factor(meta$group, levels = c("Normal", "Tumor"))
  batch <- factor(meta$batch)
  mat_corrected <- removeBatchEffect(mat_filtered, batch = batch, design = model.matrix(~group))

  label <- ifelse(meta$group == "Tumor", "Tumor", "Normal")
  list(mat = mat_corrected, label = label, meta = meta)
}

#' Load cohort 1's tumor-only raw count matrix (used by the survival
#' scripts, which score genes directly on raw counts for this cohort).
load_cohort1_tumor_matrix <- function() {
  load_expression_matrix(file.path(COHORT1_RAW_DIR, "cohort1_tumor_unstranded.csv"))
}

#' Load cohort 2's tumor-only log-CPM matrix (used by the survival
#' scripts). Batch correction is not applicable here: all cohort 2
#' tumor samples come from the same GDC batch, so there is no
#' tumor-vs-normal batch confound to remove for a tumor-only analysis.
load_cohort2_tumor_matrix <- function() {
  mat  <- as.matrix(read.csv(file.path(PROCESSED_DIR, "merged_logcpm_before_correction.csv"),
                              row.names = 1, check.names = FALSE))
  meta <- read.csv(file.path(PROCESSED_DIR, "merged_metadata.csv"))
  tumor_meta <- meta[meta$group == "Tumor", ]
  mat <- mat[, tumor_meta$sample, drop = FALSE]
  colnames(mat) <- tumor_meta$sample
  mat
}
