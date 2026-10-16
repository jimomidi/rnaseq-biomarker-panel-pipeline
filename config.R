# ============================================================
# Global configuration for the RNA-seq biomarker panel pipeline
# ============================================================
#
# Every script in R/ starts with `source("config.R")`, so all
# directory paths and analysis thresholds are defined in exactly
# one place. Run scripts from the repository root (e.g. with
# `Rscript R/03_differential_expression_cohort1.R` from the repo
# root, or by opening an RStudio project at the repo root) so
# that these relative paths resolve correctly.
#
# Nothing in this file downloads or modifies data; it only
# defines constants that the numbered scripts consume.

# ---- Reproducibility -----------------------------------------

RANDOM_SEED <- 42

# ---- Public data sources (edit for your own cohorts) ----------
#
# This pipeline is written for two independent cohorts that are
# each compared against a matched normal/control reference. The
# scripts are agnostic to which cohorts you use; set the GDC
# project identifiers and case-defining keywords below to match
# your own study.

# GDC (Genomic Data Commons) project identifiers, e.g. "TCGA-XXXX".
GDC_PROJECT_COHORT1 <- "TCGA-XXXX"
GDC_PROJECT_COHORT2 <- "TCGA-YYYY"

# Keywords used to restrict cohort 1 to a specific tumor subtype
# based on its clinical annotation (case-insensitive substring
# match against the clinical site/diagnosis column). Replace with
# the keywords that define your subtype of interest.
COHORT1_SUBTYPE_KEYWORDS <- c(
  "keyword_1", "keyword_2", "keyword_3"
)

# Keyword used to restrict cohort 2 to a specific histological
# subtype (case-insensitive substring match).
COHORT2_HISTOLOGY_KEYWORD <- "keyword"

# External normal/control reference for cohort 2, retrieved via
# the recount3 package (used because the GDC project for cohort 2
# does not provide enough matched normal samples on its own).
# `RECOUNT3_CONTROL_PROJECT_HOME` must be one of recount3's own
# source identifiers (see the recount3 documentation for the current
# list) — set it to whichever source hosts `RECOUNT3_CONTROL_PROJECT`.
RECOUNT3_CONTROL_PROJECT      <- "TISSUE_NAME"
RECOUNT3_CONTROL_PROJECT_HOME <- "data_sources/SOURCE_NAME"
RECOUNT3_ANNOTATION           <- "gencode_v26"

# ---- Directory layout -------------------------------------------

DATA_DIR      <- "data"
RAW_DIR       <- file.path(DATA_DIR, "raw")
CLINICAL_DIR  <- file.path(DATA_DIR, "clinical")
PROCESSED_DIR <- file.path(DATA_DIR, "processed")

COHORT1_RAW_DIR      <- file.path(RAW_DIR, "cohort1")
COHORT2_RAW_DIR      <- file.path(RAW_DIR, "cohort2")
COHORT1_CLINICAL_DIR <- file.path(CLINICAL_DIR, "cohort1")
COHORT2_CLINICAL_DIR <- file.path(CLINICAL_DIR, "cohort2")

RESULTS_DIR    <- "results"
DEG_DIR        <- file.path(RESULTS_DIR, "differential_expression")
SHARED_DEG_DIR <- file.path(RESULTS_DIR, "shared_deg")
PANEL_DIR      <- file.path(RESULTS_DIR, "gene_panels")
ROC_DIR        <- file.path(RESULTS_DIR, "roc_curves")
BOXPLOT_DIR    <- file.path(RESULTS_DIR, "expression_boxplots")
SURVIVAL_DIR   <- file.path(RESULTS_DIR, "survival_analysis")

# ---- Low-expression filtering ------------------------------------
# A gene is retained if log-CPM > LOW_EXPRESSION_LOGCPM in at least
# LOW_EXPRESSION_FRACTION of samples.

LOW_EXPRESSION_LOGCPM   <- 1
LOW_EXPRESSION_FRACTION <- 0.10

# ---- Differential expression thresholds ---------------------------

DEG_FDR_CUTOFF   <- 0.05
DEG_LOGFC_CUTOFF <- 1.5   # |log2 fold change| >= 1.5, applied per cohort

# Threshold applied when intersecting the two cohorts' DEG results.
# Kept as a separate constant from DEG_LOGFC_CUTOFF (even though both
# currently equal 1.5) because the two steps are conceptually distinct:
# DEG_LOGFC_CUTOFF defines "differentially expressed" within a single
# cohort, SHARED_LOGFC_CUTOFF defines "shared" once both cohorts'
# results are intersected. If you use a lenient per-cohort DEG cutoff
# to explore more candidates before intersecting, set this stricter.
SHARED_LOGFC_CUTOFF <- 1.5
SHARED_FDR_CUTOFF   <- 0.05

# ---- Gene panel selection -------------------------------------------

GENE_RANK_FDR_CUTOFF <- 0.05
PANEL_SIZE           <- 8   # genes selected per direction (up- / down-regulated)

# ---- Survival analysis ------------------------------------------------

# High/low expression groups are defined by the lower/upper score
# quantile; samples strictly between the two quantiles are excluded.
SURVIVAL_LOWER_QUANTILE <- 0.40
SURVIVAL_UPPER_QUANTILE <- 0.60

# Survival-guided panel refinement
SURVIVAL_TOP_N_GENES        <- 8
SURVIVAL_MIN_EVENTS         <- 20
SURVIVAL_SCREENING_P_CUTOFF <- 0.20

# ---- Helpers ------------------------------------------------------------

#' Create a directory (recursively) if it does not already exist.
#'
#' @param path Directory path to create.
#' @return The path, invisibly.
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}
