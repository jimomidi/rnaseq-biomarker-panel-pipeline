# ============================================================
# Dependency installation for the RNA-seq biomarker panel pipeline
# ============================================================
#
# Run once, from the repository root:
#   Rscript install_packages.R
#
# Packages are grouped by source (CRAN vs. Bioconductor) based on
# where each is published. BiocManager::install() can install both
# types, so it is used throughout for consistency.

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_packages <- c(
  "dplyr",
  "ggplot2",
  "ggrepel",
  "survival",
  "survminer",
  "pROC",
  "vegan"
)

bioc_packages <- c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "recount3",
  "edgeR",
  "limma"
)

all_packages <- c(cran_packages, bioc_packages)
missing_packages <- all_packages[!vapply(all_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0) {
  message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
  BiocManager::install(missing_packages, ask = FALSE, update = FALSE)
} else {
  message("All required packages are already installed.")
}
