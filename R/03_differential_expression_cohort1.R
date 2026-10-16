# ============================================================
# Step 3 — Differential expression analysis (cohort 1, tumor vs. normal)
# ============================================================
#
# Input : raw RNA-seq counts for cohort 1 (produced by
#         R/01_download_prepare_cohort1.R).
# Method: TMM normalization + limma-voom.
# Output: DEG table (CSV) and a volcano plot (PNG) under DEG_DIR.

source("config.R")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
for (pkg in c("limma", "edgeR", "ggplot2", "ggrepel", "dplyr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE)
}

library(limma)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(dplyr)

ensure_dir(DEG_DIR)

#' Load a single-assay expression CSV where the gene identifier is the
#' first (unnamed) column, deduplicating repeated gene names.
load_expression_matrix <- function(path) {
  m <- read.csv(path, row.names = NULL, check.names = FALSE)
  m <- m[!duplicated(m[, 1]), ]
  rownames(m) <- m[, 1]
  as.matrix(m[, -1])
}

#' Run limma-voom differential expression between two groups of raw counts.
#'
#' @param merged_counts Raw count matrix (genes x samples), tumor and
#'   normal samples combined.
#' @param group Factor with levels `c("Normal", "Tumor")`, one entry per
#'   column of `merged_counts`.
#' @return A `topTable()`-style data frame with `gene` and `significant`
#'   columns added.
run_voom_deg <- function(merged_counts, group) {
  dge_all <- calcNormFactors(DGEList(counts = merged_counts, group = group), method = "TMM")
  logcpm_all <- cpm(dge_all, log = TRUE, prior.count = 1)

  min_samples <- ceiling(LOW_EXPRESSION_FRACTION * ncol(merged_counts))
  keep <- rowSums(logcpm_all > LOW_EXPRESSION_LOGCPM) >= min_samples

  message("  Genes before low-expression filter: ", nrow(merged_counts))
  message("  Genes after low-expression filter : ", sum(keep))

  dge <- calcNormFactors(dge_all[keep, , keep.lib.sizes = FALSE], method = "TMM")

  design <- model.matrix(~group)
  colnames(design) <- c("Intercept", "TumorVsNormal")

  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)

  results <- topTable(
    fit, coef = "TumorVsNormal", number = Inf,
    sort.by = "P", adjust.method = "BH"
  )
  results$gene <- rownames(results)
  results$significant <- ifelse(
    results$adj.P.Val < DEG_FDR_CUTOFF & abs(results$logFC) >= DEG_LOGFC_CUTOFF,
    ifelse(results$logFC > 0, "Up", "Down"),
    "NS"
  )
  results
}

#' Volcano plot with the top differentially expressed genes labeled.
make_volcano_plot <- function(results, title) {
  top_genes <- results %>%
    filter(significant != "NS") %>%
    arrange(adj.P.Val) %>%
    head(20)

  col_values <- c("Up" = "#E74C3C", "Down" = "#2E86C1", "NS" = "grey70")
  legend_labels <- c(
    "Up"   = paste0("Up (n=",   sum(results$significant == "Up"),   ")"),
    "Down" = paste0("Down (n=", sum(results$significant == "Down"), ")"),
    "NS"   = paste0("NS (n=",   sum(results$significant == "NS"),   ")")
  )

  ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = significant)) +
    geom_point(data = filter(results, significant == "NS"), size = 0.8, alpha = 0.4) +
    geom_point(data = filter(results, significant != "NS"), size = 1.5, alpha = 0.8) +
    geom_vline(xintercept = c(-DEG_LOGFC_CUTOFF, DEG_LOGFC_CUTOFF),
               linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_hline(yintercept = -log10(DEG_FDR_CUTOFF),
               linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_text_repel(
      data = top_genes, aes(label = gene), size = 3, max.overlaps = 20,
      segment.size = 0.3, segment.color = "grey50", show.legend = FALSE
    ) +
    scale_color_manual(values = col_values, labels = legend_labels) +
    labs(
      title = title,
      subtitle = paste0("FDR < ", DEG_FDR_CUTOFF, "  |  |log2FC| >= ", DEG_LOGFC_CUTOFF),
      x = "log2 Fold Change (Tumor / Normal)",
      y = "-log10 (Adjusted P-value)",
      color = "Expression"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 11),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
}

# ---- Run -------------------------------------------------------------------

message(">>> Loading cohort 1 tumor and normal raw counts ...")
tumor_mat  <- load_expression_matrix(file.path(COHORT1_RAW_DIR, "cohort1_tumor_unstranded.csv"))
normal_mat <- load_expression_matrix(file.path(COHORT1_RAW_DIR, "cohort1_normal_unstranded.csv"))
message("  Tumor : ", nrow(tumor_mat), " genes x ", ncol(tumor_mat), " samples")
message("  Normal: ", nrow(normal_mat), " genes x ", ncol(normal_mat), " samples")

common_genes <- intersect(rownames(tumor_mat), rownames(normal_mat))
message(">>> Shared genes: ", length(common_genes))

merged_counts <- cbind(tumor_mat[common_genes, ], normal_mat[common_genes, ])
group <- factor(
  c(rep("Tumor", ncol(tumor_mat)), rep("Normal", ncol(normal_mat))),
  levels = c("Normal", "Tumor")
)

message(">>> Running limma-voom ...")
results <- run_voom_deg(merged_counts, group)

message("\n>>> Differential expression summary:")
message("  Genes tested   : ", nrow(results))
message("  Up-regulated   : ", sum(results$significant == "Up"))
message("  Down-regulated : ", sum(results$significant == "Down"))
message("  Not significant: ", sum(results$significant == "NS"))

write.csv(
  results,
  file = file.path(DEG_DIR, "cohort1_deg_results.csv"),
  quote = FALSE, row.names = FALSE
)
message("  [saved] cohort1_deg_results.csv")

volcano <- make_volcano_plot(results, "Volcano Plot - Cohort 1: Tumor vs Normal")
ggsave(file.path(DEG_DIR, "cohort1_volcano.png"), volcano, width = 9, height = 7, dpi = 200)
message("  [saved] cohort1_volcano.png")
