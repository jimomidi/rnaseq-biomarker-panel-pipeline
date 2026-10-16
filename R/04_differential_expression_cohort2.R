# ============================================================
# Step 4 — Differential expression analysis (cohort 2, tumor vs. normal)
# ============================================================
#
# Input : cohort 2's batch-corrected, low-expression-filtered log-CPM
#         matrix, loaded via `load_cohort2_expression()`
#         (R/utils_expression.R). That helper merges the tumor (GDC)
#         and external control (recount3) expression, then removes the
#         study-of-origin batch effect with `limma::removeBatchEffect()`
#         while protecting the tumor/normal group as a covariate — the
#         same correction is reused by gene panel selection (R/06),
#         ROC evaluation (R/07), and boxplots (R/08), so every analysis
#         that compares tumor vs. normal for cohort 2 sees the same
#         batch-corrected values.
# Method: limma directly on the (already log-transformed, already
#         batch-corrected) log-CPM matrix; voom is not used here since
#         the input is not raw counts.
# Output: DEG table (CSV), volcano plot (PNG), and a PCA plot (PNG) of
#         the batch-corrected matrix, under DEG_DIR.

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
source("R/utils_expression.R")

ensure_dir(DEG_DIR)

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

#' PCA scatter plot of the batch-corrected matrix, colored by group.
make_pca_plot <- function(mat, meta, title) {
  mat <- mat[apply(mat, 1, var) > 0, ]
  pca <- prcomp(t(mat), scale. = TRUE)
  var_explained <- round(summary(pca)$importance[2, 1:2] * 100, 1)

  df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Group = meta$group, Batch = meta$batch)

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

# ---- Run ---------------------------------------------------------------

message(">>> Loading cohort 2 batch-corrected expression matrix ...")
cohort2 <- load_cohort2_expression()
message("  Expression matrix: ", nrow(cohort2$mat), " genes x ", ncol(cohort2$mat), " samples")

group <- factor(cohort2$label, levels = c("Normal", "Tumor"))
design <- model.matrix(~group)

message(">>> Running limma ...")
fit <- lmFit(cohort2$mat, design)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)

results <- topTable(
  fit, coef = "groupTumor", number = Inf,
  sort.by = "P", adjust.method = "BH"
)
results$gene <- rownames(results)
results$significant <- ifelse(
  results$adj.P.Val < DEG_FDR_CUTOFF & abs(results$logFC) >= DEG_LOGFC_CUTOFF,
  ifelse(results$logFC > 0, "Up", "Down"),
  "NS"
)

message("\n>>> Differential expression summary:")
message("  Genes tested   : ", nrow(results))
message("  Up-regulated   : ", sum(results$significant == "Up"))
message("  Down-regulated : ", sum(results$significant == "Down"))
message("  Not significant: ", sum(results$significant == "NS"))

write.csv(
  results,
  file = file.path(DEG_DIR, "cohort2_deg_results.csv"),
  quote = FALSE, row.names = FALSE
)
message("  [saved] cohort2_deg_results.csv")

volcano <- make_volcano_plot(results, "Volcano Plot - Cohort 2: Tumor vs Normal")
ggsave(file.path(DEG_DIR, "cohort2_volcano.png"), volcano, width = 9, height = 7, dpi = 200)
message("  [saved] cohort2_volcano.png")

pca_plot <- make_pca_plot(cohort2$mat, cohort2$meta, "PCA - After removeBatchEffect")
ggsave(file.path(DEG_DIR, "cohort2_pca_batch_adjusted.png"), pca_plot, width = 7, height = 5, dpi = 200)
message("  [saved] cohort2_pca_batch_adjusted.png")
