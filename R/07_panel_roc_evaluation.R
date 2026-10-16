# ============================================================
# Step 7 — ROC evaluation of the final gene panels
# ============================================================
#
# Input : final gene panels from R/06_gene_panel_selection.R and the
#         per-cohort expression matrices (cohort 2's matrix is
#         batch-corrected; see `load_cohort2_expression()` in
#         R/utils_expression.R).
# Method: for each panel and cohort, compute a composite score as the
#         mean per-gene z-score across the panel, then evaluate its
#         ability to discriminate tumor from normal samples with ROC/AUC.
# Output: one square ROC plot per panel x cohort combination (4 total),
#         under ROC_DIR.

source("config.R")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
for (pkg in c("pROC", "ggplot2", "edgeR", "limma")) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE)
}

library(pROC)
library(ggplot2)
library(edgeR)
library(limma)
source("R/utils_expression.R")

set.seed(RANDOM_SEED)
ensure_dir(ROC_DIR)

#' Composite panel score = mean per-gene z-score across `genes`.
composite_panel_score <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  sub <- mat[g, , drop = FALSE]
  z <- t(apply(sub, 1, function(x) (x - mean(x)) / (sd(x) + 1e-10)))
  colMeans(z, na.rm = TRUE)
}

#' ROC curve for a panel's composite score against tumor/normal labels.
panel_roc <- function(mat, label, genes) {
  score <- composite_panel_score(mat, genes)
  pROC::roc(label, score, quiet = TRUE)
}

#' Square ROC plot annotated with AUC, 95% CI, and the Youden-optimal
#' sensitivity/specificity.
plot_roc_square <- function(roc_obj, title, n_genes, n_tumor, n_normal) {
  auc_val <- as.numeric(pROC::auc(roc_obj))
  ci_val  <- pROC::ci.auc(roc_obj, conf.level = 0.95)

  coords <- pROC::coords(
    roc_obj, "best", best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity")
  )
  best_sensitivity <- round(as.numeric(coords["sensitivity"]), 3)
  best_specificity <- round(as.numeric(coords["specificity"]), 3)

  roc_df <- data.frame(FPR = 1 - roc_obj$specificities, TPR = roc_obj$sensitivities)
  roc_df <- roc_df[order(roc_df$FPR, roc_df$TPR), ]

  info_label <- paste0(
    "AUC = ", round(auc_val, 4),
    "\n95% CI: ", round(ci_val[1], 3), " - ", round(ci_val[3], 3),
    "\nSensitivity = ", best_sensitivity,
    "\nSpecificity = ", best_specificity,
    "\nTumor n = ", n_tumor,
    "\nNormal n = ", n_normal,
    "\nGenes in panel = ", n_genes
  )

  ggplot(roc_df, aes(x = FPR, y = TPR)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.6) +
    geom_path(color = "#C0392B", linewidth = 1.4) +
    geom_point(aes(x = 1 - best_specificity, y = best_sensitivity), color = "#2C3E50", size = 3, shape = 18) +
    annotate("text", x = 0.52, y = 0.45, label = info_label,
             hjust = 0, vjust = 1, size = 3.8, color = "black", fontface = "bold", family = "mono") +
    scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1), expand = c(0, 0)) +
    labs(title = title, x = "False Positive Rate (1 - Specificity)", y = "True Positive Rate (Sensitivity)") +
    coord_fixed(ratio = 1) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11), aspect.ratio = 1)
}

save_square_plot <- function(plot_obj, filename) {
  ggsave(
    filename = file.path(ROC_DIR, filename), plot = plot_obj, device = "png",
    width = 5, height = 5, units = "in", dpi = 300, bg = "white", limitsize = FALSE
  )
}

# ---- Run -------------------------------------------------------------------

panel_up_up     <- read.csv(file.path(PANEL_DIR, "panel_final_up_up.csv"))$gene
panel_down_down <- read.csv(file.path(PANEL_DIR, "panel_final_down_down.csv"))$gene

message(">>> Loading cohort expression matrices ...")
cohort1 <- load_cohort1_expression()
cohort2 <- load_cohort2_expression()

roc_cohort1_up   <- panel_roc(cohort1$mat, cohort1$label, panel_up_up)
roc_cohort2_up   <- panel_roc(cohort2$mat, cohort2$label, panel_up_up)
roc_cohort1_down <- panel_roc(cohort1$mat, cohort1$label, panel_down_down)
roc_cohort2_down <- panel_roc(cohort2$mat, cohort2$label, panel_down_down)

p1 <- plot_roc_square(roc_cohort1_up,   "ROC - Cohort 1 | Up-Up Panel",
                       length(panel_up_up),   sum(cohort1$label == "Tumor"), sum(cohort1$label == "Normal"))
p2 <- plot_roc_square(roc_cohort2_up,   "ROC - Cohort 2 | Up-Up Panel",
                       length(panel_up_up),   sum(cohort2$label == "Tumor"), sum(cohort2$label == "Normal"))
p3 <- plot_roc_square(roc_cohort1_down, "ROC - Cohort 1 | Down-Down Panel",
                       length(panel_down_down), sum(cohort1$label == "Tumor"), sum(cohort1$label == "Normal"))
p4 <- plot_roc_square(roc_cohort2_down, "ROC - Cohort 2 | Down-Down Panel",
                       length(panel_down_down), sum(cohort2$label == "Tumor"), sum(cohort2$label == "Normal"))

save_square_plot(p1, "roc_cohort1_up_up.png")
save_square_plot(p2, "roc_cohort2_up_up.png")
save_square_plot(p3, "roc_cohort1_down_down.png")
save_square_plot(p4, "roc_cohort2_down_down.png")

message("\n========== ROC summary ==========")
message("Up-Up   cohort 1 AUC: ", round(as.numeric(auc(roc_cohort1_up)), 4))
message("Up-Up   cohort 2 AUC: ", round(as.numeric(auc(roc_cohort2_up)), 4))
message("Down-Down cohort 1 AUC: ", round(as.numeric(auc(roc_cohort1_down)), 4))
message("Down-Down cohort 2 AUC: ", round(as.numeric(auc(roc_cohort2_down)), 4))
message("Plots saved in: ", ROC_DIR)
