# ============================================================
# Step 6 — Gene panel selection (single-gene AUC ranking)
# ============================================================
#
# Input : shared-DEG gene lists from R/05_shared_deg_overlap.R
#         (Up_Up and Down_Down candidate genes) and the per-cohort
#         expression matrices (cohort 2's matrix is batch-corrected;
#         see `load_cohort2_expression()` in R/utils_expression.R).
# Method:
#   1. For each candidate gene, compute a single-gene tumor-vs-normal
#      AUC and a Wilcoxon rank-sum test, in each cohort separately.
#      Keep genes with the expected direction of change and FDR below
#      GENE_RANK_FDR_CUTOFF, then rank by AUC (ties broken by
#      |log fold change|, then FDR).
#   2. Walk down both cohorts' ranked lists in lock-step and collect
#      the first PANEL_SIZE genes that appear in both lists — i.e. the
#      top genes that are simultaneously well-ranked in both cohorts.
# Output: four per-cohort AUC ranking tables and two final gene panels
#         (Up_Up, Down_Down), written under PANEL_DIR.

source("config.R")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
for (pkg in c("pROC", "edgeR", "limma")) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE)
}

library(pROC)
library(edgeR)
library(limma)
source("R/utils_expression.R")

set.seed(RANDOM_SEED)
ensure_dir(PANEL_DIR)

#' Rank candidate genes by single-gene discriminative power.
#'
#' @param mat Log-CPM matrix (genes x samples).
#' @param label Character vector of "Tumor"/"Normal" labels.
#' @param candidate_genes Genes to evaluate.
#' @param direction Expected direction of change, "Up" or "Down"
#'   (genes are kept only if their observed direction matches).
#' @param fdr_cutoff FDR threshold (Wilcoxon p-values, BH-adjusted).
#' @return Data frame ranked by AUC (descending), |logFC| (descending),
#'   then adjusted p-value (ascending); includes a `Rank` column.
rank_genes_by_auc <- function(mat, label, candidate_genes, direction, fdr_cutoff) {
  genes <- intersect(candidate_genes, rownames(mat))
  message("  Candidate genes present in matrix: ", length(genes))

  rows <- lapply(genes, function(gene) {
    x <- as.numeric(mat[gene, ])
    tumor_vals  <- x[label == "Tumor"]
    normal_vals <- x[label == "Normal"]

    mean_tumor  <- mean(tumor_vals, na.rm = TRUE)
    mean_normal <- mean(normal_vals, na.rm = TRUE)
    logfc <- mean_tumor - mean_normal

    pval <- tryCatch(
      wilcox.test(tumor_vals, normal_vals, exact = FALSE)$p.value,
      error = function(e) NA
    )

    auc_raw <- tryCatch(
      as.numeric(pROC::auc(pROC::roc(label, x, quiet = TRUE))),
      error = function(e) NA
    )
    # AUC is symmetric around 0.5 depending on label encoding order;
    # report the discriminative strength regardless of that encoding.
    auc_val <- if (!is.na(auc_raw)) max(auc_raw, 1 - auc_raw) else NA

    data.frame(
      gene = gene, AUC = auc_val, P.Value = pval, logFC = logfc,
      mean_Tumor = mean_tumor, mean_Normal = mean_normal,
      stringsAsFactors = FALSE
    )
  })

  df <- do.call(rbind, rows)

  valid <- !is.na(df$P.Value)
  df$adj.P.Val <- NA
  df$adj.P.Val[valid] <- p.adjust(df$P.Value[valid], method = "BH")

  if (direction == "Up")   df <- df[!is.na(df$logFC) & df$logFC > 0, ]
  if (direction == "Down") df <- df[!is.na(df$logFC) & df$logFC < 0, ]

  df <- df[!is.na(df$AUC) & !is.na(df$adj.P.Val) & df$adj.P.Val < fdr_cutoff, ]

  df <- df[order(-df$AUC, -abs(df$logFC), df$adj.P.Val), ]
  df$Rank <- seq_len(nrow(df))
  df[, c("Rank", "gene", "AUC", "P.Value", "adj.P.Val", "logFC", "mean_Tumor", "mean_Normal")]
}

#' Select the top `n_select` genes that appear early in both ranked lists.
#'
#' Walks both lists in lock-step (comparing rank position `i` in each)
#' and keeps the first genes found to appear in both, in the order
#' they are confirmed as shared, stopping once `n_select` are found.
#'
#' @param rank1,rank2 Ranked data frames (as returned by
#'   `rank_genes_by_auc()`) with a `gene` column, best genes first.
#' @param n_select Number of genes to select.
#' @return Character vector of selected gene symbols.
select_shared_top_genes <- function(rank1, rank2, n_select, label = "") {
  message(">>> Selecting shared top genes - ", label)
  genes1 <- rank1$gene
  genes2 <- rank2$gene

  seen1 <- character(0)
  seen2 <- character(0)
  shared <- character(0)
  max_len <- max(length(genes1), length(genes2))

  for (i in seq_len(max_len)) {
    if (i <= length(genes1)) {
      g1 <- genes1[i]
      seen1 <- c(seen1, g1)
      if (g1 %in% seen2 && !g1 %in% shared) {
        shared <- c(shared, g1)
        if (length(shared) == n_select) break
      }
    }
    if (i <= length(genes2)) {
      g2 <- genes2[i]
      seen2 <- c(seen2, g2)
      if (g2 %in% seen1 && !g2 %in% shared) {
        shared <- c(shared, g2)
        if (length(shared) == n_select) break
      }
    }
  }

  if (length(shared) < n_select) {
    message("  Note: only ", length(shared), " shared genes were found.")
  }
  message("  Selected genes: ", paste(shared, collapse = ", "))
  shared
}

# ---- Run -------------------------------------------------------------------

up_genes   <- read.csv(file.path(SHARED_DEG_DIR, "shared_deg_up_up.csv"))$gene
down_genes <- read.csv(file.path(SHARED_DEG_DIR, "shared_deg_down_down.csv"))$gene

message(">>> Loading cohort expression matrices ...")
cohort1 <- load_cohort1_expression()
cohort2 <- load_cohort2_expression()
message("  Cohort 1 - Tumor: ", sum(cohort1$label == "Tumor"), " | Normal: ", sum(cohort1$label == "Normal"))
message("  Cohort 2 - Tumor: ", sum(cohort2$label == "Tumor"), " | Normal: ", sum(cohort2$label == "Normal"))

message("\n>>> Ranking candidate genes by single-gene AUC ...")
rank_cohort1_up   <- rank_genes_by_auc(cohort1$mat, cohort1$label, up_genes,   "Up",   GENE_RANK_FDR_CUTOFF)
rank_cohort1_down <- rank_genes_by_auc(cohort1$mat, cohort1$label, down_genes, "Down", GENE_RANK_FDR_CUTOFF)
rank_cohort2_up   <- rank_genes_by_auc(cohort2$mat, cohort2$label, up_genes,   "Up",   GENE_RANK_FDR_CUTOFF)
rank_cohort2_down <- rank_genes_by_auc(cohort2$mat, cohort2$label, down_genes, "Down", GENE_RANK_FDR_CUTOFF)

write.csv(rank_cohort1_up,   file.path(PANEL_DIR, "rank_cohort1_up.csv"),   row.names = FALSE)
write.csv(rank_cohort1_down, file.path(PANEL_DIR, "rank_cohort1_down.csv"), row.names = FALSE)
write.csv(rank_cohort2_up,   file.path(PANEL_DIR, "rank_cohort2_up.csv"),   row.names = FALSE)
write.csv(rank_cohort2_down, file.path(PANEL_DIR, "rank_cohort2_down.csv"), row.names = FALSE)
message("  [saved] 4 ranking tables")

panel_up_up     <- select_shared_top_genes(rank_cohort1_up,   rank_cohort2_up,   PANEL_SIZE, "Up_Up")
panel_down_down <- select_shared_top_genes(rank_cohort1_down, rank_cohort2_down, PANEL_SIZE, "Down_Down")

write.csv(data.frame(gene = panel_up_up),     file.path(PANEL_DIR, "panel_final_up_up.csv"),     row.names = FALSE)
write.csv(data.frame(gene = panel_down_down), file.path(PANEL_DIR, "panel_final_down_down.csv"), row.names = FALSE)

message("\n========== Final gene panels ==========")
message("Up_Up panel:")
print(data.frame(
  gene = panel_up_up,
  AUC_cohort1 = round(rank_cohort1_up$AUC[match(panel_up_up, rank_cohort1_up$gene)], 4),
  AUC_cohort2 = round(rank_cohort2_up$AUC[match(panel_up_up, rank_cohort2_up$gene)], 4)
))
message("Down_Down panel:")
print(data.frame(
  gene = panel_down_down,
  AUC_cohort1 = round(rank_cohort1_down$AUC[match(panel_down_down, rank_cohort1_down$gene)], 4),
  AUC_cohort2 = round(rank_cohort2_down$AUC[match(panel_down_down, rank_cohort2_down$gene)], 4)
))
message("Output directory: ", PANEL_DIR)
