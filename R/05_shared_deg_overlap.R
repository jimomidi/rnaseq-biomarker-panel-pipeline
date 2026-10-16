# ============================================================
# Step 5 — Overlap of differentially expressed genes between cohorts
# ============================================================
#
# Input : per-cohort DEG tables produced by
#         R/03_differential_expression_cohort1.R and
#         R/04_differential_expression_cohort2.R.
# Method: apply a (stricter) shared-significance threshold to each
#         cohort independently, then intersect the two gene sets and
#         classify each shared gene into one of four direction
#         patterns: Up_Up, Down_Down, Up_Down, Down_Up.
# Output: one CSV per pattern plus a combined CSV, and a volcano plot
#         per cohort colored by overlap pattern, under SHARED_DEG_DIR.

source("config.R")

library(dplyr)
library(ggplot2)

ensure_dir(SHARED_DEG_DIR)

#' Reduce a DEG table to genes passing the shared-overlap thresholds.
#'
#' @param deg_table DEG table with `gene`, `logFC`, `adj.P.Val` columns.
#' @param suffix Suffix appended to the output column names (e.g. "cohort1").
significant_subset <- function(deg_table, suffix) {
  deg_table %>%
    filter(adj.P.Val < SHARED_FDR_CUTOFF & abs(logFC) >= SHARED_LOGFC_CUTOFF) %>%
    mutate(direction = ifelse(logFC > 0, "Up", "Down")) %>%
    select(
      gene,
      !!paste0("logFC_", suffix) := logFC,
      !!paste0("FDR_", suffix) := adj.P.Val,
      !!paste0("direction_", suffix) := direction
    )
}

#' Classify each row of `df` by an existing category lookup.
add_overlap_category <- function(deg_table, shared, fdr_cutoff, logfc_cutoff) {
  deg_table %>%
    mutate(
      is_significant = adj.P.Val < fdr_cutoff & abs(logFC) >= logfc_cutoff,
      category = case_when(
        gene %in% shared$gene[shared$pattern == "Up_Up"]     ~ "Up_Up",
        gene %in% shared$gene[shared$pattern == "Down_Down"] ~ "Down_Down",
        gene %in% shared$gene[shared$pattern == "Up_Down"]   ~ "Up_Down",
        gene %in% shared$gene[shared$pattern == "Down_Up"]   ~ "Down_Up",
        is_significant                                       ~ "Sig_only",
        TRUE                                                 ~ "NS"
      ),
      category = factor(category, levels = c("Up_Up", "Down_Down", "Up_Down", "Down_Up", "Sig_only", "NS"))
    )
}

#' Volcano plot colored by shared-overlap category.
make_overlap_volcano <- function(df, title, y_break) {
  color_map <- c(
    "Up_Up" = "#E74C3C", "Down_Down" = "#3498DB", "Up_Down" = "#27AE60",
    "Down_Up" = "#9B59B6", "Sig_only" = "grey40", "NS" = "grey80"
  )
  counts <- table(df$category)
  get_n <- function(x) ifelse(is.na(counts[x]), 0, counts[x])
  legend_labels <- c(
    "Up_Up"     = paste0("Up (cohort 1) - Up (cohort 2) (n=",     get_n("Up_Up"),     ")"),
    "Down_Down" = paste0("Down (cohort 1) - Down (cohort 2) (n=", get_n("Down_Down"), ")"),
    "Up_Down"   = paste0("Up (cohort 1) - Down (cohort 2) (n=",   get_n("Up_Down"),   ")"),
    "Down_Up"   = paste0("Down (cohort 1) - Up (cohort 2) (n=",   get_n("Down_Up"),   ")"),
    "Sig_only"  = paste0("Significant, non-shared (n=", get_n("Sig_only"), ")"),
    "NS"        = paste0("Non-significant (n=",         get_n("NS"),       ")")
  )

  ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), color = category)) +
    geom_point(data = filter(df, category == "NS"), size = 0.6, alpha = 0.3) +
    geom_point(data = filter(df, category == "Sig_only"), size = 0.9, alpha = 0.6) +
    geom_point(
      data = filter(df, category %in% c("Up_Up", "Down_Down", "Up_Down", "Down_Up")),
      size = 1.8, alpha = 0.9
    ) +
    geom_vline(xintercept = c(-SHARED_LOGFC_CUTOFF, SHARED_LOGFC_CUTOFF),
               linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_hline(yintercept = -log10(SHARED_FDR_CUTOFF),
               linetype = "dashed", color = "black", linewidth = 0.4) +
    scale_color_manual(values = color_map, labels = legend_labels, name = "Expression Pattern") +
    scale_x_continuous(
      breaks = seq(floor(min(df$logFC, na.rm = TRUE)), ceiling(max(df$logFC, na.rm = TRUE)), by = 2)
    ) +
    scale_y_continuous(
      breaks = seq(0, ceiling(max(-log10(df$adj.P.Val), na.rm = TRUE) / y_break) * y_break, by = y_break)
    ) +
    labs(
      title = title,
      subtitle = paste0("FDR < ", SHARED_FDR_CUTOFF, "  |  |log2FC| >= ", SHARED_LOGFC_CUTOFF),
      x = "log2 Fold Change (Tumor / Normal)",
      y = "-log10 (Adjusted P-value)"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 11),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 9),
      panel.grid.minor = element_blank()
    )
}

# ---- Run -----------------------------------------------------------------

message(">>> Loading per-cohort DEG results ...")
cohort1_deg <- read.csv(file.path(DEG_DIR, "cohort1_deg_results.csv"))
cohort2_deg <- read.csv(file.path(DEG_DIR, "cohort2_deg_results.csv"))

cohort1_sig <- significant_subset(cohort1_deg, "cohort1")
cohort2_sig <- significant_subset(cohort2_deg, "cohort2")
message("  Cohort 1 significant genes: ", nrow(cohort1_sig))
message("  Cohort 2 significant genes: ", nrow(cohort2_sig))

shared <- inner_join(cohort1_sig, cohort2_sig, by = "gene") %>%
  mutate(pattern = paste0(direction_cohort1, "_", direction_cohort2))
message(">>> Total shared genes: ", nrow(shared))

cat("\n========== Shared DEG overlap ==========\n")
cat("Up   in cohort 1 - Up   in cohort 2 (Up_Up)    :", sum(shared$pattern == "Up_Up"),     "\n")
cat("Down in cohort 1 - Down in cohort 2 (Down_Down):", sum(shared$pattern == "Down_Down"), "\n")
cat("Up   in cohort 1 - Down in cohort 2 (Up_Down)  :", sum(shared$pattern == "Up_Down"),   "\n")
cat("Down in cohort 1 - Up   in cohort 2 (Down_Up)  :", sum(shared$pattern == "Down_Up"),   "\n")
cat("==========================================\n")

write.csv(shared, file.path(SHARED_DEG_DIR, "shared_deg_all_patterns.csv"), quote = FALSE, row.names = FALSE)
message("  [saved] shared_deg_all_patterns.csv")

for (pat in c("Up_Up", "Down_Down", "Up_Down", "Down_Up")) {
  sub <- shared %>% filter(pattern == pat)
  if (nrow(sub) > 0) {
    out_path <- file.path(SHARED_DEG_DIR, paste0("shared_deg_", tolower(pat), ".csv"))
    write.csv(sub, out_path, quote = FALSE, row.names = FALSE)
    message("  [saved] ", basename(out_path), " (", nrow(sub), " genes)")
  }
}

# ---- Overlap-colored volcano plots ----------------------------------------

cohort1_plot <- add_overlap_category(cohort1_deg, shared, SHARED_FDR_CUTOFF, SHARED_LOGFC_CUTOFF)
cohort2_plot <- add_overlap_category(cohort2_deg, shared, SHARED_FDR_CUTOFF, SHARED_LOGFC_CUTOFF)

p_cohort1 <- make_overlap_volcano(cohort1_plot, "Volcano Plot - Cohort 1: Tumor vs Normal", y_break = 10)
p_cohort2 <- make_overlap_volcano(cohort2_plot, "Volcano Plot - Cohort 2: Tumor vs Normal", y_break = 20)

ggsave(file.path(SHARED_DEG_DIR, "volcano_cohort1_overlap.png"), p_cohort1, width = 10, height = 7, dpi = 300)
ggsave(file.path(SHARED_DEG_DIR, "volcano_cohort2_overlap.png"), p_cohort2, width = 10, height = 7, dpi = 300)
message("  [saved] volcano_cohort1_overlap.png / volcano_cohort2_overlap.png")
