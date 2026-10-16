# ============================================================
# Step 8 — Expression boxplots for the final gene panels
# ============================================================
#
# Input : final gene panels from R/06_gene_panel_selection.R and the
#         per-cohort expression matrices (cohort 2's matrix is
#         batch-corrected; see `load_cohort2_expression()` in
#         R/utils_expression.R).
# Method: for each panel, draw per-gene, per-cohort boxplots of
#         tumor vs. normal expression, annotated with a
#         significance marker from a Wilcoxon rank-sum test.
# Output: one faceted boxplot PNG per panel, under BOXPLOT_DIR.

source("config.R")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
for (pkg in c("ggplot2", "edgeR", "limma")) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE)
}

library(ggplot2)
library(edgeR)
library(limma)
source("R/utils_expression.R")

set.seed(RANDOM_SEED)
ensure_dir(BOXPLOT_DIR)

#' Map a p-value to a conventional significance marker.
pvalue_to_stars <- function(p) {
  if (is.na(p))        return("ns")
  else if (p < 0.0001) return("****")
  else if (p < 0.001)  return("***")
  else if (p < 0.01)   return("**")
  else if (p < 0.05)   return("*")
  else                 return("ns")
}

#' Build a long-format data frame (gene x cohort x group) for plotting.
make_long_expression_df <- function(mat1, label1, mat2, label2, genes,
                                     cohort1_name = "Cohort1", cohort2_name = "Cohort2") {
  rows <- list()
  for (gene in genes) {
    if (gene %in% rownames(mat1)) {
      rows[[length(rows) + 1]] <- data.frame(
        gene = gene, expr = as.numeric(mat1[gene, ]), group = label1, cohort = cohort1_name
      )
    }
    if (gene %in% rownames(mat2)) {
      rows[[length(rows) + 1]] <- data.frame(
        gene = gene, expr = as.numeric(mat2[gene, ]), group = label2, cohort = cohort2_name
      )
    }
  }
  df <- do.call(rbind, rows)
  df$gene   <- factor(df$gene, levels = genes)
  df$group  <- factor(df$group, levels = c("Normal", "Tumor"))
  df$cohort <- factor(df$cohort, levels = c(cohort1_name, cohort2_name))
  df
}

#' Per gene x cohort Wilcoxon test and label placement for significance stars.
compute_significance_stars <- function(df) {
  combos <- unique(df[, c("gene", "cohort")])
  rows <- lapply(seq_len(nrow(combos)), function(i) {
    g <- as.character(combos$gene[i])
    co <- as.character(combos$cohort[i])
    sub <- df[df$gene == g & df$cohort == co, ]
    tumor_vals  <- sub$expr[sub$group == "Tumor"]
    normal_vals <- sub$expr[sub$group == "Normal"]
    pval <- tryCatch(wilcox.test(tumor_vals, normal_vals, exact = FALSE)$p.value, error = function(e) NA)

    y_max <- max(sub$expr, na.rm = TRUE)
    y_min <- min(sub$expr, na.rm = TRUE)
    span <- y_max - y_min

    data.frame(
      gene = g, cohort = co, pval = pval, star = pvalue_to_stars(pval),
      y_star = y_max + span * 0.08, y_max = y_max, y_min = y_min,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Faceted (cohort x gene) tumor-vs-normal boxplot with significance stars.
make_panel_boxplot <- function(mat1, label1, mat2, label2, genes, title,
                                cohort1_name = "Cohort1", cohort2_name = "Cohort2") {
  df    <- make_long_expression_df(mat1, label1, mat2, label2, genes, cohort1_name, cohort2_name)
  stars <- compute_significance_stars(df)
  stars$gene   <- factor(stars$gene, levels = genes)
  stars$cohort <- factor(stars$cohort, levels = c(cohort1_name, cohort2_name))

  ggplot(df, aes(x = group, y = expr, fill = group)) +
    stat_boxplot(geom = "errorbar", width = 0.3, linewidth = 0.5) +
    geom_boxplot(outlier.size = 0, outlier.alpha = 0, linewidth = 0.5, width = 0.6) +
    geom_jitter(width = 0.18, size = 1.0, alpha = 0.35, color = "gray25") +
    geom_text(data = stars, aes(x = 1.5, y = y_star, label = star),
              inherit.aes = FALSE, size = 5, fontface = "bold", color = "black") +
    facet_grid(cohort ~ gene, scales = "free_y") +
    scale_fill_manual(values = c("Normal" = "#5B9BD5", "Tumor" = "#ED7D31"), name = "Group") +
    labs(title = title, x = NULL, y = "Expression (log-CPM)") +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text.x = element_text(face = "bold", size = 9),
      strip.text.y = element_text(face = "bold", size = 10),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}

# ---- Run -------------------------------------------------------------------

panel_up_up     <- read.csv(file.path(PANEL_DIR, "panel_final_up_up.csv"))$gene
panel_down_down <- read.csv(file.path(PANEL_DIR, "panel_final_down_down.csv"))$gene

message(">>> Loading cohort expression matrices ...")
cohort1 <- load_cohort1_expression()
cohort2 <- load_cohort2_expression()

message(">>> Drawing Up-Up panel boxplot ...")
p_up <- make_panel_boxplot(
  cohort1$mat, cohort1$label, cohort2$mat, cohort2$label,
  genes = panel_up_up, title = "Expression Panel - Up-Up Genes"
)
ggsave(file.path(BOXPLOT_DIR, "boxplot_up_up.png"), p_up, width = 16, height = 6, dpi = 150, device = "png", bg = "white")

message(">>> Drawing Down-Down panel boxplot ...")
p_down <- make_panel_boxplot(
  cohort1$mat, cohort1$label, cohort2$mat, cohort2$label,
  genes = panel_down_down, title = "Expression Panel - Down-Down Genes"
)
ggsave(file.path(BOXPLOT_DIR, "boxplot_down_down.png"), p_down, width = 16, height = 6, dpi = 150, device = "png", bg = "white")

message("Plots saved in: ", BOXPLOT_DIR)
