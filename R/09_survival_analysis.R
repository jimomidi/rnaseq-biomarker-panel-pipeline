# ============================================================
# Step 9 — Survival analysis of the final gene panels
# ============================================================
#
# Input : final gene panels from R/06_gene_panel_selection.R, tumor
#         expression matrices, and per-cohort clinical data
#         (clinical.tsv, as exported from the GDC portal, placed
#         under CLINICAL_DIR/<cohort>/).
# Method:
#   1. Composite-score survival: samples are split into High/Low
#      groups by the SURVIVAL_LOWER_QUANTILE / SURVIVAL_UPPER_QUANTILE
#      of each panel's composite z-score (middle band excluded), then
#      compared with a log-rank test and a univariate Cox model.
#   2. Single-gene survival: the same High/Low split is applied to
#      each panel gene individually (on its own z-score), and the
#      resulting Kaplan-Meier curves are shown side by side per panel.
# Output: KM plots (PNG) and a combined statistics table (CSV), under
#         SURVIVAL_DIR.
#
# Note: cohort 1's composite/single-gene scores are computed directly
# on raw counts (matching the original analysis), while cohort 2 uses
# TMM log-CPM values. This asymmetry is inherited from how each
# cohort's tumor-only matrix was produced upstream (see
# R/01_download_prepare_cohort1.R and R/02_download_prepare_cohort2.R)
# and is noted here rather than silently changed; see the project
# README for details.

source("config.R")

for (pkg in c("survival", "survminer", "ggplot2", "edgeR", "dplyr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "edgeR") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install("edgeR", ask = FALSE)
    } else {
      install.packages(pkg)
    }
  }
}

library(survival)
library(survminer)
library(ggplot2)
library(edgeR)
library(dplyr)
source("R/utils_expression.R")
source("R/utils_survival.R")

set.seed(RANDOM_SEED)
ensure_dir(SURVIVAL_DIR)

#' Build a long data frame of single-gene High/Low survival groups for
#' every gene in `genes`, using the same quantile split as
#' `merge_score_with_clinical()`.
build_single_gene_survival_df <- function(expr_mat, clinical_df, genes, dataset_name) {
  per_gene <- list()
  for (gene in genes) {
    if (!gene %in% rownames(expr_mat)) {
      message("  Gene not found: ", gene, " in ", dataset_name)
      next
    }
    expr <- as.numeric(expr_mat[gene, ])
    names(expr) <- colnames(expr_mat)

    df <- data.frame(sample_id = names(expr), expr = expr, stringsAsFactors = FALSE)
    df$short_id <- shorten_sample_id(df$sample_id)

    merged <- merge(df, clinical_df, by.x = "short_id", by.y = "cases.submitter_id")
    merged <- merged[!is.na(merged$expr) & !is.na(merged$OS_months) & !is.na(merged$OS_event), ]
    if (nrow(merged) < 20) next

    merged$expr_z <- as.numeric(scale(merged$expr))
    q_low  <- quantile(merged$expr_z, SURVIVAL_LOWER_QUANTILE, na.rm = TRUE)
    q_high <- quantile(merged$expr_z, SURVIVAL_UPPER_QUANTILE, na.rm = TRUE)
    merged <- merged[merged$expr_z <= q_low | merged$expr_z >= q_high, ]

    merged$group <- factor(ifelse(merged$expr_z >= q_high, "High", "Low"), levels = c("Low", "High"))
    merged$gene <- gene

    if (length(unique(merged$group)) < 2) next
    if (sum(merged$OS_event == 1) < 3) next

    per_gene[[gene]] <- merged
  }
  if (length(per_gene) == 0) return(data.frame())
  do.call(rbind, per_gene)
}

#' Log-rank + Cox summary for each gene in a single-gene survival data frame.
compute_single_gene_stats <- function(df) {
  genes <- unique(df$gene)
  out <- list()
  for (g in genes) {
    d <- df[df$gene == g, ]
    surv_obj <- Surv(d$OS_months, d$OS_event)
    lr <- survdiff(surv_obj ~ group, data = d)
    logrank_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)

    cox <- tryCatch(coxph(surv_obj ~ group, data = d), error = function(e) NULL)
    if (is.null(cox)) next
    cs <- summary(cox)

    out[[g]] <- data.frame(
      gene = g, N = nrow(d), High_N = sum(d$group == "High"), Low_N = sum(d$group == "Low"),
      Events = sum(d$OS_event == 1), Logrank_P = logrank_p,
      Cox_HR = cs$coefficients[1, "exp(coef)"],
      Cox_HR_low = cs$conf.int[1, "lower .95"],
      Cox_HR_high = cs$conf.int[1, "upper .95"],
      Cox_P = cs$coefficients[1, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

#' Faceted (one panel per gene) Kaplan-Meier plot with per-gene HR/p labels.
make_single_gene_km_plot <- function(df, title, out_file) {
  stats_df <- compute_single_gene_stats(df)
  stats_df$label <- paste0(
    "HR = ", sprintf("%.2f", stats_df$Cox_HR), "\n",
    ifelse(stats_df$Logrank_P < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", stats_df$Logrank_P))),
    "\nN = ", stats_df$N
  )

  curves <- list()
  for (g in unique(df$gene)) {
    d <- df[df$gene == g, ]
    fit <- survfit(Surv(OS_months, OS_event) ~ group, data = d)
    s <- summary(fit)
    tmp <- data.frame(
      time = s$time, surv = s$surv, lower = s$lower, upper = s$upper,
      strata = s$strata, gene = g, stringsAsFactors = FALSE
    )
    tmp$group <- factor(ifelse(grepl("High", tmp$strata), "High", "Low"), levels = c("Low", "High"))
    curves[[g]] <- tmp
  }
  curve_df <- do.call(rbind, curves)

  p <- ggplot(curve_df, aes(x = time, y = surv, color = group, fill = group)) +
    geom_step(linewidth = 0.8) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, color = NA) +
    facet_wrap(~gene, ncol = 4) +
    geom_text(data = stats_df, aes(x = Inf, y = 0.12, label = label), inherit.aes = FALSE,
              hjust = 1.05, vjust = 0, size = 3, fontface = "bold") +
    scale_color_manual(values = c("Low" = "#2E86AB", "High" = "#E84855")) +
    scale_fill_manual(values = c("Low" = "#2E86AB", "High" = "#E84855")) +
    labs(title = title, x = "Time (months)", y = "Overall Survival Probability",
         color = "Expression", fill = "Expression") +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      strip.text = element_text(face = "bold", size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )

  ggsave(file.path(SURVIVAL_DIR, out_file), plot = p, width = 14, height = 8, dpi = 300, bg = "white")
  write.csv(stats_df, file.path(SURVIVAL_DIR, sub("\\.png$", "_stats.csv", out_file)), row.names = FALSE)
}

# ---- Run: load panels, expression, and clinical data ----------------------

panel_up_up     <- read.csv(file.path(PANEL_DIR, "panel_final_up_up.csv"))$gene
panel_down_down <- read.csv(file.path(PANEL_DIR, "panel_final_down_down.csv"))$gene

cohort1_clinical <- load_clinical(COHORT1_CLINICAL_DIR)
cohort2_clinical <- load_clinical(COHORT2_CLINICAL_DIR)

cohort1_tumor_mat <- load_cohort1_tumor_matrix()
cohort2_tumor_mat <- load_cohort2_tumor_matrix()

message("Cohort 1 tumor samples: ", ncol(cohort1_tumor_mat))
message("Cohort 2 tumor samples: ", ncol(cohort2_tumor_mat))

# ---- Composite-score survival ----------------------------------------------

message("\n>>> Composite-score survival analysis ...")
merged_cohort1_up   <- merge_score_with_clinical(cohort1_tumor_mat, cohort1_clinical, panel_up_up,     "Cohort1 Up-Up")
merged_cohort2_up   <- merge_score_with_clinical(cohort2_tumor_mat, cohort2_clinical, panel_up_up,     "Cohort2 Up-Up")
merged_cohort1_down <- merge_score_with_clinical(cohort1_tumor_mat, cohort1_clinical, panel_down_down, "Cohort1 Down-Down")
merged_cohort2_down <- merge_score_with_clinical(cohort2_tumor_mat, cohort2_clinical, panel_down_down, "Cohort2 Down-Down")

stats_df <- do.call(rbind, list(
  compute_survival_stats(merged_cohort1_up,   "Cohort1", "Up-Up"),
  compute_survival_stats(merged_cohort2_up,   "Cohort2", "Up-Up"),
  compute_survival_stats(merged_cohort1_down, "Cohort1", "Down-Down"),
  compute_survival_stats(merged_cohort2_down, "Cohort2", "Down-Down")
))
stats_df$Logrank_FDR <- p.adjust(stats_df$Logrank_P, method = "BH")
stats_df$Cox_FDR     <- p.adjust(stats_df$Cox_P, method = "BH")

write.csv(stats_df, file.path(SURVIVAL_DIR, "km_cox_survival_stats.csv"), row.names = FALSE)
print(stats_df)

p1 <- plot_km_with_hr(merged_cohort1_up,   "KM - Cohort 1 | Up-Up Panel",
                       stats_df[stats_df$Dataset == "Cohort1" & stats_df$Panel == "Up-Up", ])
p2 <- plot_km_with_hr(merged_cohort2_up,   "KM - Cohort 2 | Up-Up Panel",
                       stats_df[stats_df$Dataset == "Cohort2" & stats_df$Panel == "Up-Up", ])
p3 <- plot_km_with_hr(merged_cohort1_down, "KM - Cohort 1 | Down-Down Panel",
                       stats_df[stats_df$Dataset == "Cohort1" & stats_df$Panel == "Down-Down", ])
p4 <- plot_km_with_hr(merged_cohort2_down, "KM - Cohort 2 | Down-Down Panel",
                       stats_df[stats_df$Dataset == "Cohort2" & stats_df$Panel == "Down-Down", ])

save_km_plot(p1, "km_cohort1_up_up.png",     SURVIVAL_DIR)
save_km_plot(p2, "km_cohort2_up_up.png",     SURVIVAL_DIR)
save_km_plot(p3, "km_cohort1_down_down.png", SURVIVAL_DIR)
save_km_plot(p4, "km_cohort2_down_down.png", SURVIVAL_DIR)

# ---- Single-gene survival ---------------------------------------------------

message("\n>>> Single-gene survival analysis ...")
df_cohort1_up   <- build_single_gene_survival_df(cohort1_tumor_mat, cohort1_clinical, panel_up_up,     "Cohort1 Up-Up")
df_cohort2_up   <- build_single_gene_survival_df(cohort2_tumor_mat, cohort2_clinical, panel_up_up,     "Cohort2 Up-Up")
df_cohort1_down <- build_single_gene_survival_df(cohort1_tumor_mat, cohort1_clinical, panel_down_down, "Cohort1 Down-Down")
df_cohort2_down <- build_single_gene_survival_df(cohort2_tumor_mat, cohort2_clinical, panel_down_down, "Cohort2 Down-Down")

make_single_gene_km_plot(df_cohort1_up,   "Single-gene KM - Cohort 1 | Up-Up",     "single_gene_km_cohort1_up_up.png")
make_single_gene_km_plot(df_cohort2_up,   "Single-gene KM - Cohort 2 | Up-Up",     "single_gene_km_cohort2_up_up.png")
make_single_gene_km_plot(df_cohort1_down, "Single-gene KM - Cohort 1 | Down-Down", "single_gene_km_cohort1_down_down.png")
make_single_gene_km_plot(df_cohort2_down, "Single-gene KM - Cohort 2 | Down-Down", "single_gene_km_cohort2_down_down.png")

message("\nSurvival plots and statistics saved in: ", SURVIVAL_DIR)
