# ============================================================
# Step 10 — Survival-guided panel refinement and covariate adjustment
# ============================================================
#
# This is a secondary, exploratory extension of R/09_survival_analysis.R:
# rather than assuming the diagnostic panels from
# R/06_gene_panel_selection.R are also the best *prognostic* markers,
# it screens each panel gene individually for association with
# survival and re-weights a risk score accordingly. It then checks
# whether the composite score panels remain associated with survival
# after adjusting for available clinical covariates.
#
# Method:
#   1. Single-gene screening: fit a univariate Cox model per gene
#      (z-scored expression) and keep genes with FDR/effect size
#      passing a lenient screening threshold.
#   2. Risk score: combine the top-screened genes into a single score,
#      signing each gene's z-score by the direction of its hazard
#      ratio (HR >= 1 contributes positively, HR < 1 negatively), so
#      that a higher score always indicates higher predicted risk.
#      Samples are then split into High/Low risk groups by the score
#      median (unlike the panel-level 40/60 split used in step 9).
#   3. Multivariable Cox regression: test the composite score (from
#      step 9) and the survival-guided risk score (from step 2) as
#      continuous predictors, with and without adjustment for age,
#      gender, and tumor stage (whichever are available per cohort).
#
# Output: single-gene Cox screening tables, selected top genes,
# survival-guided KM plots, and adjusted Cox regression tables, under
# SURVIVAL_DIR / SURVIVAL_DIR/survival_guided / SURVIVAL_DIR/adjusted_cox.

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
guided_dir <- file.path(SURVIVAL_DIR, "survival_guided")
adjusted_dir <- file.path(SURVIVAL_DIR, "adjusted_cox")
ensure_dir(guided_dir)
ensure_dir(adjusted_dir)

# ---- 1. Single-gene Cox screening ------------------------------------------

#' Univariate Cox regression (z-scored expression) for each candidate
#' gene, keeping only genes with enough events in the merged cohort.
single_gene_cox_screen <- function(expr_mat, clinical_df, candidate_genes, dataset_name, panel_name) {
  genes <- intersect(candidate_genes, rownames(expr_mat))
  if (length(genes) == 0) {
    message("No candidate genes found for ", dataset_name, " ", panel_name)
    return(NULL)
  }

  expr_df <- as.data.frame(t(expr_mat[genes, , drop = FALSE]))
  expr_df$sample_id <- rownames(expr_df)
  expr_df$short_id  <- shorten_sample_id(expr_df$sample_id)

  merged <- merge(expr_df, clinical_df, by.x = "short_id", by.y = "cases.submitter_id")
  merged <- merged[!is.na(merged$OS_months) & !is.na(merged$OS_event) & merged$OS_months > 0, ]

  if (nrow(merged) == 0) {
    message("No matched samples for ", dataset_name, " ", panel_name)
    return(NULL)
  }
  if (sum(merged$OS_event == 1) < SURVIVAL_MIN_EVENTS) {
    message("Too few events for ", dataset_name, " ", panel_name)
    return(NULL)
  }

  rows <- list()
  for (g in genes) {
    x <- suppressWarnings(as.numeric(merged[[g]]))
    if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) next

    merged$gene_expr_z <- as.numeric(scale(x))
    fit <- tryCatch(coxph(Surv(OS_months, OS_event) ~ gene_expr_z, data = merged), error = function(e) NULL)
    if (is.null(fit)) next
    s <- summary(fit)

    rows[[g]] <- data.frame(
      Dataset = dataset_name, Panel = panel_name, Gene = g, N = nrow(merged),
      Events = sum(merged$OS_event == 1),
      HR = s$coefficients[1, "exp(coef)"],
      HR_low = s$conf.int[1, "lower .95"],
      HR_high = s$conf.int[1, "upper .95"],
      Cox_P = s$coefficients[1, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
  }

  res <- do.call(rbind, rows)
  if (is.null(res) || nrow(res) == 0) return(NULL)

  res$FDR <- p.adjust(res$Cox_P, method = "BH")
  res[order(res$Cox_P, -abs(log(res$HR))), ]
}

#' Keep the top `top_n` screened genes, preferring those below
#' `p_cutoff` (falling back to the full table if none pass).
select_top_survival_genes <- function(cox_table, top_n = SURVIVAL_TOP_N_GENES,
                                       p_cutoff = SURVIVAL_SCREENING_P_CUTOFF) {
  if (is.null(cox_table) || nrow(cox_table) == 0) return(NULL)

  candidates <- cox_table[!is.na(cox_table$Cox_P), ]
  below_cutoff <- candidates[candidates$Cox_P <= p_cutoff, ]
  if (nrow(below_cutoff) == 0) below_cutoff <- candidates

  below_cutoff <- below_cutoff[order(below_cutoff$Cox_P, -abs(log(below_cutoff$HR))), ]
  head(below_cutoff, top_n)
}

# ---- 2. Survival-guided risk score -----------------------------------------

#' Directional composite risk score: each gene's z-score is signed by
#' its hazard-ratio direction (HR >= 1 -> +1, HR < 1 -> -1) so that a
#' higher score always corresponds to higher predicted risk.
compute_survival_risk_score <- function(expr_mat, selected_table) {
  genes <- intersect(selected_table$Gene, rownames(expr_mat))
  if (length(genes) == 0) return(NULL)
  selected_table <- selected_table[selected_table$Gene %in% genes, ]

  sub <- expr_mat[selected_table$Gene, , drop = FALSE]
  z <- t(apply(sub, 1, function(x) {
    s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(NA, length(x)))
    (x - mean(x, na.rm = TRUE)) / s
  }))

  direction <- setNames(ifelse(selected_table$HR >= 1, 1, -1), selected_table$Gene)
  z_signed <- z * direction[rownames(z)]
  score <- colMeans(z_signed, na.rm = TRUE)

  data.frame(sample_id = names(score), score = as.numeric(score), stringsAsFactors = FALSE)
}

#' Merge a survival-guided risk score with clinical data and split into
#' High/Low groups by the score median (a straight median split, unlike
#' the 40/60 quantile split used for the diagnostic panels in step 9).
merge_survival_guided_score <- function(expr_mat, clinical_df, selected_table, dataset_name, panel_name) {
  score_df <- compute_survival_risk_score(expr_mat, selected_table)
  if (is.null(score_df)) {
    message("No score generated for ", dataset_name, " ", panel_name)
    return(NULL)
  }
  score_df$short_id <- shorten_sample_id(score_df$sample_id)

  merged <- merge(score_df, clinical_df, by.x = "short_id", by.y = "cases.submitter_id")
  merged <- merged[
    !is.na(merged$OS_months) & !is.na(merged$OS_event) &
      !is.na(merged$score) & merged$OS_months > 0,
  ]

  med <- median(merged$score, na.rm = TRUE)
  merged$group <- factor(ifelse(merged$score >= med, "High", "Low"), levels = c("Low", "High"))

  message(sprintf(
    "%s %s survival-guided - N: %d | High: %d | Low: %d | Events: %d",
    dataset_name, panel_name, nrow(merged),
    sum(merged$group == "High"), sum(merged$group == "Low"), sum(merged$OS_event == 1)
  ))
  merged
}

#' KM plot for a survival-guided risk score (computes its own log-rank
#' and Cox summary, unlike `plot_km_with_hr()` which takes them as input).
plot_km_survival_guided <- function(merged_df, title) {
  if (is.null(merged_df) || nrow(merged_df) == 0) return(NULL)

  fit <- survfit(Surv(OS_months, OS_event) ~ group, data = merged_df)
  lr <- survdiff(Surv(OS_months, OS_event) ~ group, data = merged_df)
  logrank_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)

  cox <- coxph(Surv(OS_months, OS_event) ~ group, data = merged_df)
  cox_summary <- summary(cox)
  hr <- cox_summary$coefficients[1, "exp(coef)"]
  hr_low <- cox_summary$conf.int[1, "lower .95"]
  hr_high <- cox_summary$conf.int[1, "upper .95"]

  info_text <- paste0(
    format_pvalue(logrank_p),
    "\nHR = ", sprintf("%.2f", hr), " (95% CI: ", sprintf("%.2f", hr_low), "-", sprintf("%.2f", hr_high), ")",
    "\nHigh n = ", sum(merged_df$group == "High"), " | Low n = ", sum(merged_df$group == "Low")
  )

  p_obj <- ggsurvplot(
    fit, data = merged_df, palette = c("#2E86AB", "#E84855"),
    conf.int = TRUE, conf.int.alpha = 0.12, pval = FALSE, risk.table = FALSE,
    xlab = "Time (months)", ylab = "Overall Survival Probability",
    title = title, legend.title = "Risk score", legend.labs = c("Low", "High"),
    ggtheme = theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.78, 0.90),
        legend.background = element_rect(fill = "white", color = "gray80")
      ),
    fontsize = 3.8
  )
  p_obj$plot <- p_obj$plot +
    annotate("text", x = Inf, y = Inf, label = info_text,
             hjust = 1.05, vjust = 1.15, size = 3.3, color = "black", fontface = "bold")
  p_obj
}

# ---- 3. Multivariable Cox with clinical covariates -------------------------

#' Derive age, gender, and AJCC-style stage-group covariates from
#' whichever standard GDC clinical columns are present.
prepare_clinical_covariates <- function(df) {
  if ("demographic.age_at_index" %in% colnames(df)) {
    df$age <- suppressWarnings(as.numeric(df$demographic.age_at_index))
  } else if ("diagnoses.age_at_diagnosis" %in% colnames(df)) {
    df$age <- suppressWarnings(as.numeric(df$diagnoses.age_at_diagnosis)) / 365.25
  } else {
    df$age <- NA
  }

  df$gender <- if ("demographic.gender" %in% colnames(df)) factor(df$demographic.gender) else NA

  stage_cols <- c("diagnoses.ajcc_pathologic_stage", "diagnoses.ajcc_clinical_stage", "diagnoses.tumor_stage")
  stage_col <- stage_cols[stage_cols %in% colnames(df)][1]
  df$stage_raw <- toupper(as.character(if (!is.na(stage_col)) df[[stage_col]] else NA))

  df$stage_group <- dplyr::case_when(
    grepl("STAGE I[^V]", df$stage_raw)  ~ "I",
    grepl("STAGE II[^I]", df$stage_raw) ~ "II",
    grepl("STAGE III", df$stage_raw)    ~ "III",
    grepl("STAGE IV", df$stage_raw)     ~ "IV",
    TRUE ~ NA_character_
  )
  df$stage_group <- factor(df$stage_group, levels = c("I", "II", "III", "IV"))
  df
}

#' Fit a continuous-score Cox model, then (if enough covariate data is
#' available) a second model adjusted for age/gender/stage.
run_adjusted_cox <- function(merged_df, dataset_name, panel_name) {
  if (is.null(merged_df) || nrow(merged_df) == 0) return(NULL)

  df <- prepare_clinical_covariates(merged_df)
  df$score_z <- as.numeric(scale(df$score))

  fit_raw <- coxph(Surv(OS_months, OS_event) ~ score_z, data = df)
  raw_summary <- summary(fit_raw)
  raw_res <- data.frame(
    Dataset = dataset_name, Panel = panel_name, Model = "Continuous score only",
    N = raw_summary$n, Events = raw_summary$nevent, Variable = "score_z",
    HR = raw_summary$coefficients[1, "exp(coef)"],
    HR_low = raw_summary$conf.int[1, "lower .95"],
    HR_high = raw_summary$conf.int[1, "upper .95"],
    P_value = raw_summary$coefficients[1, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )

  covariates <- character(0)
  if (sum(!is.na(df$age)) > 30) covariates <- c(covariates, "age")
  if (length(unique(na.omit(df$gender))) >= 2) covariates <- c(covariates, "gender")
  if (length(unique(na.omit(df$stage_group))) >= 2) covariates <- c(covariates, "stage_group")

  if (length(covariates) == 0) return(raw_res)

  formula_txt <- paste("Surv(OS_months, OS_event) ~ score_z +", paste(covariates, collapse = " + "))
  fit_adj <- tryCatch(coxph(as.formula(formula_txt), data = df), error = function(e) NULL)
  if (is.null(fit_adj)) return(raw_res)

  adj_summary <- summary(fit_adj)
  adj_res <- data.frame(
    Dataset = dataset_name, Panel = panel_name,
    Model = paste0("Adjusted: ", paste(covariates, collapse = " + ")),
    N = adj_summary$n, Events = adj_summary$nevent,
    Variable = rownames(adj_summary$coefficients),
    HR = adj_summary$coefficients[, "exp(coef)"],
    HR_low = adj_summary$conf.int[, "lower .95"],
    HR_high = adj_summary$conf.int[, "upper .95"],
    P_value = adj_summary$coefficients[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )

  rbind(raw_res, adj_res)
}

# ---- Run --------------------------------------------------------------------

panel_up_up     <- read.csv(file.path(PANEL_DIR, "panel_final_up_up.csv"))$gene
panel_down_down <- read.csv(file.path(PANEL_DIR, "panel_final_down_down.csv"))$gene

cohort1_clinical <- load_clinical(COHORT1_CLINICAL_DIR)
cohort2_clinical <- load_clinical(COHORT2_CLINICAL_DIR)
cohort1_tumor_mat <- load_cohort1_tumor_matrix()
cohort2_tumor_mat <- load_cohort2_tumor_matrix()

message(">>> Single-gene Cox screening ...")
cox_cohort1_up   <- single_gene_cox_screen(cohort1_tumor_mat, cohort1_clinical, panel_up_up,     "Cohort1", "Up-Up")
cox_cohort2_up   <- single_gene_cox_screen(cohort2_tumor_mat, cohort2_clinical, panel_up_up,     "Cohort2", "Up-Up")
cox_cohort1_down <- single_gene_cox_screen(cohort1_tumor_mat, cohort1_clinical, panel_down_down, "Cohort1", "Down-Down")
cox_cohort2_down <- single_gene_cox_screen(cohort2_tumor_mat, cohort2_clinical, panel_down_down, "Cohort2", "Down-Down")

all_single_cox <- do.call(rbind, list(cox_cohort1_up, cox_cohort2_up, cox_cohort1_down, cox_cohort2_down))
write.csv(all_single_cox, file.path(guided_dir, "single_gene_cox_survival_results.csv"), row.names = FALSE)

top_cohort1_up   <- select_top_survival_genes(cox_cohort1_up)
top_cohort2_up   <- select_top_survival_genes(cox_cohort2_up)
top_cohort1_down <- select_top_survival_genes(cox_cohort1_down)
top_cohort2_down <- select_top_survival_genes(cox_cohort2_down)

write.csv(top_cohort1_up,   file.path(guided_dir, "top_survival_genes_cohort1_up_up.csv"),     row.names = FALSE)
write.csv(top_cohort2_up,   file.path(guided_dir, "top_survival_genes_cohort2_up_up.csv"),     row.names = FALSE)
write.csv(top_cohort1_down, file.path(guided_dir, "top_survival_genes_cohort1_down_down.csv"), row.names = FALSE)
write.csv(top_cohort2_down, file.path(guided_dir, "top_survival_genes_cohort2_down_down.csv"), row.names = FALSE)

message(">>> Building survival-guided risk scores and KM plots ...")
merged_cohort1_up_guided   <- merge_survival_guided_score(cohort1_tumor_mat, cohort1_clinical, top_cohort1_up,   "Cohort1", "Up-Up")
merged_cohort2_up_guided   <- merge_survival_guided_score(cohort2_tumor_mat, cohort2_clinical, top_cohort2_up,   "Cohort2", "Up-Up")
merged_cohort1_down_guided <- merge_survival_guided_score(cohort1_tumor_mat, cohort1_clinical, top_cohort1_down, "Cohort1", "Down-Down")
merged_cohort2_down_guided <- merge_survival_guided_score(cohort2_tumor_mat, cohort2_clinical, top_cohort2_down, "Cohort2", "Down-Down")

save_km_plot(plot_km_survival_guided(merged_cohort1_up_guided,   "KM - Cohort 1 | Up-Up Survival-Guided Score"),
             "km_cohort1_up_up_survival_guided.png", guided_dir)
save_km_plot(plot_km_survival_guided(merged_cohort2_up_guided,   "KM - Cohort 2 | Up-Up Survival-Guided Score"),
             "km_cohort2_up_up_survival_guided.png", guided_dir)
save_km_plot(plot_km_survival_guided(merged_cohort1_down_guided, "KM - Cohort 1 | Down-Down Survival-Guided Score"),
             "km_cohort1_down_down_survival_guided.png", guided_dir)
save_km_plot(plot_km_survival_guided(merged_cohort2_down_guided, "KM - Cohort 2 | Down-Down Survival-Guided Score"),
             "km_cohort2_down_down_survival_guided.png", guided_dir)

# ---- Multivariable Cox: original panels + survival-guided scores ----------

message(">>> Multivariable Cox regression (original panels) ...")
merged_cohort1_up   <- merge_score_with_clinical(cohort1_tumor_mat, cohort1_clinical, panel_up_up,     "Cohort1 Up-Up")
merged_cohort2_up   <- merge_score_with_clinical(cohort2_tumor_mat, cohort2_clinical, panel_up_up,     "Cohort2 Up-Up")
merged_cohort1_down <- merge_score_with_clinical(cohort1_tumor_mat, cohort1_clinical, panel_down_down, "Cohort1 Down-Down")
merged_cohort2_down <- merge_score_with_clinical(cohort2_tumor_mat, cohort2_clinical, panel_down_down, "Cohort2 Down-Down")

adjusted_results <- do.call(rbind, list(
  run_adjusted_cox(merged_cohort1_up,   "Cohort1", "Up-Up original"),
  run_adjusted_cox(merged_cohort2_up,   "Cohort2", "Up-Up original"),
  run_adjusted_cox(merged_cohort1_down, "Cohort1", "Down-Down original"),
  run_adjusted_cox(merged_cohort2_down, "Cohort2", "Down-Down original")
))
adjusted_results$FDR <- p.adjust(adjusted_results$P_value, method = "BH")
write.csv(adjusted_results, file.path(adjusted_dir, "adjusted_cox_results_original_panels.csv"), row.names = FALSE)

message(">>> Multivariable Cox regression (survival-guided scores) ...")
adjusted_guided_results <- do.call(rbind, list(
  run_adjusted_cox(merged_cohort1_up_guided,   "Cohort1", "Up-Up survival-guided"),
  run_adjusted_cox(merged_cohort2_up_guided,   "Cohort2", "Up-Up survival-guided"),
  run_adjusted_cox(merged_cohort1_down_guided, "Cohort1", "Down-Down survival-guided"),
  run_adjusted_cox(merged_cohort2_down_guided, "Cohort2", "Down-Down survival-guided")
))
adjusted_guided_results$FDR <- p.adjust(adjusted_guided_results$P_value, method = "BH")
write.csv(adjusted_guided_results, file.path(adjusted_dir, "adjusted_cox_results_survival_guided.csv"), row.names = FALSE)

message("\nSurvival-guided refinement results saved in: ", SURVIVAL_DIR)
