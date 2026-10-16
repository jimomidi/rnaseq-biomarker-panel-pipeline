# ============================================================
# Shared helpers for survival analysis (R/09, R/10)
# ============================================================
#
# Requires config.R to have been sourced already, and the
# `survival`, `survminer`, and `ggplot2` packages to be loaded.

#' Load and standardize a GDC clinical.tsv file into overall-survival
#' (OS) time/event columns.
#'
#' Expects the standard GDC clinical export column names
#' (`cases.submitter_id`, `demographic.vital_status`,
#' `demographic.days_to_death`, `diagnoses.days_to_last_follow_up`).
#' OS time is days-to-death when the patient is deceased, otherwise
#' days-to-last-follow-up; converted to months (days / 30.44).
#'
#' @param clinical_dir Directory containing `clinical.tsv`.
#' @return Data frame with one row per case and added
#'   `OS_time` / `OS_event` / `OS_months` columns.
load_clinical <- function(clinical_dir) {
  path <- file.path(clinical_dir, "clinical.tsv")
  df <- read.delim(
    path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
    na.strings = c("'--", "--", "NA", "")
  )
  df <- df[!duplicated(df$cases.submitter_id), ]

  df$days_to_death    <- suppressWarnings(as.numeric(df$demographic.days_to_death))
  df$days_to_followup <- suppressWarnings(as.numeric(df$diagnoses.days_to_last_follow_up))
  df$vital_status     <- df$demographic.vital_status

  df$OS_time <- ifelse(
    !is.na(df$days_to_death) & df$days_to_death > 0,
    df$days_to_death,
    df$days_to_followup
  )
  df$OS_event  <- ifelse(tolower(df$vital_status) %in% c("dead", "deceased"), 1, 0)
  df$OS_months <- df$OS_time / 30.44

  df <- df[!is.na(df$OS_months) & df$OS_months > 0, ]

  message(sprintf(
    "  Clinical data loaded: %d patients | %d deaths | %d censored",
    nrow(df), sum(df$OS_event == 1), sum(df$OS_event == 0)
  ))
  df
}

#' Derive a case-level identifier from a longer sample barcode by
#' keeping its first 12 characters (the case-ID portion of a
#' standard TCGA-style barcode, e.g. "TCGA-XX-XXXX").
shorten_sample_id <- function(sample_ids) {
  substr(gsub("\\.", "-", sample_ids), 1, 12)
}

#' Composite panel score = mean per-gene z-score across `genes`.
#'
#' @return Data frame with `sample_id` and `score`, or `NULL` if none
#'   of `genes` are present in `expr_mat`.
compute_composite_score <- function(expr_mat, genes) {
  present_genes <- intersect(genes, rownames(expr_mat))
  if (length(present_genes) == 0) return(NULL)

  sub <- expr_mat[present_genes, , drop = FALSE]
  z <- t(apply(sub, 1, function(x) {
    s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(NA, length(x)))
    (x - mean(x, na.rm = TRUE)) / s
  }))
  score <- colMeans(z, na.rm = TRUE)

  data.frame(sample_id = names(score), score = as.numeric(score), stringsAsFactors = FALSE)
}

#' Merge a per-sample score with clinical data and split samples into
#' High/Low groups by the configured score quantiles, dropping the
#' middle band between them.
#'
#' @param expr_mat Expression matrix (genes x samples).
#' @param clinical_df Clinical data frame (see `load_clinical()`).
#' @param genes Panel genes used to compute the composite score.
#' @param dataset_name Label used in log messages.
#' @return Merged data frame with a `group` factor (`"Low"`/`"High"`),
#'   or `NULL` if no genes/samples matched.
merge_score_with_clinical <- function(expr_mat, clinical_df, genes, dataset_name) {
  score_df <- compute_composite_score(expr_mat, genes)
  if (is.null(score_df)) {
    message("  No panel genes found in ", dataset_name)
    return(NULL)
  }
  score_df$short_id <- shorten_sample_id(score_df$sample_id)

  merged <- merge(score_df, clinical_df, by.x = "short_id", by.y = "cases.submitter_id")
  if (nrow(merged) == 0) {
    message("  No matching samples in ", dataset_name)
    return(NULL)
  }

  merged <- merged[!is.na(merged$OS_months) & !is.na(merged$OS_event) & !is.na(merged$score), ]

  q_low  <- quantile(merged$score, SURVIVAL_LOWER_QUANTILE, na.rm = TRUE)
  q_high <- quantile(merged$score, SURVIVAL_UPPER_QUANTILE, na.rm = TRUE)
  merged <- merged[merged$score <= q_low | merged$score >= q_high, ]

  merged$group <- factor(ifelse(merged$score >= q_high, "High", "Low"), levels = c("Low", "High"))

  message(sprintf(
    "  %s - merged: %d | High: %d | Low: %d | Events: %d",
    dataset_name, nrow(merged), sum(merged$group == "High"),
    sum(merged$group == "Low"), sum(merged$OS_event == 1)
  ))
  merged
}

#' Log-rank test and univariate Cox proportional-hazards summary for a
#' High-vs-Low grouping.
#'
#' @return One-row data frame of summary statistics, or `NULL` if
#'   `merged_df` is empty.
compute_survival_stats <- function(merged_df, dataset_name, panel_name) {
  if (is.null(merged_df) || nrow(merged_df) == 0) return(NULL)

  surv_obj <- Surv(merged_df$OS_months, merged_df$OS_event)
  lr <- survdiff(surv_obj ~ group, data = merged_df)
  logrank_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)

  cox <- coxph(surv_obj ~ group, data = merged_df)
  cox_summary <- summary(cox)

  data.frame(
    Dataset = dataset_name,
    Panel = panel_name,
    N = nrow(merged_df),
    High_N = sum(merged_df$group == "High"),
    Low_N = sum(merged_df$group == "Low"),
    Events = sum(merged_df$OS_event == 1),
    Logrank_P = logrank_p,
    Cox_HR_High_vs_Low = cox_summary$coefficients[1, "exp(coef)"],
    Cox_HR_95CI_low = cox_summary$conf.int[1, "lower .95"],
    Cox_HR_95CI_high = cox_summary$conf.int[1, "upper .95"],
    Cox_P = cox_summary$coefficients[1, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
}

#' Format a p-value for display ("p < 0.001" or "p = 0.032").
format_pvalue <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("p < 0.001")
  paste0("p = ", sprintf("%.3f", p))
}

#' Kaplan-Meier plot for a High/Low grouping, annotated with the
#' log-rank p-value and Cox hazard ratio.
plot_km_with_hr <- function(merged_df, title, stats_row, legend_title = "Composite score") {
  if (is.null(merged_df) || nrow(merged_df) == 0) return(NULL)

  fit <- survfit(Surv(OS_months, OS_event) ~ group, data = merged_df)

  hr_text <- paste0(
    "HR = ", sprintf("%.2f", stats_row$Cox_HR_High_vs_Low),
    " (95% CI: ", sprintf("%.2f", stats_row$Cox_HR_95CI_low),
    "-", sprintf("%.2f", stats_row$Cox_HR_95CI_high), ")"
  )
  info_text <- paste0(
    format_pvalue(stats_row$Logrank_P), "\n", hr_text,
    "\nHigh n = ", stats_row$High_N, " | Low n = ", stats_row$Low_N
  )

  p_obj <- ggsurvplot(
    fit, data = merged_df,
    palette = c("#2E86AB", "#E84855"),
    conf.int = TRUE, conf.int.alpha = 0.12,
    pval = FALSE, risk.table = FALSE,
    xlab = "Time (months)", ylab = "Overall Survival Probability",
    title = title, legend.title = legend_title, legend.labs = c("Low", "High"),
    ggtheme = theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.82, 0.88),
        legend.background = element_rect(fill = "white", color = "gray80")
      ),
    fontsize = 3.8
  )

  p_obj$plot <- p_obj$plot +
    annotate("text", x = Inf, y = Inf, label = info_text,
             hjust = 1.03, vjust = 1.05, size = 3.3, color = "black", fontface = "bold")
  p_obj
}

#' Save a `ggsurvplot` object to a square PNG.
save_km_plot <- function(p_obj, filename, output_dir) {
  if (is.null(p_obj)) return(invisible(NULL))
  png(file.path(output_dir, filename), width = 2000, height = 2000, res = 250, bg = "white")
  print(p_obj)
  dev.off()
  message("  [saved] ", filename)
}
