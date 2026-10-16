# RNA-seq Biomarker Panel Pipeline

An R pipeline for discovering and validating small gene expression panels
that are reproducible across two independent RNA-seq cohorts. Starting from
raw case (tumor) vs. control (normal) RNA-seq counts, it performs
differential expression in each cohort, intersects the two DEG signatures,
ranks candidate genes by single-gene discriminative power, and assembles
compact "Up-Up" / "Down-Down" panels that are then validated with ROC
analysis, expression boxplots, and survival analysis.

The pipeline was developed against public GDC (Genomic Data Commons)
RNA-seq data and an external recount3-hosted control cohort, but the
scripts are written generically: which GDC projects, subtype filters, and
control cohort to use are all configuration values in `config.R`, not
hard-coded assumptions.

## Workflow

| Step | Script | Description |
|------|--------|-------------|
| 1 | `R/01_download_prepare_cohort1.R` | Download cohort 1 tumor (subtype-filtered) and normal RNA-seq counts from the GDC. |
| 2 | `R/02_download_prepare_cohort2.R` | Download cohort 2 tumor RNA-seq counts from the GDC, retrieve an external control expression set (recount3), merge, TMM-normalize, and run a batch-effect diagnostic (PCA + PERMANOVA). |
| 3 | `R/03_differential_expression_cohort1.R` | limma-voom differential expression, cohort 1 (tumor vs. normal). |
| 4 | `R/04_differential_expression_cohort2.R` | limma differential expression on the batch-corrected log-CPM matrix, cohort 2 (tumor vs. normal). |
| 5 | `R/05_shared_deg_overlap.R` | Intersect the two DEG signatures into four direction patterns (Up_Up, Down_Down, Up_Down, Down_Up). |
| 6 | `R/06_gene_panel_selection.R` | Rank Up_Up/Down_Down candidate genes by single-gene AUC in each cohort; select genes ranked highly in both. |
| 7 | `R/07_panel_roc_evaluation.R` | ROC/AUC of each panel's composite (mean z-score) score, per cohort. |
| 8 | `R/08_panel_expression_boxplots.R` | Tumor vs. normal expression boxplots for each panel gene, per cohort. |
| 9 | `R/09_survival_analysis.R` | Kaplan-Meier survival by panel composite score and by individual gene. |
| 10 | `R/10_survival_guided_panel_refinement.R` | Single-gene Cox screening, a survival-weighted risk score, and multivariable Cox regression adjusted for clinical covariates. |

`R/utils_expression.R` and `R/utils_survival.R` hold helper functions shared
across the numbered scripts; `config.R` centralizes every path and analysis
threshold.

Run the scripts in order (1 → 10); each one consumes files written by an
earlier step. Steps 1-2 (data download) only need to be run once per cohort.

## Requirements

* R >= 4.2
* CRAN packages: `dplyr`, `ggplot2`, `ggrepel`, `survival`, `survminer`,
  `pROC`, `vegan`
* Bioconductor packages: `TCGAbiolinks`, `SummarizedExperiment`, `recount3`,
  `edgeR`, `limma`

Install everything with:

```r
Rscript install_packages.R
```

## Installation

```bash
git clone <this-repository-url>
cd rnaseq-biomarker-panel-pipeline
Rscript install_packages.R
```

## Configuration

Before running the pipeline, edit `config.R`:

* `GDC_PROJECT_COHORT1`, `GDC_PROJECT_COHORT2` — GDC project identifiers
  for your two cohorts (e.g. TCGA project codes).
* `COHORT1_SUBTYPE_KEYWORDS` — keywords used to restrict cohort 1 to a
  specific tumor subtype, matched case-insensitively against the clinical
  site/diagnosis column.
* `COHORT2_HISTOLOGY_KEYWORD` — keyword used to restrict cohort 2 to a
  specific histological subtype.
* `RECOUNT3_CONTROL_PROJECT`, `RECOUNT3_CONTROL_PROJECT_HOME` — the
  recount3 project used as an external control/normal reference for
  cohort 2 (used because that cohort's own GDC project lacks enough
  matched normal samples). `RECOUNT3_CONTROL_PROJECT_HOME` must be one
  of recount3's own source identifiers (see the recount3 documentation).
* Analysis thresholds (low-expression filter, DEG FDR/fold-change cutoffs,
  panel size, survival split quantiles, etc.) — the defaults reproduce the
  original analysis and generally do not need to change.

## Expected inputs

* **Expression data**: downloaded automatically by steps 1-2 via
  `TCGAbiolinks` (GDC) and `recount3`. No manual download is required.
* **Clinical/survival data** (steps 9-10 only): a `clinical.tsv` file per
  cohort, in the standard tab-separated format exported by the GDC portal's
  clinical data download, placed at:
  * `data/clinical/cohort1/clinical.tsv`
  * `data/clinical/cohort2/clinical.tsv`

  At minimum, the file must contain `cases.submitter_id`,
  `demographic.vital_status`, `demographic.days_to_death`, and
  `diagnoses.days_to_last_follow_up`. The optional covariate-adjusted
  analysis in step 10 additionally uses age, gender, and AJCC stage columns
  when present.

## Usage

```bash
Rscript R/01_download_prepare_cohort1.R
Rscript R/02_download_prepare_cohort2.R
Rscript R/03_differential_expression_cohort1.R
Rscript R/04_differential_expression_cohort2.R
Rscript R/05_shared_deg_overlap.R
Rscript R/06_gene_panel_selection.R
Rscript R/07_panel_roc_evaluation.R
Rscript R/08_panel_expression_boxplots.R
Rscript R/09_survival_analysis.R
Rscript R/10_survival_guided_panel_refinement.R
```

Run every script from the repository root so the relative paths defined in
`config.R` resolve correctly.

## Outputs

All generated data and results are written outside version control (see
`.gitignore`):

```
data/
  raw/cohort1/           # GDC expression matrices and clinical data, cohort 1
  raw/cohort2/           # GDC expression matrices and clinical data, cohort 2
  clinical/cohort1/      # clinical.tsv for cohort 1 (user-supplied)
  clinical/cohort2/      # clinical.tsv for cohort 2 (user-supplied)
  processed/             # merged, TMM-normalized cohort 2 matrix + batch diagnostic plot
                          # (batch correction itself is applied per-analysis; see below)

results/
  differential_expression/   # per-cohort DEG tables and volcano plots
  shared_deg/                 # cross-cohort DEG overlap tables and volcano plots
  gene_panels/                 # AUC ranking tables and final gene panels
  roc_curves/                  # panel ROC plots
  expression_boxplots/         # panel expression boxplots
  survival_analysis/           # KM plots, Cox statistics, survival-guided refinement
```

## Reproducibility

* A fixed random seed (`RANDOM_SEED <- 42` in `config.R`) is set wherever
  randomness is used.
* All file paths are relative to the repository root; no absolute or
  machine-specific paths are used.
* All thresholds (low-expression filter, DEG significance cutoffs, shared-DEG
  cutoffs, panel size, survival group quantiles) are defined once in
  `config.R` and consumed by every script that needs them.

## Notes on methodology

* **Batch correction for cohort 2 is applied once, consistently, wherever
  it is needed.** Cohort 2's tumor and control samples come from different
  studies, so batch and biological group are confounded (every tumor
  sample is one batch, every control sample is the other). `R/02` merges
  and TMM-normalizes the two sources and runs a diagnostic (PCA +
  PERMANOVA) confirming the batch effect is present, but does not correct
  it — the actual correction (`limma::removeBatchEffect()`, with the
  tumor/normal group protected as a covariate) lives in a single place,
  `load_cohort2_expression()` in `R/utils_expression.R`, and is used by
  every analysis that compares tumor vs. normal for cohort 2: differential
  expression (step 4), gene panel ranking (step 6), ROC evaluation
  (step 7), and boxplots (step 8). Tumor-only analyses (the survival
  scripts, steps 9-10) don't need this correction, since all of cohort 2's
  tumor samples come from a single batch. Cohort 1 has native matched
  normal samples from the same study, so no batch correction is needed
  there at all.
* **Cohort 1 vs. cohort 2 preprocessing still differ in one respect.**
  Cohort 1's DEG analysis (step 3) uses raw counts with limma-voom, while
  cohort 2's uses limma directly on log-CPM (already normalized, so voom
  is not applicable). The survival scripts inherit this: cohort 1's
  composite/single-gene scores are computed on raw counts, cohort 2's on
  TMM log-CPM. If you adapt the pipeline to your own data where both
  cohorts have native matched controls, consider normalizing both the
  same way before scoring.
* **Two survival group-splitting conventions are used deliberately.** The
  panel-level survival analysis (step 9) splits samples into High/Low
  groups using the `SURVIVAL_LOWER_QUANTILE`/`SURVIVAL_UPPER_QUANTILE`
  quantiles (40%/60% by default), excluding the middle band. The
  survival-guided risk score (step 10) instead uses a plain median split.
  Both are intentional design choices from the original analysis, kept as
  documented in the step 10 script header.
* **Differential expression uses one consistent significance threshold.**
  Both cohorts define "differentially expressed" as FDR < `DEG_FDR_CUTOFF`
  and |log2 fold change| >= `DEG_LOGFC_CUTOFF` (0.05 and 1.5 by default),
  applied identically before the two DEG lists are intersected in step 5.
* **Step 10 goes beyond panel validation.** `R/10_survival_guided_panel_refinement.R`
  screens panel genes individually for a survival association and builds
  a separate, HR-weighted risk score — an exploratory extension, not a
  re-statement of the panel-level results in step 9. It's clearly marked
  as optional/secondary in its own header; skip it if you only need
  diagnostic panel discovery and validation (steps 1-9).

## Customization

* Swap in your own GDC projects, subtype keywords, and control cohort via
  `config.R` — no changes to the analysis scripts are required.
* Adjust `PANEL_SIZE`, the low-expression filter, or the DEG/shared-DEG
  significance cutoffs in `config.R` to explore different panel sizes or
  stringency.
* The survival scripts (9-10) are optional and only run if clinical data is
  available; the diagnostic panel discovery and validation (steps 1-8)
  work independently of survival data.
