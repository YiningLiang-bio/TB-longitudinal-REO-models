# Longitudinal blood-transcriptomic REO models for tuberculosis

This repository contains the reproducibility materials for a manuscript describing two locked relative-expression-ordering (REO) models of tuberculosis progression and recovery.

- **Progression model:** P1 versus Active TB; eight-pair XGBoost classifier.
- **Recovery model:** Active TB versus R late; eight-pair binomial logistic-regression classifier.

The repository follows the structure commonly used by transcriptomic biomarker studies: analysis functions are separated from frozen model assets, manuscript source tables and computational-environment records.

## Repository contents

```text
R/
  model_functions.R                  Core REO encoding and locked-model scoring
  evaluate_roc.R                     Cohort-level ROC/AUC evaluation
  cell_composition_adjustment.R      Score residualization using inferred cell fractions
  example_scoring.R                  Minimal worked example using a gene-by-sample matrix
models/
  selected_REO_pairs.csv             Ordered features for both locked models
  recovery_logistic_coefficients.csv Fitted Active TB-R late coefficients
  progression_xgboost_specification.csv
  frozen_model_bundle.rds            Frozen progression XGBoost model plus public model metadata
results/
  external_validation_auc.csv        Archived and replayed external-validation results
source_tables/
  20260922_Table1_v1.xlsx
  20260922_Table2_v1.xlsx
  20260922_Supplementary_Table_v2.xlsx
environment/
  sessionInfo.txt
```

## Input format

The scoring functions expect a numeric gene-expression matrix with HGNC gene symbols as row names and samples as columns. REO encoding is performed within each sample, so no cross-sample normalization is applied by the scoring functions. The input matrix must contain all genes used by the selected model.

For a pair written as `GENEA-GENEB`, the encoded value is:

- `+1` when `GENEA < GENEB`;
- `0` when the two values are equal;
- `-1` when `GENEA > GENEB`.

## Quick start

```r
source("R/model_functions.R")

expression_matrix <- read.csv(
  "path/to/expression_matrix.csv",
  row.names = 1,
  check.names = FALSE
)

pair_manifest <- read.csv("models/selected_REO_pairs.csv", check.names = FALSE)
coefficients <- read.csv("models/recovery_logistic_coefficients.csv", check.names = FALSE)
model_bundle <- readRDS("models/frozen_model_bundle.rds")

recovery_scores <- score_recovery_model(
  expression_matrix,
  pair_manifest,
  coefficients
)

progression_scores <- score_progression_model(
  expression_matrix,
  pair_manifest,
  model_bundle$progression_xgboost_model
)
```

The recovery function returns both `probability_R_late` and `active_TB_like_score = 1 - probability_R_late`. The progression function returns `probability_Active_TB`.

## Data availability

The study uses public transcriptomic datasets from NCBI GEO and TCGA. Accession-level roles, populations and sample counts are listed in `source_tables/20260922_Table1_v1.xlsx` and Supplementary Table S1. Raw or processed expression matrices are not redistributed here; they should be obtained from their original repositories.

## Reproducibility scope

This release provides the locked feature sets, model specifications, fitted recovery-model coefficients, frozen progression model, core scoring and evaluation code, source tables and session information. Cohort-specific preprocessing follows the procedures described in the manuscript Methods and Supplementary Table S1.

## Citation

Please cite the associated manuscript. Full citation information will be added after publication.

## License

Code is released under the MIT License. The manuscript tables remain subject to the terms of the associated article.

