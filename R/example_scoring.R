source("R/model_functions.R")

# Replace this path with a CSV containing genes in rows and samples in columns.
expression_matrix <- read.csv(
  "expression_matrix.csv",
  row.names = 1,
  check.names = FALSE
)

pair_manifest <- read.csv("models/selected_REO_pairs.csv", check.names = FALSE)
coefficients <- read.csv("models/recovery_logistic_coefficients.csv", check.names = FALSE)
model_bundle <- readRDS("models/frozen_model_bundle.rds")

recovery_scores <- score_recovery_model(expression_matrix, pair_manifest, coefficients)
progression_scores <- score_progression_model(
  expression_matrix,
  pair_manifest,
  model_bundle$progression_xgboost_model
)

write.csv(recovery_scores, "recovery_scores.csv", row.names = FALSE)
write.csv(progression_scores, "progression_scores.csv", row.names = FALSE)

