# Within-cohort sensitivity analysis for inferred immune-cell composition.

adjust_score_for_cell_composition <- function(data, score_column, fraction_columns) {
  required <- c(score_column, fraction_columns)
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns)) stop("Missing columns: ", paste(missing_columns, collapse = ", "))

  complete_index <- which(stats::complete.cases(data[, required, drop = FALSE]))
  analysis_data <- data[complete_index, required, drop = FALSE]
  if (nrow(analysis_data) <= length(fraction_columns) + 1L) {
    stop("Too few complete observations for the requested adjustment model.")
  }

  formula <- stats::reformulate(fraction_columns, response = score_column)
  fit <- stats::lm(formula, data = analysis_data)
  output <- data
  output$cell_adjusted_residual <- NA_real_
  output$cell_adjusted_residual[complete_index] <- stats::residuals(fit)
  attr(output, "adjustment_model") <- fit
  output
}

# The fraction coefficients are estimated within the analyzed cohort. Therefore,
# residual-score ROC curves are sensitivity analyses and should not be described
# as independently locked external validations.
