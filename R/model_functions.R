# Core functions for the locked tuberculosis REO models.

validate_expression_matrix <- function(expression_matrix) {
  expression_matrix <- as.matrix(expression_matrix)
  storage.mode(expression_matrix) <- "double"
  if (is.null(rownames(expression_matrix)) || anyDuplicated(rownames(expression_matrix))) {
    stop("expression_matrix must have unique HGNC gene symbols as row names.")
  }
  if (is.null(colnames(expression_matrix))) {
    colnames(expression_matrix) <- paste0("sample_", seq_len(ncol(expression_matrix)))
  }
  expression_matrix
}

encode_reo_pairs <- function(expression_matrix, pair_manifest, task) {
  expression_matrix <- validate_expression_matrix(expression_matrix)
  task_pairs <- pair_manifest[pair_manifest$task == task, , drop = FALSE]
  task_pairs <- task_pairs[order(task_pairs$feature_order), , drop = FALSE]
  if (nrow(task_pairs) == 0L) stop("No pairs were found for task: ", task)

  required_genes <- unique(c(task_pairs$gene_a, task_pairs$gene_b))
  missing_genes <- setdiff(required_genes, rownames(expression_matrix))
  if (length(missing_genes)) {
    stop("Input matrix is missing required genes: ", paste(missing_genes, collapse = ", "))
  }

  encoded <- vapply(seq_len(nrow(task_pairs)), function(index) {
    gene_a <- expression_matrix[task_pairs$gene_a[index], ]
    gene_b <- expression_matrix[task_pairs$gene_b[index], ]
    ifelse(gene_a < gene_b, 1, ifelse(gene_a == gene_b, 0, -1))
  }, numeric(ncol(expression_matrix)))

  encoded <- as.matrix(encoded)
  rownames(encoded) <- colnames(expression_matrix)
  colnames(encoded) <- task_pairs$pair
  encoded
}

score_recovery_model <- function(expression_matrix, pair_manifest, coefficients) {
  encoded <- encode_reo_pairs(expression_matrix, pair_manifest, "recovery")
  coefficient_vector <- setNames(coefficients$coefficient, coefficients$term)
  intercept <- unname(coefficient_vector["(Intercept)"])
  beta <- coefficient_vector[colnames(encoded)]
  if (is.na(intercept) || anyNA(beta)) stop("Recovery coefficients do not match the encoded features.")
  probability_r_late <- stats::plogis(intercept + drop(encoded %*% beta))
  data.frame(
    sample_id = rownames(encoded),
    probability_R_late = probability_r_late,
    active_TB_like_score = 1 - probability_r_late,
    row.names = NULL,
    check.names = FALSE
  )
}

score_progression_model <- function(expression_matrix, pair_manifest, xgboost_model) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package 'xgboost' is required for progression-model scoring.")
  }
  encoded <- encode_reo_pairs(expression_matrix, pair_manifest, "progression")
  task_pairs <- pair_manifest[pair_manifest$task == "progression", , drop = FALSE]
  task_pairs <- task_pairs[order(task_pairs$feature_order), , drop = FALSE]
  colnames(encoded) <- paste0(task_pairs$gene_a, "__lt__", task_pairs$gene_b)
  probability_active_tb <- predict(xgboost_model, xgboost::xgb.DMatrix(encoded))
  data.frame(
    sample_id = rownames(encoded),
    probability_Active_TB = as.numeric(probability_active_tb),
    row.names = NULL,
    check.names = FALSE
  )
}
