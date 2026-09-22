# ROC/AUC evaluation for a locked model score.

evaluate_locked_score <- function(labels, scores, negative_class, positive_class, bootstrap_replicates = 2000L, seed = 20260922L) {
  if (!requireNamespace("pROC", quietly = TRUE)) stop("Package 'pROC' is required.")
  keep <- stats::complete.cases(labels, scores)
  labels <- labels[keep]
  scores <- scores[keep]
  if (!all(c(negative_class, positive_class) %in% unique(labels))) {
    stop("Both requested outcome classes must be present.")
  }

  roc_object <- pROC::roc(
    response = labels,
    predictor = scores,
    levels = c(negative_class, positive_class),
    direction = "<",
    quiet = TRUE
  )
  set.seed(seed)
  confidence_interval <- pROC::ci.auc(
    roc_object,
    method = "bootstrap",
    boot.n = bootstrap_replicates,
    stratified = TRUE
  )
  data.frame(
    n = length(labels),
    negative_n = sum(labels == negative_class),
    positive_n = sum(labels == positive_class),
    auc = as.numeric(pROC::auc(roc_object)),
    ci_lower = as.numeric(confidence_interval[1]),
    ci_upper = as.numeric(confidence_interval[3]),
    check.names = FALSE
  )
}

