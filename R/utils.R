#' Check and prepare data for MAR bounds estimation
#'
#' @param data A data.frame.
#' @param Y Character name of outcome column (numeric, may be NA when missing).
#' @param A Character name of treatment column (0/1).
#' @param C Character name of missingness column (1 = missing, 0 = observed).
#' @param X Character vector of covariate column names.
#' @return List with validated vectors/matrices: Y_vec, A_vec, C_vec, X_mat, n, valid (logical for complete treatment/missingness).
#' @keywords internal
prepare_data <- function(data, Y, A, C, X) {
  if (!is.data.frame(data)) stop("data must be a data.frame")
  n <- nrow(data)
  Y_vec <- data[[Y]]
  A_vec <- data[[A]]
  C_vec <- data[[C]]
  if (is.null(Y_vec) || is.null(A_vec) || is.null(C_vec)) stop("Y, A, C must be column names in data")
  if (!all(A_vec %in% c(0, 1))) stop("A must be 0/1")
  if (!all(C_vec %in% c(0, 1))) stop("C must be 0/1")
  X_mat <- as.matrix(data[, X, drop = FALSE])
  if (ncol(X_mat) != length(X)) stop("All X columns must exist")
  list(Y_vec = Y_vec, A_vec = A_vec, C_vec = C_vec, X_mat = X_mat, n = n, X_names = X)
}

#' Clip propensity and probability estimates away from 0 and 1
#'
#' @param p Numeric vector of probabilities.
#' @param eps Small positive number (default 1e-6).
#' @return Clipped vector.
clip_probs <- function(p, eps = 1e-6) {
  pmax(eps, pmin(1 - eps, p))
}

#' Expit (inverse logit)
expit <- function(x) {
  1 / (1 + exp(-x))
}

#' Default value when argument is NULL
`%||%` <- function(x, y) if (is.null(x)) y else x
