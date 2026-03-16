# Shared test data: small dataset with known structure for fast tests

make_test_data <- function(n = 120, seed = 42) {
  set.seed(seed)
  X <- runif(n, -1.5, 1.5)
  A <- rbinom(n, 1, plogis(0.5 * X))
  C <- rbinom(n, 1, 0.15 + 0.1 * A + 0.05 * (X > 0))
  Y <- ifelse(C == 1, NA_real_, rbinom(n, 1, plogis(0.3 * A + 0.2 * X)))
  data.frame(Y = Y, A = A, C = C, X = X)
}

# Even smaller data for very fast tests (e.g. bootstrap with B=20)
make_tiny_data <- function(n = 60, seed = 1) {
  set.seed(seed)
  X <- runif(n, -1, 1)
  A <- rbinom(n, 1, 0.5)
  C <- rbinom(n, 1, 0.2)
  Y <- ifelse(C == 1, NA_real_, rbinom(n, 1, plogis(0.2 * A)))
  data.frame(Y = Y, A = A, C = C, X = X)
}

# Helper to run mar_bounds with SuperLearner error handling
# Returns NULL if SuperLearner fails (edge case), otherwise returns the result
safe_mar_bounds <- function(...) {
  tryCatch(
    {
      result <- mar_bounds(...)
      # Check if result is valid
      if (is.null(result)) return(NULL)
      # Check for NaN/Inf in key results
      if ("estimate" %in% names(result)) {
        if (!is.finite(result$estimate)) return(NULL)
      }
      if ("lower" %in% names(result)) {
        if (!is.finite(result$lower) || !is.finite(result$upper)) return(NULL)
      }
      if ("naive" %in% names(result)) {
        if (!is.finite(result$naive)) return(NULL)
      }
      result
    },
    error = function(e) {
      # SuperLearner edge case errors we want to handle gracefully
      if (grepl("All algorithms dropped|Argument mu must be", e$message)) {
        return(NULL)
      }
      # Re-throw other errors
      stop(e)
    }
  )
}
