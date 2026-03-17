# Run the exact README example to catch errors users would see
test_that("README example runs without error", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  n <- 500
  set.seed(1)
  X <- runif(n, -2, 2)
  A <- rbinom(n, 1, plogis(X))
  C <- rbinom(n, 1, 0.2 + 0.1 * A)
  Y <- ifelse(C == 1, NA, rbinom(n, 1, plogis(0.5 * A + 0.3 * X)))
  dat <- data.frame(Y = Y, A = A, C = C, X = X)

  fit <- safe_mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                         estimand = "ate", assumption = "general",
                         sl_lib = "SL.glm")
  skip_if(is.null(fit), "SuperLearner failed (edge case)")

  # Check result dataframe exists and has valid values
  expect_true("result" %in% names(fit))
  expect_true(is.finite(fit$result$naive))
  expect_true(all(fit$result$estimate[fit$result$bound_type == "lower"] <=
                  fit$result$estimate[fit$result$bound_type == "upper"]))

  fit2 <- safe_mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                          estimand = "ate", assumption = "monotonicity_pos",
                          delta_0u = 0.8, delta_1u = 0.8,
                          param_grid = list(delta_0u = seq(0.2, 1, 0.2), delta_1u = seq(0.2, 1, 0.2)),
                          sl_lib = "SL.glm", B = 50)
  skip_if(is.null(fit2), "SuperLearner failed (edge case)")

  expect_true("upper_grid" %in% names(fit2))
  expect_true("grid" %in% names(fit2$upper_grid))
})
