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
  fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                    estimand = "ate", assumption = "general")
  expect_true(is.finite(fit$naive))
  expect_true(fit$lower <= fit$upper)
  fit2 <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                     estimand = "ate", assumption = "monotonicity_pos",
                     delta_0u = 0.8, delta_1u = 0.8)
  tp <- tipping_point(fit2, bound_type = "upper", assumption = "monotonicity_pos")
  expect_true("grid" %in% names(tp))
  mb <- multiplier_bootstrap_grid(fit2,
                                  param_grid = list(delta_0u = seq(0.2, 1, 0.2), delta_1u = seq(0.2, 1, 0.2)),
                                  bound_spec = "upper", assumption = "monotonicity_pos", B = 50)
  expect_true("grid" %in% names(mb))
})
