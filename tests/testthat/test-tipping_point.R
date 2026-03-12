test_that("tipping_point returns grid and optional tipping", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_tiny_data(n = 60)
  fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                    estimand = "ate", assumption = "monotonicity_pos",
                    delta_0u = 0.8, delta_1u = 0.8, V = 2, seed = 1)
  tp <- tipping_point(fit, bound_type = "upper", assumption = "monotonicity_pos",
                      param_grid = list(delta_0u = c(0.3, 0.6, 1), delta_1u = c(0.3, 0.6, 1)))
  expect_true(all(c("grid", "tipping", "target", "bound_type", "assumption") %in% names(tp)))
  expect_equal(tp$bound_type, "upper")
  expect_true("bound_value" %in% names(tp$grid))
  expect_equal(nrow(tp$grid), 9)
  expect_true(all(is.finite(tp$grid$bound_value)))
})

test_that("tipping_point with default grid runs", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_tiny_data(n = 60)
  fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                    estimand = "ate", assumption = "bounded_delta", V = 2, seed = 1)
  tp <- tipping_point(fit, bound_type = "lower", assumption = "bounded_delta", param_grid = NULL)
  expect_true(nrow(tp$grid) >= 1)
  expect_equal(tp$target, 0)
})

test_that("tipping_point errors without phi and nuisance", {
  expect_error(
    tipping_point(list(phi = NULL, nuisance = NULL), bound_type = "upper", assumption = "monotonicity_pos"),
    "phi and nuisance"
  )
})

test_that("tipping_point with inference = pointwise returns pointwise SE and CIs", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_tiny_data(n = 60)
  fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                    estimand = "ate", assumption = "monotonicity_pos",
                    delta_0u = 0.8, delta_1u = 0.8, V = 2, seed = 1)
  tp <- tipping_point(fit, bound_type = "upper", assumption = "monotonicity_pos",
                      param_grid = list(delta_0u = c(0.5, 1), delta_1u = c(0.5, 1)),
                      inference = "pointwise", alpha = 0.05)
  expect_true(all(c("se_pointwise", "ci_lower_pointwise", "ci_upper_pointwise") %in% names(tp$grid)))
  expect_true(all(tp$grid$ci_lower_pointwise <= tp$grid$bound_value))
  expect_true(all(tp$grid$bound_value <= tp$grid$ci_upper_pointwise))
  expect_true(all(tp$grid$se_pointwise >= 0))
})

test_that("tipping_point with inference = both returns pointwise and uniform CIs", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_tiny_data(n = 60)
  fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                    estimand = "ate", assumption = "monotonicity_pos",
                    delta_0u = 0.8, delta_1u = 0.8, V = 2, seed = 1)
  tp <- tipping_point(fit, bound_type = "upper", assumption = "monotonicity_pos",
                      param_grid = list(delta_0u = c(0.5, 1), delta_1u = c(0.5, 1)),
                      inference = "both", alpha = 0.05, B = 50)
  expect_true(all(c("se_pointwise", "ci_lower_pointwise", "ci_upper_pointwise",
                    "ci_lower_uniform", "ci_upper_uniform", "se_uniform") %in% names(tp$grid)))
  expect_true(all(tp$grid$ci_lower_uniform <= tp$grid$ci_upper_uniform))
  expect_true("critical_value" %in% names(tp))
  expect_true(tp$critical_value >= 0)
  # C replaces z: uniform band should be wider than pointwise (C > z_{1-alpha/2})
  expect_true(all(tp$grid$se_uniform >= tp$grid$se_pointwise))
  z <- qnorm(1 - 0.05/2)
  expect_true(tp$critical_value >= z)
})
