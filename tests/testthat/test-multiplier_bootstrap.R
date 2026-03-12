test_that("multiplier_bootstrap_grid returns grid with estimates and CIs", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_tiny_data(n = 60)
  fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                    estimand = "ate", assumption = "monotonicity_pos",
                    delta_0u = 0.8, delta_1u = 0.8, V = 2, seed = 1)
  pg <- list(delta_0u = c(0.5, 1), delta_1u = c(0.5, 1))
  mb <- multiplier_bootstrap_grid(fit, param_grid = pg, bound_spec = "upper",
                                  assumption = "monotonicity_pos", B = 30, alpha = 0.05)
  expect_named(mb, c("grid", "boot_replicates", "critical_value", "alpha"))
  expect_equal(nrow(mb$grid), 4)
  expect_true(all(c("estimate", "se", "ci_lower", "ci_upper", "simultaneous_lower", "simultaneous_upper") %in% names(mb$grid)))
  expect_equal(dim(mb$boot_replicates), c(30, 4))
  expect_true(mb$critical_value >= 0)
  expect_equal(mb$alpha, 0.05)
  expect_true(all(mb$grid$ci_lower <= mb$grid$estimate))
  expect_true(all(mb$grid$estimate <= mb$grid$ci_upper))
  expect_true(all(mb$grid$simultaneous_lower <= mb$grid$simultaneous_upper))
})

test_that("multiplier_bootstrap_grid with point_ate assumption", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_tiny_data(n = 50)
  fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                    estimand = "ate", assumption = "point_ate",
                    delta_0 = 0.5, delta_1 = 0.5, tau = 2, V = 2, seed = 1)
  pg <- list(tau = c(1.5, 2), delta_0 = 0.5, delta_1 = 0.5)
  mb <- multiplier_bootstrap_grid(fit, param_grid = pg, bound_spec = "point",
                                  assumption = "point_ate", B = 20)
  expect_equal(nrow(mb$grid), 2)
  expect_true(all(is.finite(mb$grid$estimate)))
})

test_that("multiplier_bootstrap_grid errors without phi", {
  expect_error(
    multiplier_bootstrap_grid(list(phi = NULL), param_grid = list(delta_0u = 0.5), bound_spec = "upper", assumption = "monotonicity_pos"),
    "must contain phi"
  )
})
