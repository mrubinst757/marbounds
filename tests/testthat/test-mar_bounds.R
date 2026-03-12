test_that("mar_bounds general ATE returns expected structure", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 100)
  fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                    estimand = "ate", assumption = "general",
                    V = 2, sl_lib_prop = "SL.glm", sl_lib_miss = "SL.glm", sl_lib_outcome = "SL.glm",
                    seed = 1)
  expect_named(fit, c("naive", "lower", "upper", "se_lower", "se_upper", "nuisance", "phi", "estimand", "assumption"))
  expect_equal(fit$estimand, "ate")
  expect_equal(fit$assumption, "general")
  expect_true(fit$lower <= fit$upper)
  expect_true(all(is.finite(c(fit$naive, fit$lower, fit$upper))))
  expect_equal(nrow(fit$phi), 100)
  expect_equal(ncol(fit$phi), 6)
})

test_that("mar_bounds point_ate returns estimate and se", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 100)
  fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                    estimand = "ate", assumption = "point_ate",
                    delta_0 = 0.5, delta_1 = 0.5, tau = 2,
                    V = 2, sl_lib_prop = "SL.glm", sl_lib_miss = "SL.glm", sl_lib_outcome = "SL.glm",
                    seed = 1)
  expect_true("estimate" %in% names(fit))
  expect_true("se" %in% names(fit))
  expect_false("lower" %in% names(fit))
  expect_true(is.finite(fit$estimate))
})

test_that("mar_bounds errors on missing required columns", {
  dat <- make_test_data(n = 30)
  expect_error(mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "Missing"), "exist|undefined columns")
})

test_that("mar_bounds point_ate errors without delta_0, delta_1, tau", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 50)
  expect_error(
    mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X", estimand = "ate", assumption = "point_ate", V = 2),
    "delta_0|delta_1|tau"
  )
})

test_that("psi_naive returns EIF-based estimate, var, se, and eif", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 80)
  out <- psi_naive(dat, Y = "Y", A = "A", C = "C", X = "X")
  expect_true(all(c("estimate", "var", "se", "eif") %in% names(out)))
  expect_true(is.finite(out$estimate))
  expect_true(out$se >= 0)
  expect_equal(out$se^2, out$var, tolerance = 1e-10)
  expect_length(out$eif, 80)
})

test_that("psi_naive with precomputed nuisance matches mar_bounds naive", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 60)
  fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X", estimand = "ate", assumption = "general", V = 2, seed = 1)
  out <- psi_naive(dat, Y = "Y", A = "A", C = "C", X = "X", nuis = fit$nuisance)
  expect_equal(unname(out$estimate), unname(fit$naive), tolerance = 1e-6)
})

test_that("mar_bounds param_grid returns lower_grid and upper_grid for Psi_1 bounded_delta", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 80)
  fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                    estimand = "psi1", assumption = "bounded_delta",
                    param_grid = list(delta_0u = c(0.3, 0.6), delta_1u = c(0.3, 0.6)),
                    B = 50, seed = 1)
  expect_true("lower_grid" %in% names(fit))
  expect_true("upper_grid" %in% names(fit))
  g <- fit$lower_grid$grid
  expect_true(all(c("estimate", "se_pointwise", "se_uniform", "ci_lower_pointwise", "ci_upper_pointwise",
                    "ci_lower_uniform", "ci_upper_uniform") %in% names(g)))
  expect_equal(nrow(g), 4L)
  expect_true(is.finite(fit$lower_grid$critical_value) && fit$lower_grid$critical_value > 0)
})

test_that("mar_bounds bounded_risk_unbounded_tau matches original coefficient bounds", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 80)
  fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                    estimand = "ate", assumption = "bounded_risk_unbounded_tau",
                    delta_0u = 0.3, delta_1u = 0.3, tau_0 = 2, tau_1 = 2,
                    V = 2, seed = 1)
  expected <- marbounds:::compute_bounds(fit$phi, estimand = "ate", assumption = "bounded_risk_unbounded_tau",
                                         delta_0u = 0.3, delta_1u = 0.3, tau_0 = 2, tau_1 = 2)
  expect_equal(fit$lower, expected$lower)
  expect_equal(fit$upper, expected$upper)
})

test_that("mar_bounds param_grid accepts underscored and non-underscored delta naming", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 80)
  fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                    estimand = "psi1", assumption = "bounded_delta",
                    param_grid = list(delta0u = c(0.3, 0.6), delta1u = c(0.3, 0.6)),
                    B = 50, seed = 1)
  expect_true("lower_grid" %in% names(fit))
  expect_true("upper_grid" %in% names(fit))
  expect_equal(nrow(fit$lower_grid$grid), 4L)
})

test_that("mar_bounds param_grid works for general assumption with delta parameters", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 80)
  fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                    estimand = "ate", assumption = "general",
                    param_grid = list(delta_0u = c(0.2, 0.4), delta_1u = c(0.2, 0.4)),
                    B = 50, seed = 1)
  expect_true("lower_grid" %in% names(fit))
  expect_true("upper_grid" %in% names(fit))
  expect_equal(nrow(fit$lower_grid$grid), 4L)
})

test_that("mar_bounds uses a shared multiplier critical value for both bounds", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 80)
  fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                    estimand = "psi1", assumption = "bounded_delta",
                    param_grid = list(delta_0u = c(0.3, 0.6), delta_1u = c(0.3, 0.6)),
                    grid_bounds = "both",
                    B = 50, seed = 1)
  expect_equal(fit$lower_grid$critical_value, fit$upper_grid$critical_value)
})

test_that("mar_bounds drops constant param columns from grid output", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 80)
  fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                    estimand = "psi1", assumption = "bounded_delta",
                    param_grid = list(delta_0u = c(0.3, 0.6), delta_1u = c(0.3, 0.6), tau = 2),
                    B = 50, seed = 1)
  expect_false("tau" %in% names(fit$lower_grid$grid))
})

test_that("mar_bounds param_grid returns estimate_grid for Psi_2 point_psi2", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 80)
  fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                    estimand = "psi2", assumption = "point_psi2", delta_0 = 0.5,
                    param_grid = list(delta_0 = c(0.2, 0.5, 0.8)),
                    B = 50, seed = 1)
  expect_true("estimate_grid" %in% names(fit))
  g <- fit$estimate_grid$grid
  expect_true(all(c("estimate", "se_pointwise", "se_uniform", "ci_lower_uniform", "ci_upper_uniform") %in% names(g)))
  expect_equal(nrow(g), 3L)
})

test_that("mar_bounds warns when delta_0/delta_1 used with bounded_risk assumptions", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 80)

  # Should warn for bounded_risk
  expect_warning(
    mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
               estimand = "ate", assumption = "bounded_risk",
               delta_0 = 0.8, delta_1 = 0.8, tau_0 = 2, tau_1 = 2,
               V = 2, seed = 1),
    "delta_0 and delta_1 are interpreted as delta_0u and delta_1u"
  )

  # Should warn for bounded_risk_unbounded_tau
  expect_warning(
    mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
               estimand = "ate", assumption = "bounded_risk_unbounded_tau",
               delta_0 = 0.8, delta_1 = 0.8, tau_0 = 2, tau_1 = 2,
               V = 2, seed = 1),
    "delta_0 and delta_1 are interpreted as delta_0u and delta_1u"
  )

  # Should NOT warn for other assumptions (delta aliasing is intentional there)
  expect_no_warning(
    mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
               estimand = "ate", assumption = "general",
               delta_0 = 0.8, delta_1 = 0.8,
               V = 2, seed = 1)
  )
})
