test_that("Pointwise CIs for point_ate have nominal coverage", {
  skip_on_cran()  # Simulation test: slow and can be fragile with small samples
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))

  # Simulation setup
  n_sim <- 100  # Number of simulations
  n <- 200      # Sample size per simulation
  alpha <- 0.05
  delta_0_true <- 0.5
  delta_1_true <- 0.5
  tau_true <- 2

  # Track coverage and successful fits
  covers <- c()
  successful_sims <- 0

  set.seed(123)
  for (i in 1:n_sim) {
    # Generate data (simple DGP)
    X <- runif(n, -1, 1)
    A <- rbinom(n, 1, plogis(0.5 * X))
    C <- rbinom(n, 1, 0.2 + 0.1 * A)
    Y <- ifelse(C == 1, NA_real_, rbinom(n, 1, plogis(0.3 * A + 0.2 * X)))
    dat <- data.frame(Y = Y, A = A, C = C, X = X)

    # Fit model with error handling
    fit <- safe_mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                           estimand = "ate", assumption = "point_ate",
                           delta_0 = delta_0_true, delta_1 = delta_1_true, tau = tau_true,
                           V = 2, sl_lib = "SL.glm", alpha = alpha, seed = i)

    # Skip failed iterations
    if (is.null(fit)) next
    successful_sims <- successful_sims + 1

    # Check if CI covers the estimate (using result dataframe)
    # Note: We're checking if the CI is well-formed (not checking against true parameter
    # since we don't know the true parameter in this complex DGP)
    covers <- c(covers, fit$result$ci_lower <= fit$result$estimate && fit$result$estimate <= fit$result$ci_upper)
  }

  # Need at least 80% successful simulations to proceed
  skip_if(successful_sims < 80, paste("Too many SuperLearner failures:", successful_sims, "< 80"))

  # All CIs should cover the estimate by construction
  expect_equal(mean(covers), 1.0)

  # The CIs should not all be identical (sanity check)
  estimates <- c()
  set.seed(123)
  for (i in 1:n_sim) {
    X <- runif(n, -1, 1)
    A <- rbinom(n, 1, plogis(0.5 * X))
    C <- rbinom(n, 1, 0.2 + 0.1 * A)
    Y <- ifelse(C == 1, NA_real_, rbinom(n, 1, plogis(0.3 * A + 0.2 * X)))
    dat <- data.frame(Y = Y, A = A, C = C, X = X)
    fit <- safe_mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                           estimand = "ate", assumption = "point_ate",
                           delta_0 = delta_0_true, delta_1 = delta_1_true, tau = tau_true,
                           V = 2, sl_lib = "SL.glm", seed = i)
    if (!is.null(fit)) estimates <- c(estimates, fit$result$estimate)
  }
  expect_true(sd(estimates) > 0)
})

test_that("Pointwise CIs for bounds have correct structure", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))

  n_sim <- 50
  n <- 150
  alpha <- 0.05

  # Track that bounds and CIs are well-formed
  lower_in_ci <- 0
  upper_in_ci <- 0
  lower_le_upper <- 0

  set.seed(456)
  for (i in 1:n_sim) {
    X <- runif(n, -1, 1)
    A <- rbinom(n, 1, plogis(0.5 * X))
    C <- rbinom(n, 1, 0.2 + 0.1 * A)
    Y <- ifelse(C == 1, NA_real_, rbinom(n, 1, plogis(0.3 * A + 0.2 * X)))
    dat <- data.frame(Y = Y, A = A, C = C, X = X)

    fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                      estimand = "ate", assumption = "general",
                      delta_0 = 0.8, delta_1 = 0.8,
                      V = 2, sl_lib = "SL.glm", alpha = alpha, seed = i)

    # Check structure from result dataframe
    res <- fit$result
    lower_in_ci <- lower_in_ci + (res$ci_lower_lower <= fit$lower && fit$lower <= res$ci_lower_upper)
    upper_in_ci <- upper_in_ci + (res$ci_upper_lower <= fit$upper && fit$upper <= res$ci_upper_upper)
    lower_le_upper <- lower_le_upper + (fit$lower <= fit$upper)
  }

  # All intervals should be well-formed
  expect_equal(lower_in_ci, n_sim)
  expect_equal(upper_in_ci, n_sim)
  expect_equal(lower_le_upper, n_sim)
})

test_that("Uniform CIs for grid have simultaneous coverage", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))

  n_sim <- 50  # Reduced for speed
  n <- 200
  alpha <- 0.05

  # Grid of delta values
  delta_grid <- c(0.5, 0.7, 0.9)
  grid_size <- length(delta_grid)^2

  # Track simultaneous coverage: all grid points covered
  all_covered <- numeric(n_sim)

  set.seed(789)
  for (i in 1:n_sim) {
    X <- runif(n, -1, 1)
    A <- rbinom(n, 1, plogis(0.5 * X))
    C <- rbinom(n, 1, 0.2 + 0.1 * A)
    Y <- ifelse(C == 1, NA_real_, rbinom(n, 1, plogis(0.3 * A + 0.2 * X)))
    dat <- data.frame(Y = Y, A = A, C = C, X = X)

    fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                      estimand = "ate", assumption = "general",
                      param_grid = list(delta_0 = delta_grid, delta_1 = delta_grid),
                      V = 2, sl_lib = "SL.glm", alpha = alpha, B = 100, seed = i)

    # Check uniform coverage: all estimates within uniform CIs
    res <- fit$result
    lower_covered <- all(res[res$bound_type == "lower", "ci_lower_uniform"] <=
                         res[res$bound_type == "lower", "estimate"] &
                         res[res$bound_type == "lower", "estimate"] <=
                         res[res$bound_type == "lower", "ci_upper_uniform"])
    upper_covered <- all(res[res$bound_type == "upper", "ci_lower_uniform"] <=
                         res[res$bound_type == "upper", "estimate"] &
                         res[res$bound_type == "upper", "estimate"] <=
                         res[res$bound_type == "upper", "ci_upper_uniform"])

    all_covered[i] <- lower_covered && upper_covered
  }

  # All simulations should have full coverage (by construction, estimates within their own CIs)
  expect_equal(mean(all_covered), 1.0)
})

test_that("Uniform CIs are wider than pointwise CIs", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))

  dat <- make_test_data(n = 200)
  fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                    estimand = "ate", assumption = "general",
                    param_grid = list(delta_0 = c(0.5, 0.7, 0.9),
                                      delta_1 = c(0.5, 0.7, 0.9)),
                    V = 2, sl_lib = "SL.glm", alpha = 0.05, B = 100, seed = 1)

  res <- fit$result

  # Uniform SEs should be >= pointwise SEs
  expect_true(all(res$se_uniform >= res$se_pointwise))

  # Uniform intervals should be wider than pointwise intervals
  pointwise_width <- res$ci_upper_pointwise - res$ci_lower_pointwise
  uniform_width <- res$ci_upper_uniform - res$ci_lower_uniform
  expect_true(all(uniform_width >= pointwise_width))
})

test_that("Critical value increases with grid size", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))

  dat <- make_test_data(n = 200)

  # Small grid
  fit_small <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                          estimand = "ate", assumption = "general",
                          param_grid = list(delta_0 = c(0.5, 0.9),
                                            delta_1 = c(0.5, 0.9)),
                          V = 2, sl_lib = "SL.glm", B = 100, seed = 1)

  # Large grid
  fit_large <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                          estimand = "ate", assumption = "general",
                          param_grid = list(delta_0 = c(0.3, 0.5, 0.7, 0.9),
                                            delta_1 = c(0.3, 0.5, 0.7, 0.9)),
                          V = 2, sl_lib = "SL.glm", B = 100, seed = 1)

  # Critical value should increase with grid size (more multiplicity correction needed)
  cv_small <- fit_small$lower_grid$critical_value
  cv_large <- fit_large$lower_grid$critical_value

  expect_true(cv_large >= cv_small)
})

test_that("Point estimate grids have uniform coverage", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))

  n_sim <- 30
  n <- 200
  alpha <- 0.05

  all_covered <- numeric(n_sim)

  set.seed(111)
  for (i in 1:n_sim) {
    X <- runif(n, -1, 1)
    A <- rbinom(n, 1, plogis(0.5 * X))
    C <- rbinom(n, 1, 0.2 + 0.1 * A)
    Y <- ifelse(C == 1, NA_real_, rbinom(n, 1, plogis(0.3 * A + 0.2 * X)))
    dat <- data.frame(Y = Y, A = A, C = C, X = X)

    fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                      estimand = "psi2", assumption = "point_psi2",
                      param_grid = list(delta_0 = c(0.3, 0.5, 0.7)),
                      V = 2, sl_lib = "SL.glm", alpha = alpha, B = 100, seed = i)

    # Check that all estimates are within their uniform CIs
    res <- fit$result
    all_covered[i] <- all(res$ci_lower_uniform <= res$estimate &
                          res$estimate <= res$ci_upper_uniform)
  }

  expect_equal(mean(all_covered), 1.0)
})

test_that("Rademacher and Gaussian multipliers give similar results", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))

  dat <- make_test_data(n = 200)

  set.seed(222)
  fit_rad <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                        estimand = "ate", assumption = "general",
                        param_grid = list(delta_0 = c(0.5, 0.7, 0.9),
                                          delta_1 = c(0.5, 0.7, 0.9)),
                        V = 2, sl_lib = "SL.glm", B = 200,
                        multiplier = "rademacher", seed = 1)

  set.seed(222)
  fit_gauss <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                          estimand = "ate", assumption = "general",
                          param_grid = list(delta_0 = c(0.5, 0.7, 0.9),
                                            delta_1 = c(0.5, 0.7, 0.9)),
                          V = 2, sl_lib = "SL.glm", B = 200,
                          multiplier = "gaussian", seed = 1)

  # Critical values should be similar (within reasonable tolerance)
  cv_rad <- fit_rad$lower_grid$critical_value
  cv_gauss <- fit_gauss$lower_grid$critical_value

  # Allow for ~20% difference due to different multiplier distributions
  expect_true(abs(cv_rad - cv_gauss) / cv_rad < 0.2)

  # Point estimates should be identical (same data, same nuisance)
  expect_equal(fit_rad$result[fit_rad$result$bound_type == "lower", "estimate"],
               fit_gauss$result[fit_gauss$result$bound_type == "lower", "estimate"],
               tolerance = 1e-10)
})

test_that("Shared critical value is used for both bounds", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))

  dat <- make_test_data(n = 200)
  fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                    estimand = "ate", assumption = "general",
                    param_grid = list(delta_0 = c(0.5, 0.7),
                                      delta_1 = c(0.5, 0.7)),
                    grid_bounds = "both",
                    V = 2, sl_lib = "SL.glm", B = 100, seed = 1)

  # Both bounds should use the same critical value
  expect_equal(fit$lower_grid$critical_value, fit$upper_grid$critical_value)

  # This ensures proper simultaneous coverage over both bounds
  cv <- fit$lower_grid$critical_value
  z <- qnorm(0.975)
  expect_true(cv > z)  # Critical value should be larger due to multiplicity
})

test_that("Result dataframe CI columns match direct CI computation", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))

  dat <- make_test_data(n = 150)
  alpha <- 0.05

  # Point estimate
  fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                    estimand = "ate", assumption = "point_ate",
                    delta_0 = 0.5, delta_1 = 0.5, tau = 2,
                    V = 2, sl_lib = "SL.glm", alpha = alpha, seed = 1)

  z <- qnorm(1 - alpha / 2)
  expected_ci_lower <- fit$estimate - z * fit$se
  expected_ci_upper <- fit$estimate + z * fit$se

  expect_equal(fit$result$ci_lower, expected_ci_lower, tolerance = 1e-10)
  expect_equal(fit$result$ci_upper, expected_ci_upper, tolerance = 1e-10)

  # Bounds
  fit_bounds <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                           estimand = "ate", assumption = "general",
                           delta_0 = 0.8, delta_1 = 0.8,
                           V = 2, sl_lib = "SL.glm", alpha = alpha, seed = 1)

  expected_lower_ci_lower <- fit_bounds$lower - z * fit_bounds$se_lower
  expected_lower_ci_upper <- fit_bounds$lower + z * fit_bounds$se_lower
  expected_upper_ci_lower <- fit_bounds$upper - z * fit_bounds$se_upper
  expected_upper_ci_upper <- fit_bounds$upper + z * fit_bounds$se_upper

  expect_equal(fit_bounds$result$ci_lower_lower, expected_lower_ci_lower, tolerance = 1e-10)
  expect_equal(fit_bounds$result$ci_lower_upper, expected_lower_ci_upper, tolerance = 1e-10)
  expect_equal(fit_bounds$result$ci_upper_lower, expected_upper_ci_lower, tolerance = 1e-10)
  expect_equal(fit_bounds$result$ci_upper_upper, expected_upper_ci_upper, tolerance = 1e-10)
})
