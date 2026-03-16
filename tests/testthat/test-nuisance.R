test_that("estimate_nuisance returns list with e, pi0, pi1, mu0, mu1, fold_id", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_tiny_data(n = 80)
  nuis <- estimate_nuisance(
    X = as.matrix(dat$X), A = dat$A, C = dat$C, Y = dat$Y,
    V = 2, sl_lib_prop = "SL.glm", sl_lib_miss = "SL.glm", sl_lib_outcome = "SL.glm",
    seed = 1
  )
  expect_named(nuis, c("e", "pi0", "pi1", "mu0", "mu1", "fold_id"))
  expect_length(nuis$e, 80)
  expect_true(all(nuis$e > 0 & nuis$e < 1))
  expect_true(all(nuis$pi0 >= 0 & nuis$pi0 <= 1))
  expect_true(all(nuis$pi1 >= 0 & nuis$pi1 <= 1))
  expect_true(all(nuis$mu0 >= 0 & nuis$mu0 <= 1))
  expect_true(all(nuis$mu1 >= 0 & nuis$mu1 <= 1))
  expect_equal(length(nuis$fold_id), 80)
  expect_true(all(nuis$fold_id %in% c(1, 2)))
})

test_that("estimate_nuisance runs with default libraries", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_tiny_data(n = 100)
  nuis <- estimate_nuisance(
    X = as.matrix(dat$X), A = dat$A, C = dat$C, Y = dat$Y,
    V = 2, seed = 1
  )
  expect_length(nuis$e, 100)
})
