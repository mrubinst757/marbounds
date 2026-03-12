# marbounds

R package implementing bounds on causal effects when outcomes are missing under a mixture of informative and non-informative missingness (see paper: *Bounding causal effects with an unknown mixture of informative and non-informative missingness*).

## Features

- **Generic data interface**: Provide a `data.frame` with columns for outcome `Y`, treatment `A`, missingness indicator `C`, and covariates `X`.
- **Nuisance estimation**: SuperLearner with user-specified libraries and V-fold cross-fitting for propensity score, missingness probabilities, and outcome regressions.
- **User-specified estimand**: Average treatment effect (ATE), composite ATE (Ψ₁), or separable direct effect (Ψ₂).
- **User-specified assumptions**: General bounds, bounded proportion of informative missingness (δ), monotonicity (positive/negative), bounded outcome risk (τ), or point identification under known sensitivity parameters.
- **Legacy bounded-risk option**: Use `bounded_risk_unbounded_tau` to reproduce the original bounded-risk formula (without the min{1-μ, μ(τ−1)} tau constraint).
- **Tipping point analysis**: Find sensitivity parameter values (e.g. δ or τ) at which the bound crosses a target (e.g. 0), i.e. the “tipping point” that would explain away the naive association.
- **Multiplier bootstrap**: Simultaneous inference over a grid of sensitivity parameters via multiplier bootstrap.

## Installation

From source (in R):

```r
install.packages("path/to/marbounds", repos = NULL, type = "source")
```

Requires: `SuperLearner` (and its dependencies).

## Quick example

Install and load the package, then load SuperLearner (required for nuisance estimation). Simulated data: `Y` outcome, `A` treatment, `C` missingness (1 = missing), `X` covariates.

```r
library(marbounds)
suppressPackageStartupMessages(library(SuperLearner))

n <- 500
set.seed(1)
X <- runif(n, -2, 2)
A <- rbinom(n, 1, plogis(X))
C <- rbinom(n, 1, 0.2 + 0.1 * A)
Y <- ifelse(C == 1, NA, rbinom(n, 1, plogis(0.5 * A + 0.3 * X)))
dat <- data.frame(Y = Y, A = A, C = C, X = X)

# General bounds on ATE (no extra assumptions)
fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                  estimand = "ate", assumption = "general",
                  sl_lib_prop = "SL.glm", sl_lib_miss = "SL.glm", sl_lib_outcome = "SL.glm")
fit$naive
fit$lower
fit$upper

# Bounds under monotonicity with delta_1u = delta_0u = 0.8
fit2 <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                   estimand = "ate", assumption = "monotonicity_pos",
                   delta_0u = 0.8, delta_1u = 0.8,
                   sl_lib_prop = "SL.glm", sl_lib_miss = "SL.glm", sl_lib_outcome = "SL.glm")

# Tipping point: grid of (delta_0u, delta_1u) where upper bound crosses 0
tp <- tipping_point(fit2, bound_type = "upper", assumption = "monotonicity_pos")

# Multiplier bootstrap over a grid of deltas
mb <- multiplier_bootstrap_grid(fit2, param_grid = list(delta_0u = seq(0.2, 1, 0.2), delta_1u = seq(0.2, 1, 0.2)),
                                bound_spec = "upper", assumption = "monotonicity_pos", B = 500)
```

## Reference

Methods and notation follow the paper (main.tex / appendix) on bounding causal effects under mixed informative and non-informative missingness, with influence-function-based estimation and cross-fitting.
