# marbounds

R package implementing bounds on causal effects when outcomes are missing under a mixture of informative and non-informative missingness (see paper: *Bounding causal effects with an unknown mixture of informative and non-informative missingness* https://arxiv.org/pdf/2411.16902).

## Key features

- **Generic data interface**: Provide a `data.frame` with columns for outcome `Y`, treatment `A`, missingness indicator `C`, and covariates `X`.
- **Nuisance estimation**: SuperLearner with user-specified libraries and V-fold cross-fitting for propensity score, missingness probabilities, and outcome regressions.
- **User-specified estimand**: Average treatment effect (ATE), composite ATE (Ψ₁), or separable direct effect (Ψ₂).
- **User-specified assumptions**: General bounds, bounded proportion of informative missingness (δ), monotonicity (positive/negative), bounded outcome risk (τ), or point identification under known sensitivity parameters.
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

fit$result

fit2 <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                   assumption = "bounded_risk",
                   param_grid = list(
                     delta_0u = seq(0.5, 1, 0.1),
                     delta_1u = seq(0.5, 1, 0.1),
                     tau_1 = seq(1.5, 5, 0.5),
                     tau_0 = seq(1.5, 5, 0.5)
                   ),
                   sl_lib = "SL.glm")
# Check the dataframe output
fit2$result
```

## Reference

Methods and notation follow the paper on bounding causal effects under mixed informative and non-informative missingness, with influence-function-based estimation and cross-fitting. Note: by default, the bounded_risk option uses the method that assumes $\mu_a^\star/\mu_a \le \min(1/\mu_a, \tau_a)$, outlined in detail the Appendix. The option in the main paper may be recovered using the ``bounded_risk_unbounded_tau'' option.
