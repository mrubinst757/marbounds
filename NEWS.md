# marbounds 0.1.0

Initial release of the marbounds package for bounding causal effects under mixed informative and non-informative missingness.

## Features

* `mar_bounds()` - Main function for computing bounds on causal effects
* Support for multiple estimands: ATE, composite ATE (Ψ₁), separable direct effect (Ψ₂)
* Multiple assumption types:
  - General bounds (no assumptions)
  - Bounded proportion of informative missingness (δ)
  - Monotonicity (positive/negative)
  - Bounded outcome risk (τ)
  - Point identification under known sensitivity parameters
* SuperLearner integration for nuisance parameter estimation with V-fold cross-fitting
* Multiplier bootstrap for simultaneous inference over parameter grids
* Influence function-based estimators with asymptotic standard errors

## Reference

Rubinstein, M., Agniel, D., Han, L., Horvitz-Lennon, M., & Normand, S.-L. (2026).
Bounding causal effects with an unknown mixture of informative and non-informative missingness.
Journal of Causal Inference (Accepted). https://arxiv.org/pdf/2411.16902
