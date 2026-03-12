#' Tipping point analysis: sensitivity parameters that explain away the naive association
#'
#' For a given bound type, finds (or evaluates on a grid) values of sensitivity parameters
#' such that the bound crosses zero (e.g. lower bound = 0 or upper bound = 0), i.e. the
#' smallest amount of informative missingness that could "explain away" the observed effect.
#' Optionally returns pointwise and uniform (multiplier bootstrap) standard errors and
#' confidence intervals for the bound at each grid point.
#'
#' @param mar_result Output from \code{mar_bounds()} (must contain \code{phi} and \code{nuisance}).
#' @param bound_type "lower" or "upper" (which bound to set to 0).
#' @param assumption Bound assumption: "general", "bounded_delta", "monotonicity_pos", "monotonicity_neg", "bounded_risk", "point_ate".
#' @param param_grid Optional list of named vectors defining the grid, e.g. \code{list(delta_0u = seq(0, 1, 0.05), delta_1u = seq(0, 1, 0.05))} or \code{list(tau = seq(1, 20, 0.5))}. If NULL, a default grid is used for the assumption.
#' @param target Value to solve for (default 0); typically the bound equals this when the association is "explained away".
#' @param inference If \code{"none"} (default), no SEs or CIs. If \code{"pointwise"}, add pointwise SE and CI. If \code{"both"}, add pointwise and uniform (multiplier bootstrap) SE and CIs.
#' @param alpha Significance level for CIs (default 0.05).
#' @param B Number of multiplier bootstrap replications for uniform CIs (default 1000). Used only when \code{inference = "both"}.
#' @param multiplier \code{"rademacher"} (default) or \code{"gaussian"} for the multiplier bootstrap. Used only when \code{inference = "both"}.
#' @return List with \code{grid} (data.frame of parameter values, bound estimates, and if requested \code{se_pointwise}, \code{ci_lower_pointwise}, \code{ci_upper_pointwise}, \code{se_uniform}, \code{ci_lower_uniform}, \code{ci_upper_uniform}), \code{tipping}, \code{target}, \code{alpha}, and if \code{inference = "both"} \code{critical_value} (1-alpha quantile of max |T*|).
#' @export
tipping_point <- function(mar_result,
                          bound_type = c("lower", "upper"),
                          assumption = c("bounded_delta", "monotonicity_pos", "monotonicity_neg", "bounded_risk", "point_ate"),
                          param_grid = NULL,
                          target = 0,
                          inference = c("none", "pointwise", "both"),
                          alpha = 0.05,
                          B = 1000,
                          multiplier = c("rademacher", "gaussian")) {
  bound_type <- match.arg(bound_type)
  assumption <- match.arg(assumption)
  inference <- match.arg(inference)
  multiplier <- match.arg(multiplier)
  phi <- mar_result$phi
  nuis <- mar_result$nuisance
  n <- nrow(phi)

  if (is.null(phi) || is.null(nuis)) stop("mar_result must contain phi and nuisance (run mar_bounds with full output).")

  # Default grids
  if (is.null(param_grid)) {
    if (assumption %in% c("bounded_delta", "monotonicity_pos", "monotonicity_neg")) {
      param_grid <- list(
        delta_0u = seq(0.05, 1, by = 0.05),
        delta_1u = seq(0.05, 1, by = 0.05)
      )
    } else if (assumption == "bounded_risk") {
      param_grid <- list(
        delta_0u = 1,
        delta_1u = 1,
        tau_0 = seq(1, 10, by = 0.25),
        tau_1 = seq(1, 10, by = 0.25)
      )
    } else {
      param_grid <- list(tau = seq(1, 20, by = 0.5), delta_0 = 1, delta_1 = 1)
    }
  }

  # Build full grid
  grid_df <- expand.grid(param_grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  n_grid <- nrow(grid_df)
  bound_vals <- numeric(n_grid)

  # Helper: get coefficient vector b for a grid row
  get_b <- function(row) {
    delta_0u <- row$delta_0u %||% 1
    delta_1u <- row$delta_1u %||% 1
    tau_0 <- row$tau_0 %||% 2
    tau_1 <- row$tau_1 %||% 2
    tau <- row$tau %||% 2
    delta_0 <- row$delta_0
    delta_1 <- row$delta_1
    if (assumption == "point_ate" && !is.null(delta_0) && !is.null(delta_1)) {
      return(coef_point_ate(delta_0, delta_1, tau))
    }
    if (assumption == "monotonicity_pos") {
      if (bound_type == "lower") return(coef_mono_pos_lower_ate(delta_0u))
      return(coef_mono_pos_upper_ate(delta_1u))
    }
    if (assumption == "monotonicity_neg") {
      if (bound_type == "lower") return(coef_mono_neg_lower_ate(delta_1u))
      return(coef_mono_neg_upper_ate(delta_0u))
    }
    if (assumption == "bounded_delta") {
      if (bound_type == "lower") return(coef_delta_lower_ate(delta_0u, delta_1u))
      return(coef_delta_upper_ate(delta_0u, delta_1u))
    }
    if (assumption == "bounded_risk") {
      if (bound_type == "lower") return(coef_bounded_risk_lower_ate(delta_0u, delta_1u, tau_0, tau_1))
      return(coef_bounded_risk_upper_ate(delta_0u, delta_1u, tau_0, tau_1))
    }
    coef_point_ate(delta_0 %||% 0.5, delta_1 %||% 0.5, tau)
  }

  # Point estimates and centered EIF matrix (n x n_grid) for inference
  psi_mat <- matrix(NA, n, n_grid)
  for (j in seq_len(n_grid)) {
    row <- as.list(grid_df[j, , drop = FALSE])
    b <- get_b(row)
    bound_vals[j] <- mean(phi %*% b)
    psi_mat[, j] <- drop(phi %*% b) - bound_vals[j]
  }
  grid_df$bound_value <- bound_vals

  # Pointwise SE and CI: se = sqrt(Var(EIF)/n), CI = estimate +/- z * se
  if (inference %in% c("pointwise", "both")) {
    se_pt <- sqrt(apply(psi_mat, 2, function(x) sum(x^2) / n) / n)
    z <- stats::qnorm(1 - alpha / 2)
    grid_df$se_pointwise <- se_pt
    grid_df$ci_lower_pointwise <- bound_vals - z * se_pt
    grid_df$ci_upper_pointwise <- bound_vals + z * se_pt
  }

  # Uniform (multiplier bootstrap) CIs: C replaces z, so CI = estimate +/- C * SE_pointwise
  if (inference == "both") {
    sigma_j <- sqrt(apply(psi_mat, 2, function(x) sum(x^2) / n))
    if (multiplier == "rademacher") {
      G_fun <- function() sample(c(-1, 1), n, replace = TRUE)
    } else {
      G_fun <- function() stats::rnorm(n)
    }
    boot_mat <- matrix(NA, B, n_grid)
    for (bb in seq_len(B)) {
      G <- G_fun()
      for (j in seq_len(n_grid)) {
        boot_mat[bb, j] <- sqrt(n) * mean(G * psi_mat[, j])
      }
    }
    # C = (1-alpha) quantile of max_j |T*_j|/sigma_j (same scale as z)
    max_std <- apply(boot_mat, 1, function(r) max(abs(r) / sigma_j))
    critical_value <- stats::quantile(max_std, probs = 1 - alpha)
    grid_df$ci_lower_uniform <- bound_vals - critical_value * se_pt
    grid_df$ci_upper_uniform <- bound_vals + critical_value * se_pt
    grid_df$se_uniform <- critical_value * se_pt / z
  }

  # Tipping: find where bound crosses target
  cross_idx <- which(diff(sign(bound_vals - target)) != 0)
  tipping <- if (length(cross_idx) > 0) {
    grid_df[cross_idx[1], , drop = FALSE]
  } else {
    NULL
  }

  out <- list(
    grid = grid_df,
    tipping = tipping,
    target = target,
    bound_type = bound_type,
    assumption = assumption,
    alpha = alpha
  )
  if (inference == "both") {
    out$critical_value <- critical_value
    out$B <- B
  }
  out
}
