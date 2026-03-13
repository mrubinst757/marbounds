#' Uncentered influence functions for theta = (E(mu0), E(mu1), E(pi0), E(pi1), E(mu0*pi0), E(mu1*pi1))
#'
#' Notation: e = P(A=1|X), so e_1 = e, e_0 = 1-e. C=1 means missing.
#'
#' @param Y Outcome vector (use 0 or NA for missing; NA will be treated as 0 in I(C=0)*Y).
#' @param A Treatment 0/1.
#' @param C Missingness 1=missing, 0=observed.
#' @param e Propensity P(A=1|X).
#' @param pi0,pi1 P(C=1|X,A=0), P(C=1|X,A=1).
#' @param mu0,mu1 E(Y|X,A=0,C=0), E(Y|X,A=1,C=0).
#' @return Matrix n x 6: columns (phi_1,0, phi_1,1, phi_2,0, phi_2,1, phi_3,0, phi_3,1).
#' @keywords internal
influence_functions <- function(Y, A, C, e, pi0, pi1, mu0, mu1) {
  Y_ <- replace(Y, is.na(Y), 0)
  e0 <- 1 - e
  e1 <- e
  e0 <- clip_probs(e0)
  e1 <- clip_probs(e1)
  pi0 <- clip_probs(pi0)
  pi1 <- clip_probs(pi1)

  # phi_1,a: E(mu_a)
  p1_0 <- (C == 0 & A == 0) * (Y_ - mu0) / ((1 - pi0) * e0) + mu0
  p1_1 <- (C == 0 & A == 1) * (Y_ - mu1) / ((1 - pi1) * e1) + mu1

  # phi_2,a: E(pi_a)
  p2_0 <- (A == 0) * (C - pi0) / e0 + pi0
  p2_1 <- (A == 1) * (C - pi1) / e1 + pi1

  # phi_3,a: E(mu_a * pi_a)
  p3_0 <- p1_0 * pi0 + p2_0 * mu0 - mu0 * pi0
  p3_1 <- p1_1 * pi1 + p2_1 * mu1 - mu1 * pi1

  cbind(phi_1_0 = p1_0, phi_1_1 = p1_1, phi_2_0 = p2_0, phi_2_1 = p2_1, phi_3_0 = p3_0, phi_3_1 = p3_1)
}

#' Compute empirical mean of influence functions (point estimates of theta components)
#'
#' @param phi n x 6 matrix from influence_functions().
#' @return Vector of length 6: (m1_0, m1_1, m2_0, m2_1, m3_0, m3_1).
theta_hat <- function(phi) {
  colMeans(phi)
}

#' Empirical variance matrix of phi (for asymptotic variance of b' theta)
#'
#' @param phi n x 6 matrix.
#' @return 6 x 6 variance matrix.
phi_var <- function(phi) {
  n <- nrow(phi)
  stats::var(phi) * (n - 1) / n
}
