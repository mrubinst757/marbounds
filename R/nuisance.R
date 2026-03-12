#' Estimate nuisance functions using SuperLearner with cross-fitting
#'
#' Estimates e(X)=P(A=1|X), pi_0(X)=P(C=1|X,A=0), pi_1(X)=P(C=1|X,A=1),
#' mu_0(X)=E(Y|X,A=0,C=0), mu_1(X)=E(Y|X,A=1,C=0). Uses V-fold cross-fitting:
#' for each fold, fits nuisances on the training part and predicts on the held-out part.
#'
#' @param X Matrix of covariates.
#' @param A Treatment indicator (0/1).
#' @param C Missingness indicator (1=missing, 0=observed).
#' @param Y Outcome (NA or value when C=0); only used when C=0 for outcome models.
#' @param V Number of cross-fitting folds (default 2; the calling function may use
#'   \code{V=1} when a simple single-library SuperLearner such as \code{SL.glm} is used).
#' @param sl_lib_prop Character vector of SuperLearner library for propensity e (default "SL.glm").
#' @param sl_lib_miss Character vector of SuperLearner library for missingness pi_0, pi_1 (default "SL.glm").
#' @param sl_lib_outcome Character vector of SuperLearner library for outcome mu_0, mu_1 (default "SL.glm").
#' @param seed Optional seed for fold splits.
#' @return List with components: e, pi0, pi1, mu0, mu1 (each length n), and fold_id (fold index per row).
#' @export
estimate_nuisance <- function(X, A, C, Y,
                              V = 2,
                              sl_lib_prop = "SL.glm",
                              sl_lib_miss = "SL.glm",
                              sl_lib_outcome = "SL.glm",
                              seed = NULL) {
  if (!requireNamespace("SuperLearner", quietly = TRUE)) {
    stop("Package 'SuperLearner' is required. Install with install.packages('SuperLearner').")
  }
  n <- length(A)
  if (!is.matrix(X)) X <- as.matrix(X)
  V <- min(V, max(1L, n))
  if (is.null(seed)) seed <- 1L
  set.seed(seed)
  fold_id <- sample(rep(seq_len(V), length.out = n))

  e <- numeric(n)
  pi0 <- numeric(n)
  pi1 <- numeric(n)
  mu0 <- numeric(n)
  mu1 <- numeric(n)

  for (v in seq_len(V)) {
    train <- fold_id != v
    eval <- fold_id == v
    n_train <- sum(train)
    Xtrain <- X[train, , drop = FALSE]
    Xeval <- X[eval, , drop = FALSE]
    # Internal CV must have at least 1 obs per fold; use V=1 when n_train too small
    cv_v <- min(2L, max(1L, n_train))

    # Propensity e(X) = P(A=1|X) — skip SuperLearner if training set empty or too small
    if (n_train < 2L) {
      e[eval] <- clip_probs(mean(A))
    } else {
      fit_e <- SuperLearner::SuperLearner(
        Y = A[train],
        X = as.data.frame(Xtrain),
        newX = as.data.frame(Xeval),
        family = "binomial",
        SL.library = sl_lib_prop,
        cvControl = list(V = cv_v)
      )
      e[eval] <- clip_probs(fit_e$SL.predict)
    }

    # Missingness: fit P(C=1|X,A) by fitting separately in A=0 and A=1
    for (a in 0:1) {
      ia <- train & (A == a)
      n_ia <- sum(ia)
      if (n_ia < 2L) {
        if (a == 0) pi0[eval] <- mean(C[A == 0]) else pi1[eval] <- mean(C[A == 1])
        next
      }
      cv_v_a <- min(2L, max(1L, n_ia))
      fit_pi <- SuperLearner::SuperLearner(
        Y = C[ia],
        X = as.data.frame(X[ia, , drop = FALSE]),
        newX = as.data.frame(Xeval),
        family = "binomial",
        SL.library = sl_lib_miss,
        cvControl = list(V = cv_v_a)
      )
      if (a == 0) pi0[eval] <- clip_probs(fit_pi$SL.predict) else pi1[eval] <- clip_probs(fit_pi$SL.predict)
    }

    # Outcome: E(Y|X,A=a,C=0) for a=0,1
    for (a in 0:1) {
      ia <- train & (A == a) & (C == 0)
      n_ia <- sum(ia)
      if (n_ia < 2L) {
        if (a == 0) mu0[eval] <- mean(Y[A == 0 & C == 0], na.rm = TRUE) else mu1[eval] <- mean(Y[A == 1 & C == 0], na.rm = TRUE)
        next
      }
      cv_v_a <- min(2L, max(1L, n_ia))
      fit_mu <- SuperLearner::SuperLearner(
        Y = Y[ia],
        X = as.data.frame(X[ia, , drop = FALSE]),
        newX = as.data.frame(Xeval),
        family = "gaussian",
        SL.library = sl_lib_outcome,
        cvControl = list(V = cv_v_a)
      )
      pred <- fit_mu$SL.predict
      pred <- pmax(0, pmin(1, pred))  # bounded outcome
      if (a == 0) mu0[eval] <- pred else mu1[eval] <- pred
    }
  }

  list(e = e, pi0 = pi0, pi1 = pi1, mu0 = mu0, mu1 = mu1, fold_id = fold_id)
}
