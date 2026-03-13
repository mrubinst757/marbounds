devtools::load_all("../03-Packages/marbounds") # please keep this line for me, comment out if needed
suppressPackageStartupMessages(library(SuperLearner))
suppressPackageStartupMessages(library(tibble))

n <- 500
set.seed(1)
X <- runif(n, -2, 2)
A <- rbinom(n, 1, plogis(X))
C <- rbinom(n, 1, 0.2 + 0.1 * A)
Y <- ifelse(C == 1, NA, rbinom(n, 1, plogis(0.5 * A + 0.3 * X)))
dat <- data.frame(Y = Y, A = A, C = C, X = X)

# General bounds on ATE (no extra assumptions)
fit0 <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                  estimand = "psi2",
                  assumption = "bounded_delta",
                  delta_0 = 1,
                  sl_lib = "SL.glm",
                  family_Y = "binomial",
                  smooth_approximation = FALSE)

# General bounds on ATE: \delta = 0.8
# Testing delta aliasing: delta_0/delta_1 should be treated as delta_0u/delta_1u for bounds
fit0_a <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                   estimand = "ate",
                   assumption = "general",
                   sl_lib = "SL.glm",
                   delta_1 = 0.8,
                   delta_0 = 0.8)

c(fit0_a$naive, fit0_a$lower, fit0_a$upper)
# Check the new dataframe output
fit0_a$result

# General bounds on ATE: \delta \in {0.8, 1\}
# Testing delta aliasing in param_grid: delta_0/delta_1 should work as delta_0u/delta_1u
fit0_b <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                     estimand = "ate",
                     assumption = "general",
                     sl_lib = "SL.glm",
                     param_grid = list(
                       delta_0 = c(0.8, 1),
                       delta_1 = c(0.8, 1)
                     ))

c(fit0_b$naive, fit0_b$lower, fit0_b$upper)
# Check the combined grid result (long format) - now in 'result'
fit0_b$result

# Bounds under monotonicity with delta_1u = delta_0u = 0.8
fit1_pos <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                   assumption = "monotonicity_pos",
                   sl_lib = "SL.glm")

c(fit1_pos$naive, fit1_pos$lower, fit1_pos$upper)

fit1a_pos <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                       assumption = "monotonicity_pos",
                       delta_1u = 0.8, delta_0u = 0.8,
                       sl_lib = "SL.glm")

c(fit1a_pos$naive, fit1a_pos$lower, fit1a_pos$upper)


# Bounds under bounded_risk with tau parameters
# Testing delta aliasing: delta_0/delta_1 should work as delta_0u/delta_1u
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

fit3 <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                   estimand = "ate",
                   assumption = "point_ate",
                   delta_1 = 1, delta_0 = 1, tau = 4,
                   sl_lib = "SL.glm")

fit3$estimate
# Check the dataframe output with CI
fit3$result
args(mar_bounds)
fit4 <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                   estimand = "psi1",
                   assumption = "point_psi1",
                   delta_1 = 0.5, delta_0 = 0.5,
                   sl_lib = "SL.glm")

fit4$estimate
# Check the dataframe output with CI
fit4$result


fit5 <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                   estimand = "psi2",
                   assumption = "point_psi2",
                   delta_1 = 1, delta_0 = 0.75,
                   sl_lib = "SL.glm")

fit5$estimate
# Check the dataframe output with CI
fit5$result
