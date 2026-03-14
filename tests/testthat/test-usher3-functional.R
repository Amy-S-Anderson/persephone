# Functional tests for Usher3 model
#
# These tests verify that persephone's forward simulation produces samples
# from the same distribution as demohaz's rejection sampler. We use two-sample
# Kolmogorov-Smirnov tests to compare the empirical distributions.

library(testthat)
library(demohaz)
library(persephone)

# -----------------------------------------------------------------------------
# Shared test parameters (matching demohaz test-usher3-simulation-functional.R)
# -----------------------------------------------------------------------------

# Baseline parameter vector (LEH example):
#   th0 = [k1, k2, b_siler]  where b_siler has 5 Siler hazard parameters.
#   k1  = 0.02  : well-to-ill transition rate
#   k2  = 1.2   : mortality multiplier for ill individuals
b_siler_test <- c(.175, 1.40, .368 * .01,
                  log(.917 * .1 / (.075 * .001)) / (.917 * .1),
                  .917 * .1)
th0 <- c(2e-2, 1.2, b_siler_test)

# Extract individual parameters
k1 <- th0[1]
k2 <- th0[2]
b_siler <- th0[3:7]

# Persephone parameter vectors
mortality_param <- c(k2, b_siler)

# -----------------------------------------------------------------------------
# Forward simulation helper
#
# Runs persephone's modular simulation (apply_lesions + apply_mortality_usher3)
# until all agents have died, returning death ages and lesion status.
# -----------------------------------------------------------------------------
simulate_persephone_forward <- function(N, k1, mortality_param, dx, xmax,
                                         lesion_model = "constant",
                                         x_cut = Inf) {
  # Build lesion_param based on model

if (lesion_model == "constant") {
    lesion_param <- c(k1)
  } else if (lesion_model == "constant_to") {
    lesion_param <- c(k1, x_cut)
  } else {
    stop("Unsupported lesion_model: ", lesion_model)
  }

  # Create initial cohort: all alive, age 0, no lesions
  state <- data.frame(
    agent_id = 1:N,
    age = rep(0, N),
    lesion = rep(FALSE, N),
    dead = rep(FALSE, N),
    in_sample = rep(TRUE, N)
  )

  nsteps <- ceiling(xmax / dx)

  for (step in seq_len(nsteps)) {
    # Check if anyone is still alive
    if (all(state$dead)) break

    # Record who is alive before this step
    alive_before <- !state$dead

    # Apply lesion acquisition (before mortality, at current age)
    state <- apply_lesions(
      state = state,
      lesion_model = lesion_model,
      lesion_param = lesion_param,
      dx = dx
    )

    # Apply mortality (also advances age for survivors)
    state <- apply_mortality_usher3(
      state = state,
      mortality_model = "usher3",
      mortality_param = mortality_param,
      dx = dx
    )

    # Add uniform jitter to death ages to smear across the timestep bin.
    # This eliminates ties in the KS test.
    # Agents who died this step: were alive before, now dead
    died_this_step <- alive_before & state$dead
    n_died <- sum(died_this_step)
    if (n_died > 0) {
      state$age[died_this_step] <- state$age[died_this_step] + runif(n_died, 0, dx)
    }
  }

  # Return death ages and lesion status
  list(
    x = state$age,
    ill = state$lesion
  )
}

# -----------------------------------------------------------------------------
# Test 1: x_cut = Inf (lesion formation at all ages)
#
# Compare persephone forward simulation against demohaz rejection sampler.
# -----------------------------------------------------------------------------
test_that("persephone forward simulation matches demohaz rejection sampler (x_cut = Inf)", {
  set.seed(12345)
  N <- 20000       # sample size for each method
  dx <- 0.01       # forward simulation time step
  dx_rej <- 0.001  # rejection sampler grid spacing
  xmax <- 120      # maximum age

  # Run persephone forward simulation
  perse <- simulate_persephone_forward(
    N = N,
    k1 = k1,
    mortality_param = mortality_param,
    dx = dx,
    xmax = xmax,
    lesion_model = "constant",
    x_cut = Inf
  )

  # Run demohaz rejection sampler
  demo <- demohaz::sample_usher3(N, th0, dx_rej, xmax, x_cut = Inf)

  # Split each sample by lesion/illness status at death
  perse_well <- perse$x[!perse$ill]
  perse_ill  <- perse$x[perse$ill]
  demo_well  <- demo$x[!demo$ill]
  demo_ill   <- demo$x[demo$ill]

  # Two-sample KS tests: if both methods sample from the same distribution,
  # p-values should be non-significant (> 0.05) with high probability.
  expect_gt(ks.test(perse_well, demo_well)$p.value, 0.05)
  expect_gt(ks.test(perse_ill, demo_ill)$p.value, 0.05)

  # The fraction of individuals who are lesioned at death should agree closely
  # between the two methods.
  prop_perse <- mean(perse$ill)
  prop_demo <- mean(demo$ill)
  expect_lt(abs(prop_perse - prop_demo), 0.03)
})

# -----------------------------------------------------------------------------
# Test 2: x_cut = 6 (lesion formation only before age 6, e.g. LEH)
#
# Compare persephone forward simulation against demohaz rejection sampler.
# -----------------------------------------------------------------------------
test_that("persephone forward simulation matches demohaz rejection sampler (x_cut = 6)", {
  set.seed(67890)
  N <- 10000       # sample size for each method
  dx <- 0.01       # forward simulation time step
  dx_rej <- 0.001  # rejection sampler grid spacing
  xmax <- 120      # maximum age
  x_cut <- 6

  # Run persephone forward simulation
  perse <- simulate_persephone_forward(
    N = N,
    k1 = k1,
    mortality_param = mortality_param,
    dx = dx,
    xmax = xmax,
    lesion_model = "constant_to",
    x_cut = x_cut
  )

  # Run demohaz rejection sampler
  demo <- demohaz::sample_usher3(N, th0, dx_rej, xmax, x_cut = x_cut)

  # Split each sample by lesion/illness status at death
  perse_well <- perse$x[!perse$ill]
  perse_ill  <- perse$x[perse$ill]
  demo_well  <- demo$x[!demo$ill]
  demo_ill   <- demo$x[demo$ill]

  # Two-sample KS tests
  expect_gt(ks.test(perse_well, demo_well)$p.value, 0.05)
  expect_gt(ks.test(perse_ill, demo_ill)$p.value, 0.05)

  # Proportion lesioned at death
  prop_perse <- mean(perse$ill)
  prop_demo <- mean(demo$ill)
  expect_lt(abs(prop_perse - prop_demo), 0.03)
})

