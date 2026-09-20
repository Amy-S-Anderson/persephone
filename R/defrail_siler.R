#' Helper functions for calculating mortality risk

#'@title Siler hazard (helper)
#'@description Compute Siler hazard for a given age. This helper function is referenced in apply_mortality(). 
#'
#' @param ages Integer vector, if population is age-structured
#' @param mortality_regime Data frame with Siler parameters (a1, b1, a2, a3, b3)
#' @return Numeric hazard value
#' @keywords internal
compute_siler_risk <- function(ages, mortality_regime) { 
  if (any(mortality_regime < 0)) {
    stop('No parameters of the mortality regime can be negative')
  }
  if(mortality_regime$a3 < 1){ # These are traditional Siler parameters
    age_based_risk <- mortality_regime$a1 * exp(-mortality_regime$b1 * ages) +
      mortality_regime$a2 +
      mortality_regime$a3 * exp(mortality_regime$b3 * ages)
  }
  if(mortality_regime$a3 > 1){ # These are robust Siler parameters
    mortality_regime <- demohaz_to_trad_siler_param(mortality_regime)
    age_based_risk <- mortality_regime$a1 * exp(-mortality_regime$b1 * ages) +
      mortality_regime$a2 +
      mortality_regime$a3 * exp(mortality_regime$b3 * ages)
  }
  age_based_risk
}


#' @title Siler survivorship (helper)
#' @description Compute discrete Siler survivorship l(x). Returns the probability of surviving from birth to each age x, under the Siler hazard model. l(0) is defined as 1 (everyone is alive at birth).
#'
#' @param ages Integer vector of ages (typically 0:max_age)
#' @param mortality_regime Data frame with Siler parameters (a1, b1, a2, a3, b3)
#' @return Numeric vector of survivorship values, same length as ages
#' @keywords internal
compute_siler_survivorship <- function(ages, mortality_regime) {
  hazards        <- compute_siler_risk(ages, mortality_regime)
  survival_probs <- 1 - hazards
  cumprod(c(1, survival_probs[-length(survival_probs)]))
}



# Helper: Solves the defrailing ODE via RK4 on a fine age grid running from 0 to
# max_survivable_age. Returns the full fine-resolution solution, before
# any midpoint-offset sampling or trimming to integer ages happens --
# kept separate from defrail_siler() specifically so this core numerical
# step can be tested on its own, independent of how its output later gets
# packaged into a lookup table.
solve_defrail_ode <- function(mortality_regime, s2, max_survivable_age, step) {
  dH0_da <- function(a, H0) {
    mu_bar <- compute_siler_risk(a, mortality_regime)
    mu_bar * (1 + s2 * H0)
  }
  
  fine_ages <- seq(0, max_survivable_age, by = step)
  n         <- length(fine_ages)
  H0_fine   <- numeric(n)
  mu0_fine  <- numeric(n)
  
  H0_fine[1]  <- 0
  mu0_fine[1] <- dH0_da(0, 0)
  
  for (i in seq_len(n - 1)) {
    a  <- fine_ages[i]   # current age (start of this step)
    H  <- H0_fine[i]     # current cumulative hazard, H0(a)
    
    # k1: the slope RIGHT NOW, at the start of the interval.
    # This is exactly what Euler's method would use on its own.
    k1 <- dH0_da(a, H)
    
    # k2: the slope at the MIDPOINT of the interval (a + step/2), but
    # estimated using a value of H0 that itself only comes from a half-step
    # of Euler's method using k1. In other words: "if I trusted k1 for half
    # a step, where would I be, and what's the slope there?"
    # step / 2 appears twice here for the same reason: once to move a to
    # the midpoint age, and once to move H half a step forward using k1.
    k2 <- dH0_da(a + step / 2,   H + step * k1 / 2)
    
    # k3: ANOTHER slope estimate at the same midpoint age, but this time
    # using k2 (not k1) to estimate H at the midpoint. This refines the
    # midpoint estimate -- k2 was a first guess at the midpoint slope using
    # a fairly crude H estimate; k3 redoes it using a better H estimate.
    k3 <- dH0_da(a + step / 2,   H + step * k2 / 2)
    
    # k4: the slope at the END of the interval (a + step), using H
    # projected a full step forward via k3.
    k4 <- dH0_da(a + step,       H + step * k3)
    
    # The actual update: a WEIGHTED AVERAGE of the four slope estimates,
    # not a simple average. The weights (1, 2, 2, 1), summing to 6 (hence
    # dividing by 6), aren't arbitrary -- they come from matching this
    # formula to a Taylor series expansion of the true solution, such that
    # errors up to 4th order cancel out. Intuitively: the two midpoint
    # estimates (k2, k3) are considered more representative of the
    # interval's "typical" slope than the endpoint estimates (k1, k4), so
    # they get double weight -- this is the same weighting pattern as
    # Simpson's rule for numerical integration, which isn't a coincidence.
    H0_fine[i + 1]  <- H + (step / 6) * (k1 + 2*k2 + 2*k3 + k4)
    
    # Once we know H0 at the new age, we can get mu0 there too: mu0 IS
    # dH0_da evaluated at that point (see the function's own definition --
    # mu0(a) = dH0/da), so this just calls the same right-hand-side function
    # again at the new (a + step, H0) pair, rather than deriving it
    # differently.
    mu0_fine[i + 1] <- dH0_da(a + step, H0_fine[i + 1])
  }
  
  
  
  list(fine_ages = fine_ages, H0_fine = H0_fine, mu0_fine = mu0_fine)
}







#' Defrail a Siler mortality regime
#'
#' Given a Siler mortality regime whose parameters were fit to observed
#' (population-aggregate) mortality data, compute the underlying individual-
#' level baseline hazard mu_0(a) that, when mixed over a Gamma frailty
#' distribution with mean 1 and variance frailty_variance, reproduces the
#' observed hazard. Returns a lookup table of mu_0 by integer age.
#'
#' The Gamma-frailty mixture relationship is:
#'   mu_bar(a) = mu_0(a) / (1 + sigma^2 * H_0(a))
#'   Reading it in words: the population-average hazard you'd observe at age a 
#'   equals the individual-level baseline hazard, divided by a "selection discount" 
#'   that grows over time. The denominator, 1 + σ²·H₀(a), is exactly that discount —
#'    it starts at 1 (no discount at birth, since H₀(0)=0, nobody has died yet, 
#'    no selection has happened), and grows as H₀(a) accumulates, meaning the 
#'    observed population hazard falls further and further behind the true 
#'    individual hazard as the survivors skew healthier.
#'    σ² (the frailty variance) controls how strong this effect is.
#'    
#'    
#' The Gamma-frailty mixture equation for mortality hazard implies the ODE:
#'   dH_0/da = mu_bar(a) * (1 + sigma^2 * H_0(a)),  H_0(0) = 0
#' solved numerically via RK4, then mu_0(a) = dH_0/da.
#' mu_bar(a) is the Siler value at age a. sigma is frailty variance. 
#' H_0(a) is the target of estimation. dH_0/da is its derivative: the slope, or instantaneous hazard at age a. 
#'
#' @param mortality_regime Data frame with Siler parameters (a1, b1, a2, a3, b3)
#' @param frailty_variance Numeric scalar, variance of the Gamma frailty
#'   distribution (mean is fixed at 1). Set to 0 to recover the observed
#'   hazard unchanged.
#' @param max_age Integer, oldest age to compute (default 120)
#' @param step Numeric, integration step size in years (default 0.01; finer
#'   steps improve accuracy at the cost of speed)
#'
#' @return A data frame with columns:
#'   \item{age}{Integer ages 0:max_age}
#'   \item{mu0}{Defrailted baseline hazard at each age}
#' @export
defrail_siler <- function(mortality_regime,
                          frailty_variance,
                          max_age             = 110,
                          max_survivable_age  = 110,
                          step                = 0.01) {
  
  
  # If frailty_variance == 0, or if frailty_variance is not specified, the defrailed hazard is just the observed hazard
  # --- σ² == 0 branch (currently samples raw Siler hazard at exact integer ages) ---
  if (frailty_variance == 0 || is.null(frailty_variance)) {
    ages <- 0:max_age
    return(data.frame(
      age = ages,
      mu0 = compute_siler_risk(ages + 0.5, mortality_regime)  # <- 0.5 = year midpoint offset to address the discretization of hazards problem, which introduces error into cumulative hazards by only summing hazards from the beginning of every year. 
      # this helper function (compute_siler_risk) is in mortality.R
    ))
  }
  
  if (frailty_variance < 0) {
    stop("frailty_variance must be non-negative")
  }
  
  
  
  
  # Case when frailty variance > 0:
  s2 <- frailty_variance
  # pass to helper function, solve_defrail_ode
  ode <- solve_defrail_ode(mortality_regime, s2, max_survivable_age, step)
  fine_ages <- ode$fine_ages
  H0_fine   <- ode$H0_fine
  mu0_fine  <- ode$mu0_fine
  n         <- length(fine_ages)
  
  # Guard: check that defrailed hazard kills off the population by age 110
  # Guard: check survivorship at max_survivable_age (110)
  max_survival_at_cap <- 0.001
  idx_cap             <- round(max_survivable_age / step) + 1
  H0_at_cap           <- H0_fine[min(idx_cap, n)]
  survival_at_cap     <- exp(-H0_at_cap)
  
  if (survival_at_cap > max_survival_at_cap) {
    stop(sprintf(
      paste0("The combination of this mortality regime and frailty_variance = %.4f ",
             "produces a defrailed baseline hazard under which %.2f%% of agents ",
             "would survive to age %d",
             "Reduce frailty_variance so that the low-frailty tail of the population ",
             "does not produce implausibly long lifespans."),
      frailty_variance,
      survival_at_cap * 100,
      max_survivable_age
    ))
  }
  
  #   # Sample mu0 at integer ages for the lookup table. Clip output table to max_survivable_age
  
  
  # --- σ² > 0 branch: sample mu0_fine at the midpoint of each age-year, not its start ---
  integer_ages     <- 0:min(max_age, max_survivable_age)
  midpoint_indices <- round((integer_ages + 0.5) / step) + 1   # <- offset by half a year
  midpoint_indices <- pmin(midpoint_indices, n)
  
  data.frame(
    age = integer_ages,
    mu0 = mu0_fine[midpoint_indices]
  )
  
}



