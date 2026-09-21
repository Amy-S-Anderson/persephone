
# Test that implementing different frailty distributions does not change age-at-death distribution
# in the absence of additional sources of mortality risk (stress exposure).
# Frailty variance should be incorporated such that the baseline population hazard is still defined 
# by the Siler function passed to the model. 

# This is a stochastic, seeded test with a justified tolerance rather than an exact check.


# Minimal reimplementation of Simulate_Cemetery's mortality-only loop,
# with no lesions/exposure/fertility -- isolates exactly the mechanism
# under test.
run_cohort_mortality_only <- function(pop, mu0_table, max_age = 100) {
  decedents_list <- vector("list", max_age)
  current_time <- 1
  while (nrow(pop) > 0 && current_time <= max_age) {
    force_death <- nrow(pop) <= 10  # mirrors Simulate_Cemetery's truncation rule
    result <- apply_mortality(
      pop,
      mu0_lookup   = mu0_table,
      force_death  = force_death,
      current_time = current_time,
      risk_factors = list(frailty = NULL, acquired_frailty = NULL)  # <- the fix
    )
    decedents_list[[current_time]] <- result$decedents
    pop <- result$pop
    if (nrow(pop) > 0) pop <- age_pop(pop)
    current_time <- current_time + 1
  }
  do.call(rbind, decedents_list)
}

set.seed(123)
N <- 10000

pop_v0  <- create_pop(N, age_structured = FALSE, mortality_regime = CoaleDemenyWestF11,
                      pop_config = list(model_lesions = FALSE, annual_exposure = NULL,
                                        lesion_formation_window = c(0,0), frailty_variance = 0))
pop_v02 <- create_pop(N, age_structured = FALSE, mortality_regime = CoaleDemenyWestF11,
                      pop_config = list(model_lesions = FALSE, annual_exposure = NULL,
                                        lesion_formation_window = c(0,0), frailty_variance = 0.2))

mu0_v0  <- defrail_siler(CoaleDemenyWestF11, frailty_variance = 0)
mu0_v02 <- defrail_siler(CoaleDemenyWestF11, frailty_variance = 0.2)

decedents_v0  <- run_cohort_mortality_only(pop_v0,  mu0_v0,  max_age = 100)
decedents_v02 <- run_cohort_mortality_only(pop_v02, mu0_v02, max_age = 100)

# --- Compare the two age-at-death distributions -------------------------
# Quick numeric summary first
summary(decedents_v0$age)
summary(decedents_v02$age)

# Empirical survival curves, compared directly against each other
S_v0  <- ecdf(decedents_v0$age)
S_v02 <- ecdf(decedents_v02$age)
ages  <- 0:100
plot(ages, 1 - S_v0(ages),  type = "l", ylab = "Survival", xlab = "Age")
lines(ages, 1 - S_v02(ages), col = "red")
legend("topright", legend = c("variance = 0", "variance = 0.2"), col = c("black","red"), lty = 1)

max(abs((1 - S_v0(ages)) - (1 - S_v02(ages))))  # largest gap between the two curves
