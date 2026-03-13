



# --------------------------------------------------------------------------------------------------------------
      #### PERSEPHONE ####
# --------------------------------------------------------------------------------------------------------------




#### This script contains the function command for Persephone, an agent-based model in which individual agents

#    - Are born as a starting cohort of a specified size

#    - Face a stable annual risk of forming a skeletal lesion within a specified age range, outside of which ages they cannot form a new lesion or lose an existing lesion.

#    - individuals with a skeletal lesion *may* experience a different age-specific risk of dying than individuals without skeletal lesions. When relative_mortality_risk = 1, then all individuals experience the same mortality risks, regardless of their skeletal lesion status. 

#    - Age-specific mortality risks are determined by population-specific parameter values from a Siler function




### Helper: Create the initial cohort data frame
create_cohort <- function(cohort_size) {
  data.frame(agent_id = 1:cohort_size, # each person gets a unique ID
             age = 0,                   # newborns
             lesion = 0,                # no one is born with lesions
             dead = FALSE,
             in_sample = TRUE) # all agents start in the sample
}

### Helper: Age up all living agents to the current timestep
age_cohort <- function(cohort, k) {
  Alive <- which(!cohort$dead)
  cohort$age[Alive] <- k
  cohort
}

### Helper: Compute Siler hazard for a given age
compute_siler_risk <- function(k, mortality_regime) {
  mortality_regime$a1 * exp(-mortality_regime$b1 * k) +
    mortality_regime$a2 +
    mortality_regime$a3 * exp(mortality_regime$b3 * k)
}

### Helper: For a single agent, roll for lesion formation
form_lesion <- function(cohort, i, formation_window_opens, formation_window_closes, lesion_formation_rate) {
  Stress <- runif(1, 0, 1)  # chance of being exposed to lesion-causing stressor
  cohort$lesion[i] <- ifelse(cohort$age[i] >= formation_window_opens &
                               cohort$age[i] <= formation_window_closes &
                               Stress <= lesion_formation_rate, 1, cohort$lesion[i])
  cohort
}

### Helper: For a single agent, roll for death
apply_mortality <- function(cohort, i, age_based_risk, mortality_risk_type, relative_mortality_risk) {
  death_dice <- runif(1, 0, 1)
  cohort$dead[i] <- ifelse(cohort$lesion[i] == 0 & death_dice < age_based_risk, TRUE,
                           ifelse(cohort$lesion[i] == 1 & mortality_risk_type == "proportional" & death_dice < age_based_risk * relative_mortality_risk, TRUE,
                                  ifelse(cohort$lesion[i] == 1 & mortality_risk_type == "time_decreasing" & death_dice < age_based_risk * relative_mortality_risk / ((cohort$age[i] / 10) + relative_mortality_risk), TRUE,
                                         ifelse(cohort$lesion[i] == 1 & mortality_risk_type == "time_increasing" & death_dice < age_based_risk * ((cohort$age[i] / 10) + relative_mortality_risk) / relative_mortality_risk, TRUE,
                                                cohort$dead[i]))))
  cohort
}

### Helper: Record survivor snapshot for the current timestep
record_survivors <- function(cohort, k) {
  n_alive <- sum(!cohort$dead & cohort$age == k)
  n_lesion <- sum(!cohort$dead & cohort$lesion == 1 & cohort$age == k)
  data.frame(Age = k,
             Alive = n_alive,
             Lesion = n_lesion,
             Lesion_perc = ifelse(n_alive == 0, NA, round(n_lesion / n_alive * 100, 1)))
}


### Function: Persephone ABM -- Formation of childhood skeletal lesions, with potential for lesion-related mortality

Simulate_Cemetery <- function(cohort_size, # starting population, a named object created using the Generate.Cohort function
                                     lesion_formation_rate, # The annual probability of developing a lesion, A number between 0 and 1
                                     formation_window_opens = 0, # Age at which skeletal lesions can start forming
                                     formation_window_closes, # the oldest age at which a new lesion can form
                                     mortality_risk_type = "proportional", # One of the following: proportional, time-decreasing, time-increasing
                                     relative_mortality_risk = 1, # A number between 0 and Inf, but most likely between 1 and 5. Describes the risk of individuals with lesions, a multiplier of the risk experienced by individuals without lesions
                                     mortality_regime){ # A vector of Siler parameters based on/named for a population
  cohort <- create_cohort(cohort_size)
  
  k <- 0 # Initialize time counter
  
  # set up tables for model output: survivor data
  Alive_sum <- data.frame(Age = integer(),
                          Alive = integer(),
                          Lesion = integer(),
                          Lesion_perc = numeric())
  
  # As long as more than 10 people are alive
  while(sum(!cohort$dead) >= 10){
    k <- k + 1  # Increment time
    cohort <- age_cohort(cohort, k)
    Alive <- which(!cohort$dead)

    # calculate immediate mortality risk for individuals based on their current age
    age_based_risk <- compute_siler_risk(k, mortality_regime)
    
    for(i in Alive){
      cohort <- form_lesion(cohort, i, formation_window_opens, formation_window_closes, lesion_formation_rate)
      cohort <- apply_mortality(cohort, i, age_based_risk, mortality_risk_type, relative_mortality_risk)
    }

    # Update summary table for Survivors
    Alive_sum <- rbind(Alive_sum, record_survivors(cohort, k))
  }
  
  # Once 10 or fewer people are left alive
  k <- k + 1  # Increment time
  Alive <- which(!cohort$dead)
  cohort$age[Alive] <- k
  cohort <- cohort %>% select(-"dead") # They all enter the cemetery
  
  
  # Model output: data frame of all individuals with lesion status and age at death; Frequency table of Cemetery ages and lesion counts; Frequency table of Survivors at each age, and their lesion counts.
  output <- list(individual_outcomes = cohort, survivors = Alive_sum)
  return(output)
}






test <- Simulate_Cemetery(cohort_size = 500,
                                 lesion_formation_rate = .10,
                                 formation_window_closes = 5,
                                 mortality_regime = CoaleDemenyWestF5,
                                 mortality_risk_type = "proportional"
                                 )
