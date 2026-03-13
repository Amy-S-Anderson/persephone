



# --------------------------------------------------------------------------------------------------------------
      #### PERSEPHONE ####
# --------------------------------------------------------------------------------------------------------------




#### This script contains the function command for Persephone, an agent-based model in which individual agents

#    - Are born as a starting cohort of a specified size

#    - Face a stable annual risk of forming a skeletal lesion within a specified age range, outside of which ages they cannot form a new lesion or lose an existing lesion.

#    - individuals with a skeletal lesion *may* experience a different age-specific risk of dying than individuals without skeletal lesions. When relative_mortality_risk = 1, then all individuals experience the same mortality risks, regardless of their skeletal lesion status. 

#    - Age-specific mortality risks are determined by population-specific parameter values from a Siler function




### Function: Persephone ABM -- Formation of childhood skeletal lesions, with potential for lesion-related mortality

Simulate_Cemetery <- function(cohort_size, # starting population, a named object created using the Generate.Cohort function
                                     lesion_formation_rate, # The annual probability of developing a lesion, A number between 0 and 1
                                     formation_window_opens = 0, # Age at which skeletal lesions can start forming
                                     formation_window_closes, # the oldest age at which a new lesion can form
                                     mortality_risk_type = "proportional", # One of the following: proportional, time-decreasing, time-increasing
                                     relative_mortality_risk = 1, # A number between 0 and Inf, but most likely between 1 and 5. Describes the risk of individuals with lesions, a multiplier of the risk experienced by individuals without lesions
                                     mortality_regime){ # A vector of Siler parameters based on/named for a population
  cohort <- data.frame(agent_id = 1:cohort_size, # each person gets a unique ID
                       Age = 0, # newborns
                       Lesion = 0, # no one is born with lesions
                       Dead = "No")
  
  k <- 0 # Initialize time counter
  
  # set up tables for model output: survivor data
  Alive_sum <- data.frame(Age = integer(),
                          Alive = integer(),
                          Lesion = integer(),
                          Lesion_perc = numeric())
  
  # As long as more than 10 people are alive
  while(sum(cohort$Dead == "No") >= 10){
    k <- k + 1  # Increment time
    Alive <- which(cohort$Dead == "No")
    cohort$Age[Alive] <- k  
    
    # calculate immediate mortality risk for individuals based on their current age
    age_based_risk <- (mortality_regime$a1 * exp(-mortality_regime$b1 * k) + 
                            mortality_regime$a2 + 
                            mortality_regime$a3 * exp(mortality_regime$b3 * k))
    
    for(i in Alive){
      Stress <- runif(1, 0, 1)  # chance of being exposed to lesion-causing stressor
      # New lesions form
      cohort$Lesion[i] <- ifelse(cohort$Age[i] >= formation_window_opens & # if an individual is in the right age range 
                                   cohort$Age[i] <= formation_window_closes & 
                                   Stress <= lesion_formation_rate, 1, cohort$Lesion[i]) # and is exposed to lesion-causing forces.
      
      # Now, draw a random number between 0 and 1. (Pick a card, any card...)
      death_dice <- runif(1, 0, 1)
      
      # Dead, without a lesion
      cohort$Dead[i] <- ifelse(cohort$Lesion[i] == 0 & death_dice < age_based_risk, "Yes", 
                               # Dead WITH a lesion
                               ifelse(cohort$Lesion[i] == 1 & mortality_risk_type == "proportional" & death_dice < age_based_risk * relative_mortality_risk, "Yes", 
                                      ifelse(cohort$Lesion[i] == 1 & mortality_risk_type == "time_decreasing" & death_dice < age_based_risk * relative_mortality_risk / ((cohort$Age[i] / 10) + relative_mortality_risk), "Yes", 
                                             ifelse(cohort$Lesion[i] == 1 & mortality_risk_type == "time_increasing" & death_dice < age_based_risk * ((cohort$Age[i] / 10) + relative_mortality_risk) / relative_mortality_risk, "Yes", 
                                                    cohort$Dead[i])))) #... or still alive
    }
    
        # Update summary table for Survivors
    Alive_sum <- rbind(Alive_sum, data.frame(Age = k,
                                             Alive = sum(cohort$Dead == "No" & cohort$Age == k),
                                             Lesion = sum(cohort$Dead == "No" & cohort$Lesion == 1 & cohort$Age == k),
                                             Lesion_perc = ifelse(sum(cohort$Dead == "No" & cohort$Age == k) == 0, NA, round(sum(cohort$Dead == "No" & cohort$Lesion == 1 & cohort$Age == k) / sum(cohort$Dead == "No" & cohort$Age == k) * 100, 1))))
  }
  
  # Once 10 or fewer people are left alive
  k <- k + 1  # Increment time
  Alive <- which(cohort$Dead == "No")
  cohort$Age[Alive] <- k  
  cohort <- cohort %>% select(-"Dead") # They all enter the cemetery
  
  
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
