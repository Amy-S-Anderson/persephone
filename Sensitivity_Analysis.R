

#### Known Unknowns: Sensitivity Analysis ####

#### This script 
# 1. runs a parameter sweep: Sweeping across lesion exposure rate (1-10% annually in each of the first ten years) and mortality regime (Coale and Demeny West model life tables), for (1) risk-free lesions and (2) risk-doubling lesions
# 2. plots the results of the sweep in a series of multi-panel line plots
# 3. 


### Probability of exposure to lesion-causing conditions (% annual probability of exposure in each of the first ten years) ###
# - 1, 2, 3, 4, 5, 6, 7, 8, 9, 10.

### Mortality curve (Coale & Demeny West model life tables for females) ###
# - CDW3, 5, 7, 9, 11, 13, 15




#### Load libraries ####
library(here) # for setting the working directory relative to the project file.

# packages for data wrangling
library(tidyverse)
library(dplyr)
library(purrr) # for working with lists

# packages for survival analysis
library(survival)
library(survminer)
library(survRM2) # for calculating restricted mean survival time

# packages for data visualization
library(ggplot2)
library(cowplot)

# packages for handling model output
library(jsonlite)
library(data.table)




# ---------------
# 1. Set up parallel processing to reduce run time.
# ---------------
library(doParallel)
library(foreach)

# Detect available cores and register
ncores <- parallel::detectCores() - 1  # leave one free
cl <- makeCluster(ncores)
registerDoParallel(cl)





# ---------------
#### MORTALITY REGIMES ####
# 2. Define a series of Siler functions fit to Coale & Demeny West model life tables for females.
# ---------------

# Fast mortality
CoaleDemenyWestF3 <- data.frame(
  a1= 0.558,
  b1= 1.05,
  a2= 0.01225, 
  a3= 0.000520, 
  b3= 0.0727,
  name="CoaleDemenyWestF3") 

CoaleDemenyWestF5 <- data.frame(
  a1= 0.457,
  b1= 1.07,
  a2= 0.01037, 
  a3= 0.000359, 
  b3= 0.0763,
  name="CoaleDemenyWestF5") 

CoaleDemenyWestF11 <- data.frame(
  a1= 0.256,
  b1= 1.17,
  a2= 0.00596, 
  a3= 0.000133, 
  b3= 0.086,
  name="CoaleDemenyWestF11") 

CoaleDemenyWestF15 <- data.frame(
  a1= 0.175,
  b1= 1.40,
  a2= 0.00368, 
  a3= 0.000075, 
  b3= 0.0917,
  name="CoaleDemenyWestF15") 

CoaleDemenyWestF17 <- data.frame(
  a1= 0.14,
  b1= 1.57,
  a2= 0.00265, 
  a3= 0.000056, 
  b3= 0.0949,
  name="CoaleDemenyWestF17") 

# Slow mortality
CoaleDemenyWestF21 <- data.frame(
  a1= 0.091,
  b1= 2.78,
  a2= 0.00092, 
  a3= 0.000025, 
  b3= 0.1033,
  name="CoaleDemenyWestF21") 



#### Load Anderson's agent-based model ####
source(here("Known Unknowns/Model_Core_Simulate_Cemetery.R")) # <- See this file for model details. 
#### and functions for running and reading a model sweep ####
source(here("Known Unknowns/Sweep_Utility_Functions.R"))
#### and functions plotting results of model sweep ####
source(here("Known Unknowns/Sweep_Plotting_Functions.R"))






# -----------------------------------------------------------------------------------------------------------------------------------------

# ---------------
#### SETUP ####
# 3. Define the parameter values to sweep across 
# ---------------

# Define default parameters
params <- get_default_params()

# Sweep configuration
reps <- 100 # The number of iterations to run each model scenario. Identical parameters, different random seeds for each run. 
lesion_exposure_rates <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1)
mortality_regimes <- list(
  CoaleDemenyWestF3,
  CoaleDemenyWestF5,
  CoaleDemenyWestF11,
  CoaleDemenyWestF15,
  CoaleDemenyWestF17,
  CoaleDemenyWestF21
)




# ---------------
#### RUN AND SAVE SWEEPS ####
# 4. Run many iterations of the model, systematically changing key parameters (defined above)
# ---------------
## These loops will take some time to execute. 
 
 # Run parameter sweep in parallel and return a single combined data frame
# sweep_results_rmr1 <- foreach(regime = mortality_regimes,
#                                .packages = c("dplyr", "data.table", "jsonlite"),
#                                .combine = rbind) %dopar% {
#                                  
#                                  # Run sweep for this mortality regime
#                                  result_df <- run_sweep(
#                                    mortality = regime,
#                                    lesion_rates = lesion_exposure_rates,
#                                    reps = reps,
#                                    rmr = 1 # risk-free lesions
#                                  )
#                                  
#                                  # Add a column for the mortality regime name
#                                  result_df$mortality <- regime$name
#                                  
#                                  # Return the data frame for this iteration
#                                  return(result_df)
#                                }
# 
# sweep_results_rmr2 <- foreach(regime = mortality_regimes,
#                                .packages = c("dplyr", "data.table", "jsonlite"),
#                                .combine = rbind) %dopar% {
#                                  
#                                  # Run sweep for this mortality regime
#                                  result_df <- run_sweep(
#                                    mortality = regime,
#                                    lesion_rates = lesion_exposure_rates,
#                                    reps = reps,
#                                    rmr = 2 # risk-doubling lesions
#                                  )
#                                  
#                                  # Add a column for the mortality regime name
#                                  result_df$mortality <- regime$name
#                                  
#                                  # Return the data frame for this iteration
#                                  return(result_df)
#                                }
#  


# ----------------------------------------------------------
#### If you have run the sweep in a previous work session and would like to read in the data for plotting now, run the lines below: 
sweep_results_rmr1 <- foreach(regime = mortality_regimes,
                              .packages = c("dplyr", "data.table", "jsonlite"),
                              .combine = rbind) %dopar% {

                                # Run sweep for this mortality regime
                                result_df <- read_sweep(
                                  mortality = regime,
                                  lesion_rates = lesion_exposure_rates,
                                  reps = reps,
                                  rmr = 1 # risk-free lesions
                                )

                                # Add a column for the mortality regime name
                                result_df$mortality <- regime$name

                                # Return the data frame for this iteration
                                return(result_df)
                              }

sweep_results_rmr2 <- foreach(regime = mortality_regimes,
                              .packages = c("dplyr", "data.table", "jsonlite"),
                              .combine = rbind) %dopar% {

                                # Run sweep for this mortality regime
                                result_df <- read_sweep(
                                  mortality = regime,
                                  lesion_rates = lesion_exposure_rates,
                                  reps = reps,
                                  rmr = 2 # risk-free lesions
                                )

                                # Add a column for the mortality regime name
                                result_df$mortality <- regime$name

                                # Return the data frame for this iteration
                                return(result_df)
                              }
# ----------------------------------------------------------


# label the relative_mortality_risk (rmr) in each data frame
sweep_results_rmr1$rmr <- 1
sweep_results_rmr2$rmr <- 2

# Combine all results into a single data frame
sweep_data <- rbind(sweep_results_rmr1, sweep_results_rmr2) %>%
  # recode mortality as a factor and specify the order of the levels from fast to slow mortality schedules. 
  mutate(mortality = factor(mortality, levels = c("CoaleDemenyWestF3", "CoaleDemenyWestF5", "CoaleDemenyWestF11", "CoaleDemenyWestF15", "CoaleDemenyWestF17", "CoaleDemenyWestF21")),
         mortality = recode(mortality,
                            "CoaleDemenyWestF3" = "CDW Level 3",
                            "CoaleDemenyWestF5" = "CDW Level 5",
                            "CoaleDemenyWestF11" = "CDW Level 11",
                            "CoaleDemenyWestF15" = "CDW Level 15",
                            "CoaleDemenyWestF17" = "CDW Level 17",
                            "CoaleDemenyWestF21" = "CDW Level 21"))




#### COMBINE & PLOT ####

# Convert to summary data
plot_data <- lesions_to_percents(sweep_data, group_vars = c("mortality", "lesion_formation_rate", "rmr", "rep")) 



# ---------------
# Figure S1
# ---------------
# Line plot showing the effect of varying lesion formation rate and mortality regime on age-specific lesion frequency
lineplot1 <- ggplot(data = plot_data %>% filter(rmr == 1), aes(x = interval_midpoint, y = Lesion_Percent)) +
  geom_line(aes(group = interaction(rep, lesion_formation_rate), 
                color = as.factor(lesion_formation_rate)), 
            alpha = 0.05) +
  facet_wrap(~ mortality, ncol = 2) +  # One panel per mortality regime
  labs(
    title = "Lesion Percent by Age, Colored by Exposure Rate",
    x = "Age (years)",
    y = "Percent with Lesions",
    color = "Exposure Rate"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.position = "none")



lineplot2 <- ggplot(data = plot_data %>% filter(rmr == 2), aes(x = interval_midpoint, y = Lesion_Percent)) +
  geom_line(aes(group = interaction(rep, lesion_formation_rate), 
                color = as.factor(lesion_formation_rate)), 
            alpha = 0.05) +
  facet_wrap(~ mortality, ncol = 2) +  # One panel per mortality regime
  labs(
    x = "Age (years)",
    y = "Percent with Lesions",
    color = "Exposure Rate",
    caption = "Each line = one simulation replicate"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.position = "bottom")

legend <- get_legend(lineplot2)
lineplot2 <- lineplot2 +
  theme(legend.position = "none")


ggsave(filename =  paste0(here("Known Unknowns", "Figures"), "/Figure_S1.pdf"), 
       plot = plot_grid(lineplot1,
                        lineplot2,
                        legend,
                        ncol = 1,
                        rel_heights = c(1,1, 0.1),
                        labels = c("A", "B")),
       width = 7,
       height = 10,
       units = "in",
       dpi = 1200)




# ---------------
# Figure S2
# ---------------
# Same as main Figure 4, but without truncated x axis. 

# See `Main Analysis.Rmd`





# ---------------
# Figure S3
# ---------------
# Combined plot of all worlds (same as below, just superimposed on the same grid): Figure S3
all_lineplot <- ggplot(data = plot_data, aes(x = interval_midpoint, y = Lesion_Percent)) +
  geom_line(aes(group = interaction(rmr, rep), color = as.factor(rmr)), alpha = 0.05) +
  scale_color_manual(values = c("black", "magenta")) +
  facet_grid(
    rows = vars(mortality), 
    cols = vars(lesion_formation_rate), 
    labeller = labeller(
      lesion_formation_rate = function(x) paste("Exposure Rate:", x)
    )
  ) +
  labs(
    title = "Lesion Percent by Age, Across Mortality Regimes and Exposure Rates",
    x = "Age (years)",
    y = "Percent with Lesions",
    caption = "Each line = one simulation replicate",
    color = "Relative Mortality Risk"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 8.5),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

# image size: full page, landscape, with 1/2 inch for image caption text
ggsave(filename =  paste0(here("Known Unknowns", "Figures"), "/Figure_S3.pdf"), 
       plot = all_lineplot,
       width = 13.5,
       height = 8,
       units = "in",
       dpi = 1200)




# Generate plots of lesion frequency in each run of the model

# ---------------
# Figure S4
# ---------------
#...when lesions are risk-free.
full_lineplot1 <- sweep_lineplot(plot_data %>% filter(rmr == 1), plot_color = "black") 

ggsave(filename = paste0(here("Known Unknowns", "Figures"), "/Figure_S4.pdf"), 
       plot = full_lineplot1,
       width = 13.5,
       height = 8,
       units = "in",
       dpi = 1200)



# ---------------
# Figure S5
# ---------------
#...when lesions double your risk. 
full_lineplot2 <- sweep_lineplot(plot_data %>% filter(rmr == 2), plot_color = "magenta") 

ggsave(filename =  paste0(here("Known Unknowns", "Figures"), "/Figure_S5.pdf"), 
       plot = full_lineplot2,
       width = 13.5,
       height = 8,
       units = "in",
       dpi = 1200)

###############################################################################################################


# -------------------------
#### SURVIVAL ANALYSES ####
# -------------------------

# See Sweep_Utility_Functions file for details of the run_survival_analysis() function.
# See Plotting_Functions file for details of the plot_survival_curves() function. 

# Survival Analysis: First Run of the risk-free model

survival1 <- run_survival_analysis(sweep_data %>% filter(rmr == 1), parallel = TRUE) 
# Survival Analysis: First Run of the risk-doubling model
survival2 <- run_survival_analysis(sweep_data %>% filter(rmr == 2), parallel = TRUE)


# generate faceted plots of survival curves for one run of each combination of parameter values. 



# ---------------
# Figure S6
# ---------------
rmr1_survival_plots <- plot_survival_curves(survival1$survival_data %>%
                                              filter(rep == 1))

ggsave(filename =  paste0(here("Known Unknowns", "Figures"), "/Figure_S6.pdf"), 
       plot = rmr1_survival_plots,
       width = 15,
       height = 8,
       units = "in",
       dpi = 1200)



# ---------------
# Figure S7
# ---------------
rmr2_survival_plots <- plot_survival_curves(survival2$survival_data%>%
                                              filter(rep == 2))

ggsave(filename =  paste0(here("Known Unknowns", "Figures"), "/Figure_S7.pdf"), 
       plot = rmr2_survival_plots,
       width = 15,
       height = 8,
       units = "in",
       dpi = 1200)




# ---------------
# Figure S8
# ---------------

# summarise The % of significant logrank tests for risk-free lesion cohorts
logrank_rmr1 <- survival1$logrank_results %>%
  mutate(significant = if_else(p_value < 0.05, 1, 0)) %>%
  group_by(mortality, lesion_formation_rate) %>%
  summarise(percent_significant = (sum(significant, na.rm = T) / reps) * 100, .groups = "keep")  

# summarise significant logrank tests for risk-DOUBLING lesion cohorts
logrank_rmr2 <- survival2$logrank_results %>%
  mutate(significant = if_else(p_value < 0.05, 1, 0)) %>%
  group_by(mortality, lesion_formation_rate) %>%
  summarise(percent_significant = (sum(significant, na.rm = T) / reps) * 100, .groups = "keep") 


# plot the results
logrank_rmr1_plot <- ggplot(logrank_rmr1, aes(x = lesion_formation_rate, y = percent_significant)) +
  geom_line(aes(group = mortality, color = mortality)) +
  labs(x = "Exposure Rate", y = "% of runs with Significant Log-Rank test", color = "Mortality") +
  theme_bw() +
  scale_x_continuous(limits = c(0.005,.1), expand = c(0,0), 
                     breaks = c(0.01, 0.03, 0.05, 0.07, 0.09)) +
  theme(legend.position = "none")

logrank_rmr2_plot <- ggplot(logrank_rmr2, aes(x = lesion_formation_rate, y = percent_significant)) +
  geom_line(aes(group = mortality, color = mortality)) +
  labs(x = "Exposure Rate", y = "% of runs with Significant Log-Rank test", color = "Mortality") +
  theme_bw() +
  scale_x_continuous(limits = c(0.005,.1), expand = c(0,0), 
                     breaks = c(0.01, 0.03, 0.05, 0.07, 0.09)) +
  theme(legend.position = "bottom")

legend <- get_legend(logrank_rmr2_plot)
logrank_rmr2_plot <- logrank_rmr2_plot +
  theme(legend.position = "none")

logrank_significance_plot <- plot_grid(logrank_rmr1_plot, 
          logrank_rmr2_plot, 
          legend, 
          ncol = 1,
          rel_heights = c(1,1, 0.2)) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))


ggsave(filename = paste0(here("Known Unknowns", "Figures"), "/Figure_S8.pdf"), 
       plot = logrank_significance_plot,
       width = 5,
       height = 7,
       units = "in",
       dpi = 1200)




###############################################################################################################



# Summary table of Mortality Schedules (Table S1)

life_expectancies <- sweep_data %>%
  filter(rmr == 1) %>%
  mutate(adult = if_else(Age > 15, 1, 0)) %>%
  group_by(mortality) %>%
  summarise(life_expectancy_at_birth = round(mean(Age), 1),
            juvenile_mortality = round(sum(Age < 15) / n() * 100, 1),
            adult_life_expectancy = round(mean(Age[adult == 1]))) 

names(life_expectancies) <- c("Mortality regime", "Life expectancy at birth", "Juvenile mortality (%)", "Mean adult age at death")


life_expectancies

# ---------------
# Table S1
# ---------------
# write.csv(life_expectancies, file = paste0(here("Known Unknowns", "Tables"), "/Table_S1.csv"))








