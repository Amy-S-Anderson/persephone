


# --------------------------------------------------------------------------------------------------------------
#### PLOTTING FUNCTIONS ####
# --------------------------------------------------------------------------------------------------------------

# These functions are called in the Sensitivity_Analysis.R file to generate some of the plots in the supplemental materials. 




# -------------------------
# Function to calculate percent of individuals at each age in the cemetery with skeletal lesions
# -------------------------

lesions_to_percents <- function(cemetery_df, group_vars) { # this function takes a data frame and a c() of named variables. 
  
  # Define age interval labels and their midpoints
  age_lookup <- data.frame(
    Age_Interval = factor(
      c("0-1", "2-5", "6-9", "10-14", "15-19",
        "20-29", "30-39", "40-49", "50-59", "60+"),
      levels = c("0-1", "2-5", "6-9", "10-14", "15-19",
                 "20-29", "30-39", "40-49", "50-59", "60+")
    ),
    interval_midpoint = c(0.5, 3.5, 7.5, 12.5, 17.5,
                          25, 35, 45, 55, 75)
  )
  
  # Create age intervals and recode mortality levels
  plot_prep_data <- cemetery_df %>%
    mutate(
      Age_Interval = factor(case_when(
        Age < 2 ~ "0-1",
        Age >= 2 & Age < 6 ~ "2-5",
        Age >= 6 & Age < 10 ~ "6-9",
        Age >= 10 & Age < 15 ~ "10-14",
        Age >= 15 & Age < 20 ~ "15-19",
        Age >= 20 & Age < 30 ~ "20-29",
        Age >= 30 & Age < 40 ~ "30-39",
        Age >= 40 & Age < 50 ~ "40-49",
        Age >= 50 & Age < 60 ~ "50-59",
        Age >= 60 ~ "60+"
      ), levels = levels(age_lookup$Age_Interval))
    ) %>%
    # Group by specified variables and Age_Interval
    group_by(across(all_of(c(group_vars, "Age_Interval")))) %>%
    summarise(
      Lesion_Percent = sum(Lesion) / n() * 100,
      .groups = "drop"
    ) %>%
    # Join midpoints only where Age_Interval exists
    left_join(age_lookup, by = "Age_Interval")
  
  return(plot_prep_data)
}






# -------------------------
# Make a line plot for this lesion_exposure parameter sweep:
# -------------------------

sweep_lineplot <- function(plot_data, plot_color){
  ggplot(data = plot_data, aes(x = interval_midpoint, y = Lesion_Percent)) +
    geom_line(aes(group = rep), color = plot_color, alpha = 0.05) +
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
      caption = "Each line = one simulation replicate"
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 8.5),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}





# -------------------------
# Plot AVERAGE OUTCOMES for each combination scenario of lesion exposure and mortality
# -------------------------

scenario_lineplot <- function(plot_summary_data){
  ggplot(data = plot_summary_data, aes(x = Age_Interval, y = Lesion_Percent)) +
    geom_line(aes( group = lesion_formation_rate,
                   color = as.factor(lesion_formation_rate))) +
    facet_wrap(~ mortality, ncol = 2) +  # One panel per mortality regime
    labs(
      title = "Lesion Percent by Age, Colored by Exposure Rate",
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
    guides(colour = guide_legend(override.aes = list(alpha = 1)))
  
}





# -------------------------
# Plot Survival curves for a single run of each parameter combination
# -------------------------

plot_survival_curves <- function(survival_data){
  lesion_colors = c("black", '#b9c28d')
  plots <- ggplot(survival_data, aes(x = time, y = survival)) +
    geom_line(aes(group = group, color = group)) +
    scale_color_manual(values = lesion_colors, labels = c("Absent", "Present")) +
    facet_grid(mortality ~ lesion_formation_rate,
               labeller = labeller(
                 lesion_formation_rate = function(x) paste("Exposure Rate:", x))) +
    labs(color = "Skeletal Lesions") +
    theme_bw()
  
  return(plots)
}




