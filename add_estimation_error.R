
 
 #' Add structured noise to known ages
 #'
 #' @param ages Numeric vector of true ages (0-100)
 #' @param sd_per_year Rate at which SD increases per year after age 20 (default calibrated so SD ≈ 7.5 at age 60)
 #' @param sd_at_20 SD of random error at age 20 (default = 2, giving ~±4 yr range)
 #' @param bias_start Age at which systematic bias begins (default = 60)
 #' @param bias_at_90 Mean systematic error at age 90 (default = 25)
 #' @param seed Optional random seed for reproducibility
 #' @return Numeric vector of ages with noise added
 
 add_age_noise <- function(ages,
                           sd_per_year  = 0.1375,
                           sd_at_20     = 2,
                           bias_start   = 60,
                           bias_at_90   = 25,
                           seed         = NULL) {
   
   if (!is.null(seed)) set.seed(seed)
   
   # --- Random error: SD is piecewise linear in age ---
   # SD = 0 at age 0, rises linearly to sd_at_20 at age 20,
   # then continues rising at sd_per_year beyond age 20.
   sd_vec <- case_when(
     ages <= 0  ~ 0,
     ages <= 20 ~ (ages / 20) * sd_at_20,
     TRUE       ~ sd_at_20 + (ages - 20) * sd_per_year
   )
   
   random_error <- rnorm(length(ages), mean = 0, sd = sd_vec)
   
   # --- Systematic error: linear ramp starting at bias_start ---
   # Bias = 0 up to bias_start, then rises linearly.
   # Calibrated so bias = bias_at_90 at age 90.
   bias_slope <- bias_at_90 / (90 - bias_start)
   
   systematic_error <- case_when(
     ages <= bias_start ~ 0,
     TRUE               ~ (ages - bias_start) * bias_slope
   )
   
   ages + random_error - systematic_error
 } 

 