
# -------------------------
##### Function to Run a Survival Analysis on each run of the model ####
# This function is used in both the main analysis and the sensitivity analysis.
# It returns a list of two items:
#    1) the survival data for plotting a survival curve for each run.
#    2) the results of a log-rank test for each run.
# it uses parallel processing for a faster run time, since it is designed to run survival analysis on a huge number of data sets.
# -------------------------

#' @importFrom foreach foreach %dopar%
#' @export
run_survival_analysis <- function(sweep_data, parallel = TRUE, workers = NULL) {
  # Ensure dead column exists (event indicator for survival analysis — all are dead)
  sweep_data$dead <- 1

  # Split data into groups once (avoids repeated filter calls in loop)
  grouped_data <- sweep_data %>%
    group_by(mortality, lesion_formation_rate, rep) %>%
    group_split()

  # Set up parallel backend
  if (parallel) {
    if (is.null(workers)) workers <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(workers)
    doParallel::registerDoParallel(cl)
  }

  # Parallel loop
  results <- foreach(d = grouped_data,
                     .packages = c("dplyr", "survival")) %dopar% {

                       combo <- d[1, c("mortality", "lesion_formation_rate", "rep"), drop = FALSE]

                       # Skip if columns are missing
                       if (!all(c("age", "dead", "lesion") %in% names(d))) {
                         return(NULL)
                       }

                       # Survival object + model
                       surv_obj <- survival::Surv(time = d$age, event = d$dead)
                       fit <- survival::survfit(surv_obj ~ lesion, data = d)

                       surv_summary <- summary(fit)

                       surv_df <- data.frame(
                         time = surv_summary$time,
                         survival = surv_summary$surv,
                         group = as.character(surv_summary$strata),
                         mortality = combo$mortality,
                         lesion_formation_rate = combo$lesion_formation_rate,
                         rep = combo$rep
                       )

                       # Log-rank test
                       logrank_test <- survival::survdiff(surv_obj ~ lesion, data = d)
                       p_value <- 1 - pchisq(logrank_test$chisq, df = length(logrank_test$n) - 1)

                       list(
                         survival_data = surv_df,
                         logrank_results = data.frame(
                           mortality = combo$mortality,
                           lesion_formation_rate = combo$lesion_formation_rate,
                           rep = combo$rep,
                           p_value = p_value
                         )
                       )
                     }

  # Stop parallel backend
  if (parallel) parallel::stopCluster(cl)

  # Clean results
  results <- Filter(Negate(is.null), results)

  # Bind results
  survival_data <- bind_rows(lapply(results, `[[`, "survival_data"))
  logrank_results <- bind_rows(lapply(results, `[[`, "logrank_results"))

  return(list(
    survival_data = survival_data,
    logrank_results = logrank_results
  ))
}
