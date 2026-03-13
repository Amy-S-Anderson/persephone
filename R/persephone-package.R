#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom stats pchisq runif
#' @importFrom utils write.csv
"_PACKAGE"

#' Launch the Persephone Shiny app
#'
#' @export
run_app <- function() {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required to run the app. Install it with install.packages('shiny').")
  }
  shiny::runApp(system.file("app", package = "persephone"))
}

# Suppress R CMD check NOTEs for NSE column references used in dplyr/ggplot2
utils::globalVariables(c(
  "Age_Interval", "CoaleDemenyWestF5", "Lesion_Percent",
  "d", "group", "interval_midpoint", "lesion",
  "lesion_formation_rate", "mortality", "rep", "survival", "time"
))
