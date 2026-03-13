#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom stats pchisq runif
#' @importFrom utils write.csv
"_PACKAGE"

# Suppress R CMD check NOTEs for NSE column references used in dplyr/ggplot2
utils::globalVariables(c(
  "Age_Interval", "CoaleDemenyWestF5", "Lesion_Percent",
  "d", "group", "interval_midpoint", "lesion",
  "lesion_formation_rate", "mortality", "rep", "survival", "time"
))
