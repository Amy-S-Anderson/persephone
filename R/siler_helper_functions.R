

# helper functions for siler models

#' Compute baseline Siler death probability
#' @param age Numeric. Age in years.
#' @param regime Data frame with Siler parameters (a1, b1, a2, a3, b3).
#' @return Numeric death probability.
#' @export
baseline_death_prob <- function(age, regime) {
  regime$a1 * exp(-regime$b1 * age) +
    regime$a2 +
    regime$a3 * exp(regime$b3 * age)
}

#' Compute lesion-modified Siler death probability
#' @param age Numeric. Age in years.
#' @param regime Data frame with Siler parameters (a1, b1, a2, a3, b3).
#' @param risk_type Character. One of "proportional", "time_decreasing", "time_increasing".
#' @param rmr Numeric. Relative mortality risk multiplier.
#' @return Numeric death probability.
#' @export
lesion_death_prob <- function(age, regime, risk_type, rmr) {
  base <- baseline_death_prob(age, regime)

  if (risk_type == "proportional") {
    out <- base * rmr
  } else if (risk_type == "time_decreasing") {
    out <- base * rmr / ((age / 10) + rmr)
  } else if (risk_type == "time_increasing") {
    out <- base * ((age / 10) + rmr) / rmr
  } else {
    stop("Unknown risk_type")
  }

  out
}

#' Compute discrete age-at-death distribution from a probability function
#' @param prob_fun Function taking age and returning death probability.
#' @param age_max Integer. Maximum age to compute. Default 100.
#' @return Data frame with columns age, qx, px, S_start, dx.
#' @export
make_discrete_death_distribution <- function(prob_fun, age_max = 100) {
  ages <- 1:age_max
  qx <- sapply(ages, prob_fun)

  # enforce valid probabilities
  qx <- pmin(pmax(qx, 0), 1)

  px <- 1 - qx

  # survival to start of each age
  # S_start[1] = 1
  S_start <- c(1, cumprod(px))[1:length(ages)]

  # probability of death at each age
  dx <- S_start * qx

  data.frame(
    age = ages,
    qx = qx,
    px = px,
    S_start = S_start,
    dx = dx
  )
}


