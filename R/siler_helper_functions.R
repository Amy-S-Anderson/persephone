


# helper functions for siler models

baseline_death_prob <- function(age, regime) {
  regime$a1 * exp(-regime$b1 * age) +
    regime$a2 +
    regime$a3 * exp(regime$b3 * age)
}

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


