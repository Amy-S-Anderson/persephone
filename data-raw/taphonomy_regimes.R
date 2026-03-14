# Each regime is a one-row data frame with columns: a1, b1, a2, a3, b3, name

no_decay <- data.frame(
  a1 = 0,
  b1 = 0,
  a2 = 0,
  a3 = 0,
  b3 = 0,
  name = "no_decay")

weak_decay <- data.frame(
  a1 = 20 * 0.05,
  b1 = 0.29,
  a2 = 0.8 * 0.05,
  a3 = 0.006 * 0.05,
  b3 = 0.05,
  name = "weak_decay")

moderate_decay <- data.frame(
  a1 = 20 * 0.1,
  b1 = 0.29,
  a2 = 0.8 * 0.1,
  a3 = 0.006 * 0.1,
  b3 = 0.05,
  name = "moderate_decay")

strong_decay <- data.frame(
  a1 = 20 * 0.5,
  b1 = 0.29,
  a2 = 0.8 * 0.5,
  a3 = 0.006 * 0.5,
  b3 = 0.05,
  name = "strong_decay")

usethis::use_data(
  no_decay,
  weak_decay,
  moderate_decay,
  strong_decay,
  overwrite = TRUE
)
