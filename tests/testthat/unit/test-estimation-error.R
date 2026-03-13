# Unit tests for estimation error module

# -----------------------------------------------------------------------------
# Test fixtures
# -----------------------------------------------------------------------------

# Default error parameters: c(sd_per_year, sd_at_20, bias_start, bias_at_90)
error_param_default <- c(0.1375, 2, 60, 25)

# Zero-error parameters (no random error, no bias)
error_param_zero <- c(0, 0, 100, 0)

create_test_state <- function(n = 100) {
  data.frame(
    agent_id = 1:n,
    age = as.numeric(sample(0:80, n, replace = TRUE)),
    lesion = sample(c(TRUE, FALSE), n, replace = TRUE),
    dead = rep(TRUE, n),
    was_deposited = rep(TRUE, n),
    in_sample = rep(TRUE, n)
  )
}

# -----------------------------------------------------------------------------
# Input validation tests
# -----------------------------------------------------------------------------

test_that("estimation error module raises error when error_model is NULL", {

  state <- create_test_state()

  expect_error(
    apply_estimation_error(
      state = state,
      error_model = NULL,
      error_param = error_param_default
    )
  )
})

test_that("estimation error module raises error when error_model is unsupported", {
  state <- create_test_state()

  expect_error(
    apply_estimation_error(
      state = state,
      error_model = "unknown",
      error_param = error_param_default
    )
  )

  expect_error(
    apply_estimation_error(
      state = state,
      error_model = "constant",
      error_param = error_param_default
    )
  )
})

test_that("bespoke_increasing_sd model requires length-4 parameter vector", {
  state <- create_test_state()

  expect_error(
    apply_estimation_error(
      state = state,
      error_model = "bespoke_increasing_sd",
      error_param = c(0.1, 2, 60)
    )
  )

  expect_error(
    apply_estimation_error(
      state = state,
      error_model = "bespoke_increasing_sd",
      error_param = c(0.1, 2, 60, 25, 10)
    )
  )

  expect_error(
    apply_estimation_error(
      state = state,
      error_model = "bespoke_increasing_sd",
      error_param = c()
    )
  )
})

# -----------------------------------------------------------------------------
# Successful execution tests
# -----------------------------------------------------------------------------

test_that("bespoke_increasing_sd model executes successfully", {
  state <- create_test_state()

  expect_no_error(
    apply_estimation_error(
      state = state,
      error_model = "bespoke_increasing_sd",
      error_param = error_param_default
    )
  )
})

test_that("estimation error module returns a dataframe", {
  state <- create_test_state()

  result <- apply_estimation_error(
    state = state,
    error_model = "bespoke_increasing_sd",
    error_param = error_param_default
  )

  expect_s3_class(result, "data.frame")
})

# -----------------------------------------------------------------------------
# Output structure tests
# -----------------------------------------------------------------------------

test_that("output state contains estimated_age column", {
  state <- create_test_state()

  result <- apply_estimation_error(
    state = state,
    error_model = "bespoke_increasing_sd",
    error_param = error_param_default
  )

  expect_true("estimated_age" %in% names(result))
})

test_that("output state contains all original columns plus estimated_age", {
  state <- create_test_state()

  result <- apply_estimation_error(
    state = state,
    error_model = "bespoke_increasing_sd",
    error_param = error_param_default
  )

  # All original columns should be present
  expect_true(all(names(state) %in% names(result)))
  # Plus estimated_age

  expect_equal(ncol(result), ncol(state) + 1)
})

test_that("output state has same number of rows as input state", {
  state <- create_test_state(n = 150)

  result <- apply_estimation_error(
    state = state,
    error_model = "bespoke_increasing_sd",
    error_param = error_param_default
  )

  expect_equal(nrow(result), nrow(state))
})

test_that("estimated_age column is numeric", {
  state <- create_test_state()

  result <- apply_estimation_error(
    state = state,
    error_model = "bespoke_increasing_sd",
    error_param = error_param_default
  )

  expect_type(result$estimated_age, "double")
})

# -----------------------------------------------------------------------------
# Eligibility tests (who gets estimated_age calculated)
# -----------------------------------------------------------------------------

test_that("agents in_sample TRUE get estimated_age calculated (not NA)", {
  state <- data.frame(
    agent_id = 1:10,
    age = as.numeric(rep(30, 10)),
    lesion = rep(FALSE, 10),
    dead = rep(TRUE, 10),
    was_deposited = rep(TRUE, 10),
    in_sample = rep(TRUE, 10)
  )

  result <- apply_estimation_error(
    state = state,
    error_model = "bespoke_increasing_sd",
    error_param = error_param_default
  )

  # All agents are in_sample, so all should have estimated_age values

  expect_true(all(!is.na(result$estimated_age)))
})

test_that("agents in_sample FALSE get estimated_age NA", {
  state <- data.frame(
    agent_id = 1:10,
    age = as.numeric(rep(30, 10)),
    lesion = rep(FALSE, 10),
    dead = rep(TRUE, 10),
    was_deposited = rep(TRUE, 10),
    in_sample = rep(FALSE, 10)
  )

  result <- apply_estimation_error(
    state = state,
    error_model = "bespoke_increasing_sd",
    error_param = error_param_default
  )

  # No agents are in_sample, so all should have NA

  expect_true(all(is.na(result$estimated_age)))
})

test_that("mixed in_sample population: only in_sample TRUE get values", {
  state <- data.frame(
    agent_id = 1:10,
    age = as.numeric(rep(30, 10)),
    lesion = rep(FALSE, 10),
    dead = rep(TRUE, 10),
    was_deposited = rep(TRUE, 10),
    in_sample = c(rep(TRUE, 5), rep(FALSE, 5))
  )

  result <- apply_estimation_error(
    state = state,
    error_model = "bespoke_increasing_sd",
    error_param = error_param_default
  )

  # First 5 should have values, last 5 should be NA
  expect_true(all(!is.na(result$estimated_age[1:5])))
  expect_true(all(is.na(result$estimated_age[6:10])))
})

# -----------------------------------------------------------------------------
# State preservation tests (columns that should NOT change)
# -----------------------------------------------------------------------------

test_that("agent ages are not modified by estimation error module", {
  state <- create_test_state(n = 100)

  result <- apply_estimation_error(
    state = state,
    error_model = "bespoke_increasing_sd",
    error_param = error_param_default
  )

  expect_equal(result$age, state$age)
})

test_that("lesion status is not modified by estimation error module", {
  state <- create_test_state(n = 100)

  result <- apply_estimation_error(
    state = state,
    error_model = "bespoke_increasing_sd",
    error_param = error_param_default
  )

  expect_equal(result$lesion, state$lesion)
})

test_that("dead status is not modified by estimation error module", {
  state <- create_test_state(n = 100)

  result <- apply_estimation_error(
    state = state,
    error_model = "bespoke_increasing_sd",
    error_param = error_param_default
  )

  expect_equal(result$dead, state$dead)
})

test_that("was_deposited is not modified by estimation error module", {
  state <- create_test_state(n = 100)
  state$was_deposited <- sample(c(TRUE, FALSE), 100, replace = TRUE)

  result <- apply_estimation_error(
    state = state,
    error_model = "bespoke_increasing_sd",
    error_param = error_param_default
  )

  expect_equal(result$was_deposited, state$was_deposited)
})

test_that("in_sample is not modified by estimation error module", {
  state <- create_test_state(n = 100)
  state$in_sample <- sample(c(TRUE, FALSE), 100, replace = TRUE)

  result <- apply_estimation_error(
    state = state,
    error_model = "bespoke_increasing_sd",
    error_param = error_param_default
  )

  expect_equal(result$in_sample, state$in_sample)
})

test_that("agent_id is not modified by estimation error module", {
  state <- create_test_state(n = 100)

  result <- apply_estimation_error(
    state = state,
    error_model = "bespoke_increasing_sd",
    error_param = error_param_default
  )

  expect_equal(result$agent_id, state$agent_id)
})

# -----------------------------------------------------------------------------
# Error model logic tests
# -----------------------------------------------------------------------------

test_that("estimated_age equals true age when all error parameters are zero", {
  state <- data.frame(
    agent_id = 1:100,
    age = as.numeric(rep(30, 100)),
    lesion = rep(FALSE, 100),
    dead = rep(TRUE, 100),
    was_deposited = rep(TRUE, 100),
    in_sample = rep(TRUE, 100)
  )

  # Zero SD parameters and bias_start beyond any age means no error
  zero_param <- c(0, 0, 200, 0)  # sd_per_year=0, sd_at_20=0, bias_start=200, bias_at_90=0

  result <- apply_estimation_error(
    state = state,
    error_model = "bespoke_increasing_sd",
    error_param = zero_param
  )

  # With zero error, estimated_age should equal true age

  expect_equal(result$estimated_age, state$age)
})

test_that("systematic bias causes underestimation for ages above bias_start", {
  set.seed(12345)
  state <- data.frame(
    agent_id = 1:1000,
    age = as.numeric(rep(80, 1000)),  # Well above bias_start of 60
    lesion = rep(FALSE, 1000),
    dead = rep(TRUE, 1000),
    was_deposited = rep(TRUE, 1000),
    in_sample = rep(TRUE, 1000)
  )

  # Use zero random error but keep systematic bias
  # bias_at_90 = 25 means at age 90, bias = 25 years
  # At age 80: bias = (80 - 60) * (25 / (90 - 60)) = 20 * (25/30) ≈ 16.67
  bias_only_param <- c(0, 0, 60, 25)

  result <- apply_estimation_error(
    state = state,
    error_model = "bespoke_increasing_sd",
    error_param = bias_only_param
  )

  # All estimated ages should be less than true ages due to systematic bias
  expect_true(all(result$estimated_age < state$age))

  # Mean underestimation should be approximately 16.67 years
  mean_error <- mean(state$age - result$estimated_age)
  expect_equal(mean_error, 20 * (25 / 30), tolerance = 0.01)
})

test_that("no systematic bias for ages at or below bias_start", {
  set.seed(12345)
  state <- data.frame(
    agent_id = 1:100,
    age = as.numeric(rep(50, 100)),  # Below bias_start of 60
    lesion = rep(FALSE, 100),
    dead = rep(TRUE, 100),
    was_deposited = rep(TRUE, 100),
    in_sample = rep(TRUE, 100)
  )

  # Use zero random error, bias_start = 60
  no_random_param <- c(0, 0, 60, 25)

  result <- apply_estimation_error(
    state = state,
    error_model = "bespoke_increasing_sd",
    error_param = no_random_param
  )

  # With no random error and age below bias_start, estimated should equal true

  expect_equal(result$estimated_age, state$age)
})

test_that("random error SD is zero at age zero", {
  set.seed(12345)
  state <- data.frame(
    agent_id = 1:100,
    age = as.numeric(rep(0, 100)),
    lesion = rep(FALSE, 100),
    dead = rep(TRUE, 100),
    was_deposited = rep(TRUE, 100),
    in_sample = rep(TRUE, 100)
  )

  # Standard parameters but age = 0 means SD = 0
  result <- apply_estimation_error(
    state = state,
    error_model = "bespoke_increasing_sd",
    error_param = error_param_default
  )

  # With SD = 0 at age 0 and no bias (age < bias_start), all should equal 0

  expect_equal(result$estimated_age, state$age)
})

test_that("random error variance increases with age", {
  set.seed(12345)

  # Young agents (age 10) - lower SD expected
  state_young <- data.frame(
    agent_id = 1:1000,
    age = as.numeric(rep(10, 1000)),
    lesion = rep(FALSE, 1000),
    dead = rep(TRUE, 1000),
    was_deposited = rep(TRUE, 1000),
    in_sample = rep(TRUE, 1000)
  )

  # Older agents (age 40) - higher SD expected
  state_old <- data.frame(
    agent_id = 1:1000,
    age = as.numeric(rep(40, 1000)),
    lesion = rep(FALSE, 1000),
    dead = rep(TRUE, 1000),
    was_deposited = rep(TRUE, 1000),
    in_sample = rep(TRUE, 1000)
  )

  # Use parameters with random error but no bias
  random_only_param <- c(0.1375, 2, 200, 0)  # bias_start very high

  set.seed(12345)
  result_young <- apply_estimation_error(
    state = state_young,
    error_model = "bespoke_increasing_sd",
    error_param = random_only_param
  )

  set.seed(12345)
  result_old <- apply_estimation_error(
    state = state_old,
    error_model = "bespoke_increasing_sd",
    error_param = random_only_param
  )

  # Variance of errors should be higher for older agents
  var_young <- var(result_young$estimated_age - state_young$age)
  var_old <- var(result_old$estimated_age - state_old$age)

  expect_gt(var_old, var_young)
})

# -----------------------------------------------------------------------------
# Edge case tests
# -----------------------------------------------------------------------------

test_that("handles empty state (zero rows)", {
  state <- data.frame(
    agent_id = integer(0),
    age = numeric(0),
    lesion = logical(0),
    dead = logical(0),
    was_deposited = logical(0),
    in_sample = logical(0)
  )

  result <- apply_estimation_error(
    state = state,
    error_model = "bespoke_increasing_sd",
    error_param = error_param_default
  )

  expect_equal(nrow(result), 0)
  expect_true("estimated_age" %in% names(result))
})

test_that("handles state with no eligible agents (all in_sample FALSE)", {
  state <- data.frame(
    agent_id = 1:10,
    age = as.numeric(rep(30, 10)),
    lesion = rep(FALSE, 10),
    dead = rep(TRUE, 10),
    was_deposited = rep(TRUE, 10),
    in_sample = rep(FALSE, 10)
  )

  expect_no_error(
    result <- apply_estimation_error(
      state = state,
      error_model = "bespoke_increasing_sd",
      error_param = error_param_default
    )
  )

  expect_true(all(is.na(result$estimated_age)))
})

test_that("handles agents at extreme ages (e.g., 100+)", {
  state <- data.frame(
    agent_id = 1:10,
    age = as.numeric(c(100, 105, 110, 115, 120, 90, 95, 85, 80, 75)),
    lesion = rep(FALSE, 10),
    dead = rep(TRUE, 10),
    was_deposited = rep(TRUE, 10),
    in_sample = rep(TRUE, 10)
  )

  expect_no_error(
    result <- apply_estimation_error(
      state = state,
      error_model = "bespoke_increasing_sd",
      error_param = error_param_default
    )
  )

  # All should have estimated ages (not NA)
  expect_true(all(!is.na(result$estimated_age)))
})

test_that("handles negative ages gracefully", {
  # Edge case: what if age is negative (shouldn't happen, but defensive)
  state <- data.frame(
    agent_id = 1:5,
    age = as.numeric(c(-1, 0, 1, 2, 3)),
    lesion = rep(FALSE, 5),
    dead = rep(TRUE, 5),
    was_deposited = rep(TRUE, 5),
    in_sample = rep(TRUE, 5)
  )

  expect_no_error(
    result <- apply_estimation_error(
      state = state,
      error_model = "bespoke_increasing_sd",
      error_param = error_param_default
    )
  )

  # Should return values (not error out)
  expect_true(all(!is.na(result$estimated_age)))
})

# -----------------------------------------------------------------------------
# Reproducibility tests
# -----------------------------------------------------------------------------

test_that("results are reproducible with set.seed before call", {
  state <- create_test_state(n = 100)

  set.seed(42)
  result1 <- apply_estimation_error(
    state = state,
    error_model = "bespoke_increasing_sd",
    error_param = error_param_default
  )

  set.seed(42)
  result2 <- apply_estimation_error(
    state = state,
    error_model = "bespoke_increasing_sd",
    error_param = error_param_default
  )

  expect_equal(result1$estimated_age, result2$estimated_age)
})

test_that("different seeds produce different results", {
  state <- data.frame(
    agent_id = 1:100,
    age = as.numeric(rep(30, 100)),
    lesion = rep(FALSE, 100),
    dead = rep(TRUE, 100),
    was_deposited = rep(TRUE, 100),
    in_sample = rep(TRUE, 100)
  )

  set.seed(42)
  result1 <- apply_estimation_error(
    state = state,
    error_model = "bespoke_increasing_sd",
    error_param = error_param_default
  )

  set.seed(999)
  result2 <- apply_estimation_error(
    state = state,
    error_model = "bespoke_increasing_sd",
    error_param = error_param_default
  )

  # Results should differ

  expect_false(all(result1$estimated_age == result2$estimated_age))
})

