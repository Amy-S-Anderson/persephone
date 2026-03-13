# Unit tests for preservation module

# -----------------------------------------------------------------------------
# Test fixtures
# -----------------------------------------------------------------------------

# Siler parameters (Gage & Dyke 1986, Table 2, Level 15)
a_siler_test <- c(0.175, 1.40, 0.00368, 0.000075, 0.0917)
b_siler_test <- demohaz::trad_to_demohaz_siler_param(a_siler_test)

# preservation_param = b_siler (length 5)
preservation_param_test <- b_siler_test

create_test_state <- function(n = 100) {
  data.frame(
    agent_id = 1:n,
    age = as.numeric(sample(0:80, n, replace = TRUE)),
    lesion = sample(c(TRUE, FALSE), n, replace = TRUE),
    dead = rep(TRUE, n),
    was_deposited = rep(TRUE, n),  # All deposited for preservation tests
    in_sample = rep(TRUE, n)
  )
}

# -----------------------------------------------------------------------------
# Input validation tests
# -----------------------------------------------------------------------------

test_that("preservation module raises error when preservation_model is NULL", {
  state <- create_test_state()
  dx <- 1

  expect_error(
    apply_preservation(
      state = state,
      preservation_model = NULL,
      preservation_param = preservation_param_test,
      dx = dx
    )
  )
})

test_that("preservation module raises error when preservation_model is unsupported", {
  state <- create_test_state()
  dx <- 1

  expect_error(
    apply_preservation(
      state = state,
      preservation_model = "unknown",
      preservation_param = preservation_param_test,
      dx = dx
    )
  )

  expect_error(
    apply_preservation(
      state = state,
      preservation_model = "exponential",
      preservation_param = preservation_param_test,
      dx = dx
    )
  )
})

test_that("siler model requires length-5 parameter vector", {
  state <- create_test_state()
  dx <- 1

  expect_error(
    apply_preservation(
      state = state,
      preservation_model = "siler",
      preservation_param = c(0.1, 0.2, 0.3),  # wrong length
      dx = dx
    )
  )

  expect_error(
    apply_preservation(
      state = state,
      preservation_model = "siler",
      preservation_param = c(),  # empty
      dx = dx
    )
  )
})

# -----------------------------------------------------------------------------
# Successful execution tests
# -----------------------------------------------------------------------------

test_that("siler model executes successfully", {
  state <- create_test_state()
  dx <- 1

  expect_no_error(
    apply_preservation(
      state = state,
      preservation_model = "siler",
      preservation_param = preservation_param_test,
      dx = dx
    )
  )
})

test_that("preservation module returns a dataframe", {
  state <- create_test_state()
  dx <- 1

  result <- apply_preservation(
    state = state,
    preservation_model = "siler",
    preservation_param = preservation_param_test,
    dx = dx
  )

  expect_s3_class(result, "data.frame")
})

# -----------------------------------------------------------------------------
# Output structure tests
# -----------------------------------------------------------------------------

test_that("output state contains same columns as input state", {
  state <- create_test_state()
  dx <- 1

  result <- apply_preservation(
    state = state,
    preservation_model = "siler",
    preservation_param = preservation_param_test,
    dx = dx
  )

  expect_equal(sort(names(result)), sort(names(state)))
})

test_that("output state has same dimensions as input state", {
  state <- create_test_state(n = 150)
  dx <- 1

  result <- apply_preservation(
    state = state,
    preservation_model = "siler",
    preservation_param = preservation_param_test,
    dx = dx
  )

  expect_equal(nrow(result), nrow(state))
  expect_equal(ncol(result), ncol(state))
  expect_equal(dim(result), dim(state))
})

test_that("output state has expected column types", {
  state <- create_test_state()
  dx <- 1

  result <- apply_preservation(
    state = state,
    preservation_model = "siler",
    preservation_param = preservation_param_test,
    dx = dx
  )

  expect_type(result$agent_id, "integer")
  expect_type(result$age, "double")
  expect_type(result$lesion, "logical")
  expect_type(result$dead, "logical")
  expect_type(result$was_deposited, "logical")
  expect_type(result$in_sample, "logical")
})

# -----------------------------------------------------------------------------
# Preservation logic tests
# -----------------------------------------------------------------------------

test_that("non-deposited agents have in_sample set to FALSE", {
  state <- data.frame(
    agent_id = 1:10,
    age = as.numeric(rep(30, 10)),
    lesion = rep(FALSE, 10),
    dead = rep(TRUE, 10),
    was_deposited = c(rep(TRUE, 5), rep(FALSE, 5)),
    in_sample = rep(TRUE, 10)  # Start all TRUE
  )
  dx <- 1

  result <- apply_preservation(
    state = state,
    preservation_model = "siler",
    preservation_param = preservation_param_test,
    dx = dx
  )

  # Non-deposited agents must have in_sample = FALSE
  expect_true(all(result$in_sample[6:10] == FALSE))
})

test_that("agents already in_sample FALSE remain FALSE", {
  state <- data.frame(
    agent_id = 1:10,
    age = as.numeric(rep(30, 10)),
    lesion = rep(FALSE, 10),
    dead = rep(TRUE, 10),
    was_deposited = rep(TRUE, 10),
    in_sample = rep(FALSE, 10)  # All start FALSE
  )
  dx <- 1

  result <- apply_preservation(
    state = state,
    preservation_model = "siler",
    preservation_param = preservation_param_test,
    dx = dx
  )

  # All should remain FALSE

  expect_true(all(result$in_sample == FALSE))
})

test_that("zero hazard preserves all deposited agents", {
  state <- data.frame(
    agent_id = 1:100,
    age = as.numeric(rep(30, 100)),
    lesion = rep(FALSE, 100),
    dead = rep(TRUE, 100),
    was_deposited = rep(TRUE, 100),
    in_sample = rep(TRUE, 100)
  )
  dx <- 1
  # Zero Siler hazard: all parameters zero
  zero_param <- c(0, 1, 0, 0, 0)

  result <- apply_preservation(
    state = state,
    preservation_model = "siler",
    preservation_param = zero_param,
    dx = dx
  )

  # All should be preserved (in_sample = TRUE)
  expect_true(all(result$in_sample == TRUE))
})

test_that("deposited agents are subject to Siler filtering", {
  set.seed(12345)
  state <- data.frame(
    agent_id = 1:1000,
    age = as.numeric(rep(30, 1000)),
    lesion = rep(FALSE, 1000),
    dead = rep(TRUE, 1000),
    was_deposited = rep(TRUE, 1000),
    in_sample = rep(TRUE, 1000)
  )
  dx <- 1

  result <- apply_preservation(
    state = state,
    preservation_model = "siler",
    preservation_param = preservation_param_test,
    dx = dx
  )

  # Some should be filtered out (stochastic, but with 1000 agents we expect some FALSE)
  # The Siler hazard at age 30 is non-zero, so some filtering should occur
  expect_true(any(result$in_sample == FALSE))
  expect_true(any(result$in_sample == TRUE))
})

# -----------------------------------------------------------------------------
# State preservation tests
# -----------------------------------------------------------------------------

test_that("agent ages are not modified by preservation module", {
  state <- create_test_state(n = 100)
  dx <- 1

  result <- apply_preservation(
    state = state,
    preservation_model = "siler",
    preservation_param = preservation_param_test,
    dx = dx
  )

  expect_equal(result$age, state$age)
})

test_that("lesion status is not modified by preservation module", {
  state <- create_test_state(n = 100)
  dx <- 1

  result <- apply_preservation(
    state = state,
    preservation_model = "siler",
    preservation_param = preservation_param_test,
    dx = dx
  )

  expect_equal(result$lesion, state$lesion)
})

test_that("dead status is not modified by preservation module", {
  state <- create_test_state(n = 100)
  dx <- 1

  result <- apply_preservation(
    state = state,
    preservation_model = "siler",
    preservation_param = preservation_param_test,
    dx = dx
  )

  expect_equal(result$dead, state$dead)
})

test_that("was_deposited is not modified by preservation module", {
  state <- create_test_state(n = 100)
  # Mix of deposited and not
  state$was_deposited <- sample(c(TRUE, FALSE), 100, replace = TRUE)
  dx <- 1

  result <- apply_preservation(
    state = state,
    preservation_model = "siler",
    preservation_param = preservation_param_test,
    dx = dx
  )

  expect_equal(result$was_deposited, state$was_deposited)
})

