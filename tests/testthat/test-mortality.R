# Unit tests for mortality module

# -----------------------------------------------------------------------------
# Test fixtures
# -----------------------------------------------------------------------------

# Siler parameters (Gage & Dyke 1986, Table 2, Level 15)
a_siler_test <- c(0.175, 1.40, 0.00368, 0.000075, 0.0917)
b_siler_test <- demohaz::trad_to_demohaz_siler_param(a_siler_test)

# mortality_param = c(k2, b_siler) where k2 = 1.2
mortality_param_test <- c(1.2, b_siler_test)

create_test_state <- function(n = 100) {
  data.frame(
    agent_id = 1:n,
    age = as.numeric(sample(0:80, n, replace = TRUE)),
    lesion = sample(c(TRUE, FALSE), n, replace = TRUE),
    dead = rep(FALSE, n),
    in_sample = rep(TRUE, n)
  )
}

# -----------------------------------------------------------------------------
# Input validation tests
# -----------------------------------------------------------------------------

test_that("mortality module raises error when mortality_model is not 'usher3'", {
  state <- create_test_state()
  mortality_param <- mortality_param_test
  dx <- 5

  expect_error(
    apply_mortality_usher3(
      state = state,
      mortality_model = "gompertz",
      mortality_param = mortality_param,
      dx = dx
    ),
    regexp = "usher3"
  )

  expect_error(
    apply_mortality_usher3(
      state = state,
      mortality_model = "siler",
      mortality_param = mortality_param,
      dx = dx
    ),
    regexp = "usher3"
  )

  expect_error(
    apply_mortality_usher3(
      state = state,
      mortality_model = NULL,
      mortality_param = mortality_param,
      dx = dx
    )
  )
})

# -----------------------------------------------------------------------------
# Successful execution tests
# -----------------------------------------------------------------------------

test_that("mortality module executes successfully with usher3 model", {
  state <- create_test_state()
  mortality_param <- mortality_param_test
  dx <- 5

  expect_no_error(
    apply_mortality_usher3(
      state = state,
      mortality_model = "usher3",
      mortality_param = mortality_param,
      dx = dx
    )
  )
})

test_that("mortality module returns a dataframe", {
  state <- create_test_state()
  mortality_param <- mortality_param_test
  dx <- 5

  result <- apply_mortality_usher3(
    state = state,
    mortality_model = "usher3",
    mortality_param = mortality_param,
    dx = dx
  )

  expect_s3_class(result, "data.frame")
})

# -----------------------------------------------------------------------------
# Output structure tests
# -----------------------------------------------------------------------------

test_that("output state contains same columns as input state", {
  state <- create_test_state()
  mortality_param <- mortality_param_test
  dx <- 5

  result <- apply_mortality_usher3(
    state = state,
    mortality_model = "usher3",
    mortality_param = mortality_param,
    dx = dx
  )

  expect_equal(sort(names(result)), sort(names(state)))
})

test_that("output state has same dimensions as input state", {
  state <- create_test_state(n = 150)
  mortality_param <- mortality_param_test
  dx <- 5

  result <- apply_mortality_usher3(
    state = state,
    mortality_model = "usher3",
    mortality_param = mortality_param,
    dx = dx
  )

  expect_equal(nrow(result), nrow(state))
  expect_equal(ncol(result), ncol(state))
  expect_equal(dim(result), dim(state))
})
  
test_that("output state has expected column types", {
  state <- create_test_state()
  mortality_param <- mortality_param_test
  dx <- 5

  result <- apply_mortality_usher3(
    state = state,
    mortality_model = "usher3",
    mortality_param = mortality_param,
    dx = dx
  )

  expect_type(result$agent_id, "integer")
  expect_type(result$age, "double")
  expect_type(result$lesion, "logical")
  expect_type(result$dead, "logical")
  expect_type(result$in_sample, "logical")
})

# -----------------------------------------------------------------------------
# Age increment (dx) tests
# -----------------------------------------------------------------------------

test_that("surviving agents have age incremented by dx", {
  state <- create_test_state(n = 100)
  mortality_param <- mortality_param_test
  dx <- 5

  result <- apply_mortality_usher3(
    state = state,
    mortality_model = "usher3",
    mortality_param = mortality_param,
    dx = dx
  )

  # Survivors (dead == FALSE) should have age increased by dx
  survivors <- result$dead == FALSE
  if (any(survivors)) {
    expect_equal(result$age[survivors], state$age[survivors] + dx)
  }
})

test_that("dead agents retain their original age", {
  state <- create_test_state(n = 100)
  mortality_param <- mortality_param_test
  dx <- 5

  result <- apply_mortality_usher3(
    state = state,
    mortality_model = "usher3",
    mortality_param = mortality_param,
    dx = dx
  )

  # Dead agents should have their original age (not incremented)
  dead <- result$dead == TRUE
  if (any(dead)) {
    expect_equal(result$age[dead], state$age[dead])
  }
})

test_that("age increment works correctly with dx = 1", {
  state <- create_test_state(n = 100)
  mortality_param <- mortality_param_test
  dx <- 1

  result <- apply_mortality_usher3(
    state = state,
    mortality_model = "usher3",
    mortality_param = mortality_param,
    dx = dx
  )

  survivors <- result$dead == FALSE
  dead <- result$dead == TRUE
  
  if (any(survivors)) {
    expect_equal(result$age[survivors], state$age[survivors] + 1)
  }
  if (any(dead)) {
    expect_equal(result$age[dead], state$age[dead])
  }
})

test_that("age increment works correctly with dx = 10", {
  state <- create_test_state(n = 100)
  mortality_param <- mortality_param_test
  dx <- 10

  result <- apply_mortality_usher3(
    state = state,
    mortality_model = "usher3",
    mortality_param = mortality_param,
    dx = dx
  )

  survivors <- result$dead == FALSE
  dead <- result$dead == TRUE
  
  if (any(survivors)) {
    expect_equal(result$age[survivors], state$age[survivors] + 10)
  }
  if (any(dead)) {
    expect_equal(result$age[dead], state$age[dead])
  }
})

test_that("age increment works correctly with fractional dx", {
  state <- create_test_state(n = 100)
  mortality_param <- mortality_param_test
  dx <- 0.5

  result <- apply_mortality_usher3(
    state = state,
    mortality_model = "usher3",
    mortality_param = mortality_param,
    dx = dx
  )

  survivors <- result$dead == FALSE
  dead <- result$dead == TRUE
  
  if (any(survivors)) {
    expect_equal(result$age[survivors], state$age[survivors] + 0.5)
  }
  if (any(dead)) {
    expect_equal(result$age[dead], state$age[dead])
  }
})

test_that("age increment works correctly with small dx values", {
  state <- create_test_state(n = 100)
  mortality_param <- mortality_param_test
  dx <- 0.1

  result <- apply_mortality_usher3(
    state = state,
    mortality_model = "usher3",
    mortality_param = mortality_param,
    dx = dx
  )

  survivors <- result$dead == FALSE
  dead <- result$dead == TRUE
  
  if (any(survivors)) {
    expect_equal(result$age[survivors], state$age[survivors] + 0.1)
  }
  if (any(dead)) {
    expect_equal(result$age[dead], state$age[dead])
  }
})

test_that("module works with multiple dx values", {
  state <- create_test_state(n = 100)
  mortality_param <- mortality_param_test
  
  dx_values <- c(1, 2, 5, 0.5, 10)
  
  for (dx in dx_values) {
    result <- apply_mortality_usher3(
      state = state,
      mortality_model = "usher3",
      mortality_param = mortality_param,
      dx = dx
    )
    
    survivors <- result$dead == FALSE
    dead <- result$dead == TRUE
    
    if (any(survivors)) {
      expect_equal(
        result$age[survivors], 
        state$age[survivors] + dx,
        info = paste("Survivors failed for dx =", dx)
      )
    }
    
    if (any(dead)) {
      expect_equal(
        result$age[dead], 
        state$age[dead],
        info = paste("Dead agents failed for dx =", dx)
      )
    }
  }
})
