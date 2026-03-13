# Unit tests for lesion module

# -----------------------------------------------------------------------------
# Test fixtures
# -----------------------------------------------------------------------------

# Lesion parameters for each model
lesion_param_constant <- c(0.05)  # k1 only
lesion_param_constant_to <- c(0.05, 6)  # k1, age_end
lesion_param_constant_from <- c(0.05, 18)  # k1, age_start
lesion_param_constant_interval <- c(0.05, 2, 6)  # k1, age_start, age_end

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

test_that("lesion module raises error when lesion_model is NULL", {
  state <- create_test_state()
  dx <- 1

  expect_error(
    apply_lesions(
      state = state,
      lesion_model = NULL,
      lesion_param = lesion_param_constant,
      dx = dx
    )
  )
})

test_that("lesion module raises error when lesion_model is unsupported", {
  state <- create_test_state()
  dx <- 1

  expect_error(
    apply_lesions(
      state = state,
      lesion_model = "gompertz",
      lesion_param = lesion_param_constant,
      dx = dx
    )
  )

  expect_error(
    apply_lesions(
      state = state,
      lesion_model = "unknown",
      lesion_param = lesion_param_constant,
      dx = dx
    )
  )
})

test_that("constant model requires length-1 parameter vector", {
  state <- create_test_state()
  dx <- 1

  expect_error(
    apply_lesions(
      state = state,
      lesion_model = "constant",
      lesion_param = c(0.05, 6),  # wrong length
      dx = dx
    )
  )

  expect_error(
    apply_lesions(
      state = state,
      lesion_model = "constant",
      lesion_param = c(),  # empty
      dx = dx
    )
  )
})

test_that("constant_to model requires length-2 parameter vector", {
  state <- create_test_state()
  dx <- 1

  expect_error(
    apply_lesions(
      state = state,
      lesion_model = "constant_to",
      lesion_param = c(0.05),  # too short
      dx = dx
    )
  )

  expect_error(
    apply_lesions(
      state = state,
      lesion_model = "constant_to",
      lesion_param = c(0.05, 2, 6),  # too long
      dx = dx
    )
  )
})

test_that("constant_from model requires length-2 parameter vector", {
  state <- create_test_state()
  dx <- 1

  expect_error(
    apply_lesions(
      state = state,
      lesion_model = "constant_from",
      lesion_param = c(0.05),  # too short
      dx = dx
    )
  )

  expect_error(
    apply_lesions(
      state = state,
      lesion_model = "constant_from",
      lesion_param = c(0.05, 2, 6),  # too long
      dx = dx
    )
  )
})

test_that("constant_interval model requires length-3 parameter vector", {
  state <- create_test_state()
  dx <- 1

  expect_error(
    apply_lesions(
      state = state,
      lesion_model = "constant_interval",
      lesion_param = c(0.05, 6),  # too short
      dx = dx
    )
  )

  expect_error(
    apply_lesions(
      state = state,
      lesion_model = "constant_interval",
      lesion_param = c(0.05),  # too short
      dx = dx
    )
  )
})

# -----------------------------------------------------------------------------
# Successful execution tests
# -----------------------------------------------------------------------------

test_that("constant model executes successfully", {
  state <- create_test_state()
  dx <- 1

  expect_no_error(
    apply_lesions(
      state = state,
      lesion_model = "constant",
      lesion_param = lesion_param_constant,
      dx = dx
    )
  )
})

test_that("constant_to model executes successfully", {
  state <- create_test_state()
  dx <- 1

  expect_no_error(
    apply_lesions(
      state = state,
      lesion_model = "constant_to",
      lesion_param = lesion_param_constant_to,
      dx = dx
    )
  )
})

test_that("constant_from model executes successfully", {
  state <- create_test_state()
  dx <- 1

  expect_no_error(
    apply_lesions(
      state = state,
      lesion_model = "constant_from",
      lesion_param = lesion_param_constant_from,
      dx = dx
    )
  )
})

test_that("constant_interval model executes successfully", {
  state <- create_test_state()
  dx <- 1

  expect_no_error(
    apply_lesions(
      state = state,
      lesion_model = "constant_interval",
      lesion_param = lesion_param_constant_interval,
      dx = dx
    )
  )
})

test_that("lesion module returns a dataframe", {
  state <- create_test_state()
  dx <- 1

  result <- apply_lesions(
    state = state,
    lesion_model = "constant",
    lesion_param = lesion_param_constant,
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

  result <- apply_lesions(
    state = state,
    lesion_model = "constant",
    lesion_param = lesion_param_constant,
    dx = dx
  )

  expect_equal(sort(names(result)), sort(names(state)))
})

test_that("output state has same dimensions as input state", {
  state <- create_test_state(n = 150)
  dx <- 1

  result <- apply_lesions(
    state = state,
    lesion_model = "constant",
    lesion_param = lesion_param_constant,
    dx = dx
  )

  expect_equal(nrow(result), nrow(state))
  expect_equal(ncol(result), ncol(state))
  expect_equal(dim(result), dim(state))
})

test_that("output state has expected column types", {
  state <- create_test_state()
  dx <- 1

  result <- apply_lesions(
    state = state,
    lesion_model = "constant",
    lesion_param = lesion_param_constant,
    dx = dx
  )

  expect_type(result$agent_id, "integer")
  expect_type(result$age, "double")
  expect_type(result$lesion, "logical")
  expect_type(result$dead, "logical")
  expect_type(result$in_sample, "logical")
})

# -----------------------------------------------------------------------------
# Lesion acquisition logic tests
# -----------------------------------------------------------------------------

test_that("agents with existing lesions retain lesion status", {
  state <- create_test_state(n = 100)
  dx <- 1

  # Track which agents already have lesions
  had_lesion <- state$lesion

  result <- apply_lesions(
    state = state,
    lesion_model = "constant",
    lesion_param = lesion_param_constant,
    dx = dx
  )

  # All agents who had lesions should still have them
  expect_true(all(result$lesion[had_lesion] == TRUE))
})

test_that("dead agents are not processed for lesion acquisition", {
  state <- create_test_state(n = 100)
  # Mark some agents as dead and without lesions
  state$dead[1:20] <- TRUE
  state$lesion[1:20] <- FALSE
  dx <- 1

  result <- apply_lesions(
    state = state,
    lesion_model = "constant",
    lesion_param = lesion_param_constant,
    dx = dx
  )

  # Dead agents should not have acquired lesions
  expect_true(all(result$lesion[1:20] == FALSE))
})

test_that("lesion acquisition respects age window for constant_to model", {
  # Create state with controlled ages
  state <- data.frame(
    agent_id = 1:100,
    age = c(rep(3, 50), rep(10, 50)),  # 50 agents at age 3, 50 at age 10
    lesion = rep(FALSE, 100),
    dead = rep(FALSE, 100),
    in_sample = rep(TRUE, 100)
  )
  dx <- 1
  # Window ends at age 6
  lesion_param <- c(0.99, 6)  # High k1 to ensure some transitions

  result <- apply_lesions(
    state = state,
    lesion_model = "constant_to",
    lesion_param = lesion_param,
    dx = dx
  )

  # Agents at age 10 (outside window) should not have acquired lesions
  expect_true(all(result$lesion[51:100] == FALSE))
})

test_that("lesion acquisition respects age window for constant_from model", {
  # Create state with controlled ages
  state <- data.frame(
    agent_id = 1:100,
    age = c(rep(10, 50), rep(25, 50)),  # 50 agents at age 10, 50 at age 25
    lesion = rep(FALSE, 100),
    dead = rep(FALSE, 100),
    in_sample = rep(TRUE, 100)
  )
  dx <- 1
  # Window starts at age 18
  lesion_param <- c(0.99, 18)  # High k1 to ensure some transitions

  result <- apply_lesions(
    state = state,
    lesion_model = "constant_from",
    lesion_param = lesion_param,
    dx = dx
  )

  # Agents at age 10 (outside window) should not have acquired lesions
  expect_true(all(result$lesion[1:50] == FALSE))
})

test_that("lesion acquisition respects age window for constant_interval model", {
  # Create state with controlled ages
  state <- data.frame(
    agent_id = 1:150,
    age = c(rep(1, 50), rep(4, 50), rep(10, 50)),  # ages 1, 4, and 10
    lesion = rep(FALSE, 150),
    dead = rep(FALSE, 150),
    in_sample = rep(TRUE, 150)
  )
  dx <- 1
  # Window is [2, 6)
  lesion_param <- c(0.99, 2, 6)  # High k1 to ensure some transitions

  result <- apply_lesions(
    state = state,
    lesion_model = "constant_interval",
    lesion_param = lesion_param,
    dx = dx
  )

  # Agents at age 1 (before window) should not have acquired lesions
  expect_true(all(result$lesion[1:50] == FALSE))
  # Agents at age 10 (after window) should not have acquired lesions
  expect_true(all(result$lesion[101:150] == FALSE))
})

# -----------------------------------------------------------------------------
# Age unchanged tests
# -----------------------------------------------------------------------------

test_that("agent ages are not modified by lesion module", {
  state <- create_test_state(n = 100)
  dx <- 5

  result <- apply_lesions(
    state = state,
    lesion_model = "constant",
    lesion_param = lesion_param_constant,
    dx = dx
  )

  # Ages should be unchanged - lesion module does not modify age
  expect_equal(result$age, state$age)
})

# -----------------------------------------------------------------------------
# Transition probability tests
# -----------------------------------------------------------------------------

test_that("no transitions occur when k1 is zero", {
  state <- data.frame(
    agent_id = 1:100,
    age = as.numeric(rep(5, 100)),
    lesion = rep(FALSE, 100),
    dead = rep(FALSE, 100),
    in_sample = rep(TRUE, 100)
  )
  dx <- 1
  lesion_param <- c(0)  # k1 = 0

  result <- apply_lesions(
    state = state,
    lesion_model = "constant",
    lesion_param = lesion_param,
    dx = dx
  )

  # No one should acquire a lesion
  expect_true(all(result$lesion == FALSE))
})

test_that("no transitions when k1 is zero for constant_to model", {
  state <- data.frame(
    agent_id = 1:100,
    age = as.numeric(rep(3, 100)),  # within window [0, 6)
    lesion = rep(FALSE, 100),
    dead = rep(FALSE, 100),
    in_sample = rep(TRUE, 100)
  )
  dx <- 1
  lesion_param <- c(0, 6)  # k1 = 0, window ends at 6

  result <- apply_lesions(
    state = state,
    lesion_model = "constant_to",
    lesion_param = lesion_param,
    dx = dx
  )

  expect_true(all(result$lesion == FALSE))
})

# -----------------------------------------------------------------------------
# Window boundary tests
# -----------------------------------------------------------------------------

test_that("constant_to: agents exactly at age_end do not acquire lesions", {
  state <- data.frame(
    agent_id = 1:100,
    age = as.numeric(rep(6, 100)),  # exactly at boundary
    lesion = rep(FALSE, 100),
    dead = rep(FALSE, 100),
    in_sample = rep(TRUE, 100)
  )
  dx <- 1
  lesion_param <- c(0.99, 6)  # high k1, window ends at 6

  result <- apply_lesions(
    state = state,
    lesion_model = "constant_to",
    lesion_param = lesion_param,
    dx = dx
  )

  # Agents at age 6 are at or past the boundary, no transitions
  expect_true(all(result$lesion == FALSE))
})

test_that("constant_from: agents exactly at age_start can acquire lesions", {
  set.seed(12345)
  state <- data.frame(
    agent_id = 1:1000,
    age = as.numeric(rep(18, 1000)),  # exactly at boundary
    lesion = rep(FALSE, 1000),
    dead = rep(FALSE, 1000),
    in_sample = rep(TRUE, 1000)
  )
  dx <- 1
  lesion_param <- c(0.99, 18)  # high k1, window starts at 18

  result <- apply_lesions(
    state = state,
    lesion_model = "constant_from",
    lesion_param = lesion_param,
    dx = dx
  )

  # Agents at age 18 are eligible, most should transition
  expect_gt(mean(result$lesion), 0.9)
})

test_that("constant_from: agents just below age_start cannot acquire lesions", {
  state <- data.frame(
    agent_id = 1:100,
    age = as.numeric(rep(17.9, 100)),  # just below boundary
    lesion = rep(FALSE, 100),
    dead = rep(FALSE, 100),
    in_sample = rep(TRUE, 100)
  )
  dx <- 1
  lesion_param <- c(0.99, 18)  # high k1, window starts at 18

  result <- apply_lesions(
    state = state,
    lesion_model = "constant_from",
    lesion_param = lesion_param,
    dx = dx
  )

  # Agents below 18 should not transition
  expect_true(all(result$lesion == FALSE))
})

test_that("constant_interval: agents at boundaries behave correctly", {
  # Test lower boundary (included)
  set.seed(12345)
  state_lower <- data.frame(
    agent_id = 1:1000,
    age = as.numeric(rep(2, 1000)),  # exactly at age_start
    lesion = rep(FALSE, 1000),
    dead = rep(FALSE, 1000),
    in_sample = rep(TRUE, 1000)
  )
  lesion_param <- c(0.99, 2, 6)  # window [2, 6)

  result_lower <- apply_lesions(
    state = state_lower,
    lesion_model = "constant_interval",
    lesion_param = lesion_param,
    dx = 1
  )

  # At age 2 (lower boundary), transitions should occur
  expect_gt(mean(result_lower$lesion), 0.9)

  # Test upper boundary (excluded)
  state_upper <- data.frame(
    agent_id = 1:100,
    age = as.numeric(rep(6, 100)),  # exactly at age_end
    lesion = rep(FALSE, 100),
    dead = rep(FALSE, 100),
    in_sample = rep(TRUE, 100)
  )

  result_upper <- apply_lesions(
    state = state_upper,
    lesion_model = "constant_interval",
    lesion_param = lesion_param,
    dx = 1
  )

  # At age 6 (upper boundary), no transitions
  expect_true(all(result_upper$lesion == FALSE))
})

# -----------------------------------------------------------------------------
# Combined eligibility tests
# -----------------------------------------------------------------------------

test_that("only alive agents without lesions in window can acquire lesions", {
  set.seed(12345)
  state <- data.frame(
    agent_id = 1:8,
    age = as.numeric(c(3, 3, 3, 3, 10, 10, 3, 3)),
    lesion = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
    dead = c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE),
    in_sample = rep(TRUE, 8)
  )
  # Agent 1: alive, no lesion, in window -> CAN transition
  # Agent 2: dead, no lesion, in window -> CANNOT (dead)
  # Agent 3: alive, has lesion, in window -> CANNOT (already has)
  # Agent 4: dead, has lesion, in window -> CANNOT (dead + already has)
  # Agent 5: alive, no lesion, outside window -> CANNOT (outside)
  # Agent 6: dead, no lesion, outside window -> CANNOT (dead + outside)
  # Agent 7-8: alive, no lesion, in window -> CAN transition

  lesion_param <- c(0.99, 6)  # high k1, window ends at 6

  result <- apply_lesions(
    state = state,
    lesion_model = "constant_to",
    lesion_param = lesion_param,
    dx = 1
  )

  # Agent 2: should remain FALSE (dead)
  expect_false(result$lesion[2])
  # Agent 3: should remain TRUE (already had)
  expect_true(result$lesion[3])
  # Agent 4: should remain TRUE (already had)
  expect_true(result$lesion[4])
  # Agent 5: should remain FALSE (outside window)
  expect_false(result$lesion[5])
  # Agent 6: should remain FALSE (dead + outside)
  expect_false(result$lesion[6])
})

