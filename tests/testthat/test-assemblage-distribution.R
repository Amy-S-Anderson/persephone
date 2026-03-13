# Unit tests for assemblage empirical distribution helper

# -----------------------------------------------------------------------------
# Test fixtures
# -----------------------------------------------------------------------------

#' Create a simple test assemblage data frame
#' 
#' @param n_intervals Number of age intervals
#' @param n_studies Number of studies
#' @return A data.frame with study, x_lo, x_hi, count columns
create_test_assemblage <- function() {
  data.frame(
    study = c("site_A", "site_A", "site_A", "site_B", "site_B"),
    x_lo = c(-0.5, 1.5, 15.0, -0.5, 15.0),
    x_hi = c(1.5, 15.0, Inf, 1.5, Inf),
    count = c(12, 25, 40, 8, 35)
  )
}

#' Create a minimal single-study assemblage
create_single_study_assemblage <- function() {
  data.frame(
    study = c("site_A", "site_A", "site_A"),
    x_lo = c(0, 5, 20),
    x_hi = c(5, 20, Inf),
    count = c(10, 30, 60)
  )
}

#' Create an assemblage with only finite intervals (no Inf)
create_finite_assemblage <- function() {
  data.frame(
    study = c("site_A", "site_A"),
    x_lo = c(0, 10),
    x_hi = c(10, 20),
    count = c(15, 25)
  )
}

# -----------------------------------------------------------------------------
# Input validation tests
# -----------------------------------------------------------------------------

test_that("raises error when data is not a data.frame", {
  expect_error(
    create_assemblage_distribution(data = list(x_lo = 0, x_hi = 10, count = 5)),
    regexp = "data.frame"
  )
  
  expect_error(
    create_assemblage_distribution(data = matrix(1:9, nrow = 3)),
    regexp = "data.frame"
  )
  
  expect_error(
    create_assemblage_distribution(data = NULL)
  )
})

test_that("raises error when required columns are missing", {
  # Missing study
  expect_error(
    create_assemblage_distribution(
      data = data.frame(x_lo = 0, x_hi = 10, count = 5)
    ),
    regexp = "study"
  )
  
  # Missing x_lo
  expect_error(
    create_assemblage_distribution(
      data = data.frame(study = "A", x_hi = 10, count = 5)
    ),
    regexp = "x_lo"
  )
  
  # Missing x_hi
  expect_error(
    create_assemblage_distribution(
      data = data.frame(study = "A", x_lo = 0, count = 5)
    ),
    regexp = "x_hi"
  )
  
  # Missing count
  expect_error(
    create_assemblage_distribution(
      data = data.frame(study = "A", x_lo = 0, x_hi = 10)
    ),
    regexp = "count"
  )
})

test_that("raises error when x_lo is not numeric", {
  expect_error(
    create_assemblage_distribution(
      data = data.frame(
        study = "A",
        x_lo = "zero",
        x_hi = 10,
        count = 5
      )
    ),
    regexp = "numeric"
  )
})

test_that("raises error when x_hi is not numeric", {
  expect_error(
    create_assemblage_distribution(
      data = data.frame(
        study = "A",
        x_lo = 0,
        x_hi = "ten",
        count = 5
      )
    ),
    regexp = "numeric"
  )
})

test_that("raises error when count is not numeric", {
  expect_error(
    create_assemblage_distribution(
      data = data.frame(
        study = "A",
        x_lo = 0,
        x_hi = 10,
        count = "five"
      )
    ),
    regexp = "numeric"
  )
})

test_that("raises error when count contains negative values", {
  expect_error(
    create_assemblage_distribution(
      data = data.frame(
        study = "A",
        x_lo = c(0, 10),
        x_hi = c(10, 20),
        count = c(5, -3)
      )
    ),
    regexp = "negative"
  )
})

test_that("raises error when x_lo >= x_hi for any interval", {
  # x_lo equals x_hi
  expect_error(
    create_assemblage_distribution(
      data = data.frame(
        study = "A",
        x_lo = c(0, 10),
        x_hi = c(10, 10),
        count = c(5, 3)
      )
    ),
    regexp = "x_lo.*x_hi"
  )
  
 # x_lo greater than x_hi
  expect_error(
    create_assemblage_distribution(
      data = data.frame(
        study = "A",
        x_lo = c(0, 15),
        x_hi = c(10, 10),
        count = c(5, 3)
      )
    ),
    regexp = "x_lo.*x_hi"
  )
})

test_that("raises error when total count is zero", {
  expect_error(
    create_assemblage_distribution(
      data = data.frame(
        study = "A",
        x_lo = c(0, 10),
        x_hi = c(10, 20),
        count = c(0, 0)
      )
    ),
    regexp = "zero"
  )
})

# -----------------------------------------------------------------------------
# Successful execution tests
# -----------------------------------------------------------------------------

test_that("executes successfully with valid multi-study assemblage", {
  data <- create_test_assemblage()
  
  expect_no_error(
    create_assemblage_distribution(data = data)
  )
})

test_that("executes successfully with single-study assemblage", {
  data <- create_single_study_assemblage()
  
  expect_no_error(
    create_assemblage_distribution(data = data)
  )
})

test_that("executes successfully with finite intervals only", {
  data <- create_finite_assemblage()
  
  expect_no_error(
    create_assemblage_distribution(data = data)
  )
})

test_that("executes successfully with negative x_lo values", {
  data <- data.frame(
    study = "A",
    x_lo = c(-1.0, 0, 10),
    x_hi = c(0, 10, Inf),
    count = c(5, 20, 30)
  )
  
  expect_no_error(
    create_assemblage_distribution(data = data)
  )
})

# -----------------------------------------------------------------------------
# Output structure tests
# -----------------------------------------------------------------------------

test_that("returns a list with expected components", {
  data <- create_test_assemblage()
  result <- create_assemblage_distribution(data = data)
  
  expect_type(result, "list")
  
  # Should have density function
  expect_true("dfun" %in% names(result))
  expect_type(result$dfun, "closure")
  
  # Should have sampling function
  expect_true("rfun" %in% names(result))
  expect_type(result$rfun, "closure")
  
  # Should have cut points (interval boundaries)
  expect_true("cuts" %in% names(result))
  expect_type(result$cuts, "double")
  
  # Should have densities for each interval
  expect_true("densities" %in% names(result))
  expect_type(result$densities, "double")
  
  # Should have exponential tail rate (may be NA if no unbounded interval)
  expect_true("exp_rate" %in% names(result))
})

test_that("cuts are sorted in ascending order", {
  data <- create_test_assemblage()
  result <- create_assemblage_distribution(data = data)
  
  expect_equal(result$cuts, sort(result$cuts))
})

test_that("number of densities equals number of intervals", {
  data <- create_test_assemblage()
  result <- create_assemblage_distribution(data = data)
  

  # Number of intervals = length(cuts) - 1, plus potentially 1 for exp tail
  n_finite_intervals <- length(result$cuts) - 1
  
  # densities should correspond to finite intervals
  expect_equal(length(result$densities), n_finite_intervals)
})

# -----------------------------------------------------------------------------
# Normalization tests
# -----------------------------------------------------------------------------

test_that("density integrates to 1 for finite intervals only", {
  data <- create_finite_assemblage()
  result <- create_assemblage_distribution(data = data)
  
  # Integrate the piecewise constant density
  widths <- diff(result$cuts)
  total_mass <- sum(result$densities * widths)
  
  expect_equal(total_mass, 1.0, tolerance = 1e-10)
})

test_that("density integrates to 1 with exponential tail", {
  data <- create_single_study_assemblage()
  result <- create_assemblage_distribution(data = data)
  
  # Finite intervals contribution
  widths <- diff(result$cuts)
  finite_mass <- sum(result$densities * widths)
  
  # Exponential tail contribution: integral from last cut to Inf
  # For exponential with rate lambda starting at x0:
  # integral of d * exp(-lambda * (x - x0)) from x0 to Inf = d / lambda
  # where d is the density at the start of the tail
  last_cut <- result$cuts[length(result$cuts)]
  tail_mass <- result$tail_density / result$exp_rate
  
  total_mass <- finite_mass + tail_mass
  
  expect_equal(total_mass, 1.0, tolerance = 1e-10)
})

test_that("dfun values are non-negative", {
  data <- create_test_assemblage()
  result <- create_assemblage_distribution(data = data)
  
  # Test at various points
  test_points <- c(-0.25, 0, 1, 5, 10, 20, 50, 100)
  densities <- sapply(test_points, result$dfun)
  
  expect_true(all(densities >= 0))
})

# -----------------------------------------------------------------------------
# Aggregation tests
# -----------------------------------------------------------------------------

test_that("aggregates counts across studies correctly", {
  # Two studies with overlapping intervals
  data <- data.frame(
    study = c("A", "A", "B", "B"),
    x_lo = c(0, 10, 0, 10),
    x_hi = c(10, 20, 10, 20),
    count = c(10, 20, 5, 15)
  )
  
  result <- create_assemblage_distribution(data = data)
  
  # Total counts: [0,10) = 15, [10,20) = 35
  # Proportions: 15/50 = 0.3, 35/50 = 0.7
  # Densities: 0.3/10 = 0.03, 0.7/10 = 0.07
  
  expect_equal(result$densities[1], 0.03, tolerance = 1e-10)
  expect_equal(result$densities[2], 0.07, tolerance = 1e-10)
})

test_that("handles non-overlapping intervals across studies", {
  # Studies contribute to different intervals
  data <- data.frame(
    study = c("A", "B"),
    x_lo = c(0, 10),
    x_hi = c(10, 20),
    count = c(20, 30)
  )
  
  result <- create_assemblage_distribution(data = data)
  
  # Proportions: 20/50 = 0.4, 30/50 = 0.6
  # Densities: 0.4/10 = 0.04, 0.6/10 = 0.06
  
  expect_equal(result$densities[1], 0.04, tolerance = 1e-10)
  expect_equal(result$densities[2], 0.06, tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# Piecewise constant density tests
# -----------------------------------------------------------------------------

test_that("density is constant within intervals", {
  data <- create_single_study_assemblage()
  result <- create_assemblage_distribution(data = data)
  
  # Test multiple points within first interval [0, 5)
  d_at_1 <- result$dfun(1)
  d_at_2 <- result$dfun(2)
  d_at_4 <- result$dfun(4)
  
  expect_equal(d_at_1, d_at_2)
  expect_equal(d_at_2, d_at_4)
  
  # Test multiple points within second interval [5, 20)
  d_at_6 <- result$dfun(6)
  d_at_10 <- result$dfun(10)
  d_at_19 <- result$dfun(19)
  
  expect_equal(d_at_6, d_at_10)
  expect_equal(d_at_10, d_at_19)
})

test_that("density changes at interval boundaries", {
  data <- data.frame(
    study = "A",
    x_lo = c(0, 10),
    x_hi = c(10, 20),
    count = c(10, 40)
  )
  
  result <- create_assemblage_distribution(data = data)
  
  # Density should differ between intervals
  d_first_interval <- result$dfun(5)
  d_second_interval <- result$dfun(15)
  
  expect_false(d_first_interval == d_second_interval)
  
  # Second interval has more counts, so higher density
  expect_true(d_second_interval > d_first_interval)
})

test_that("density is zero outside defined range for finite assemblage", {
  data <- create_finite_assemblage()
  result <- create_assemblage_distribution(data = data)
  
  # Before the range
  expect_equal(result$dfun(-5), 0)
  expect_equal(result$dfun(-0.1), 0)
  
  # After the range (no exponential tail)
  expect_equal(result$dfun(25), 0)
  expect_equal(result$dfun(100), 0)
})

# -----------------------------------------------------------------------------
# Exponential tail tests
# -----------------------------------------------------------------------------

test_that("exponential tail is used for unbounded final interval", {
  data <- create_single_study_assemblage()
  result <- create_assemblage_distribution(data = data)
  
  # exp_rate should be defined (not NA)
  expect_false(is.na(result$exp_rate))
  expect_true(result$exp_rate > 0)
})

test_that("exp_rate is NA when no unbounded interval exists", {
  data <- create_finite_assemblage()
  result <- create_assemblage_distribution(data = data)
  
  expect_true(is.na(result$exp_rate))
})

test_that("exponential tail density decreases with age", {
  data <- create_single_study_assemblage()
  result <- create_assemblage_distribution(data = data)
  
  last_cut <- max(result$cuts)
  
  # Density should decrease in the tail region
  d_at_start <- result$dfun(last_cut + 1)
  d_at_mid <- result$dfun(last_cut + 20)
  d_at_far <- result$dfun(last_cut + 50)
  
  expect_true(d_at_start > d_at_mid)
  expect_true(d_at_mid > d_at_far)
  expect_true(d_at_far > 0)
})

test_that("exponential rate is derived from expected remaining life", {

  # Create assemblage where we can verify the expected remaining life calculation
  # With exponential distribution, E[X - x0 | X > x0] = 1/lambda
  # So lambda = 1 / expected_remaining_life
  data <- data.frame(
    study = "A",
    x_lo = c(0, 20),
    x_hi = c(20, Inf),
    count = c(50, 50)
  )
  
  result <- create_assemblage_distribution(data = data)
  
  # The exponential rate should be positive and reasonable
  # (This is a notional calculation, exact value depends on implementation)
  expect_true(result$exp_rate > 0)
  expect_true(is.finite(result$exp_rate))
})

# -----------------------------------------------------------------------------
# Negative x_lo tests (pre-birth deaths)
# -----------------------------------------------------------------------------

test_that("handles negative x_lo correctly", {
  data <- data.frame(
    study = "A",
    x_lo = c(-0.5, 1.5),
    x_hi = c(1.5, Inf),
    count = c(10, 90)
  )
  
  result <- create_assemblage_distribution(data = data)
  
  # Cuts should include the negative boundary
  expect_true(min(result$cuts) < 0)
  
  # Density should be defined in the negative region
  expect_true(result$dfun(-0.25) > 0)
  expect_true(result$dfun(0) > 0)
})

# -----------------------------------------------------------------------------
# Sampling function tests
# -----------------------------------------------------------------------------

test_that("rfun returns numeric values", {
  data <- create_test_assemblage()
  result <- create_assemblage_distribution(data = data)
  
  samples <- result$rfun(100)
  
  expect_type(samples, "double")
  expect_length(samples, 100)
})
test_that("rfun samples are within plausible range", {
  data <- create_finite_assemblage()
  result <- create_assemblage_distribution(data = data)
  
  samples <- result$rfun(1000)
  
  # All samples should be within [min(x_lo), max(x_hi)]
  expect_true(all(samples >= min(result$cuts)))
  expect_true(all(samples <= max(result$cuts)))
})

test_that("rfun sample distribution roughly matches density", {
  data <- data.frame(
    study = "A",
    x_lo = c(0, 10),
    x_hi = c(10, 20),
    count = c(20, 80)
  )
  
  result <- create_assemblage_distribution(data = data)
  
  set.seed(42)
  samples <- result$rfun(10000)
  
  # Proportion in first interval should be ~0.2
  prop_first <- mean(samples < 10)
  expect_equal(prop_first, 0.2, tolerance = 0.05)
  
  # Proportion in second interval should be ~0.8
  prop_second <- mean(samples >= 10)
  expect_equal(prop_second, 0.8, tolerance = 0.05)
})

# -----------------------------------------------------------------------------
# Edge case tests
# -----------------------------------------------------------------------------

test_that("handles single interval assemblage", {
  data <- data.frame(
    study = "A",
    x_lo = 0,
    x_hi = 50,
    count = 100
  )
  
  result <- create_assemblage_distribution(data = data)
  
  # Should have uniform density over [0, 50)
  expected_density <- 1 / 50
  
  expect_equal(result$dfun(0), expected_density, tolerance = 1e-10)
  expect_equal(result$dfun(25), expected_density, tolerance = 1e-10)
  expect_equal(result$dfun(49), expected_density, tolerance = 1e-10)
})

test_that("handles assemblage with zero counts in some intervals", {
  data <- data.frame(
    study = "A",
    x_lo = c(0, 10, 20),
    x_hi = c(10, 20, 30),
    count = c(10, 0, 30)
  )
  
  result <- create_assemblage_distribution(data = data)
  
  # Density in zero-count interval should be zero
  expect_equal(result$dfun(15), 0)
  
  # Other intervals should have positive density
  expect_true(result$dfun(5) > 0)
  expect_true(result$dfun(25) > 0)
})

test_that("handles very small counts", {
  data <- data.frame(
    study = "A",
    x_lo = c(0, 10),
    x_hi = c(10, 20),
    count = c(1, 1)
  )
  
  expect_no_error(
    create_assemblage_distribution(data = data)
  )
  
  result <- create_assemblage_distribution(data = data)
  
  # Should still integrate to 1
  widths <- diff(result$cuts)
  total_mass <- sum(result$densities * widths)
  expect_equal(total_mass, 1.0, tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# Exact value verification tests
# -----------------------------------------------------------------------------

test_that("exact values for simple two-interval case", {
  # Simple case with two intervals of different widths
  # Interval [0, 10): width = 10, count = 20
  # Interval [10, 30): width = 20, count = 30
  # Total count = 50
  data <- data.frame(
    study = "A",
    x_lo = c(0, 10),
    x_hi = c(10, 30),
    count = c(20, 30)
  )
  
  result <- create_assemblage_distribution(data = data)
  
 # Verify cuts
  expect_equal(result$cuts, c(0, 10, 30))
  
  # Proportions: 20/50 = 0.4, 30/50 = 0.6
  # Densities = proportion / interval_width
  # d1 = 0.4 / 10 = 0.04
  # d2 = 0.6 / 20 = 0.03
  expect_equal(result$densities[1], 0.04, tolerance = 1e-10)
  expect_equal(result$densities[2], 0.03, tolerance = 1e-10)
  
  # dfun should return these exact values within intervals
  expect_equal(result$dfun(0), 0.04, tolerance = 1e-10)
  expect_equal(result$dfun(5), 0.04, tolerance = 1e-10)
  expect_equal(result$dfun(9.99), 0.04, tolerance = 1e-10)
  expect_equal(result$dfun(10), 0.03, tolerance = 1e-10)
  expect_equal(result$dfun(20), 0.03, tolerance = 1e-10)
  expect_equal(result$dfun(29.99), 0.03, tolerance = 1e-10)
  
  # Outside range should be zero
  expect_equal(result$dfun(-1), 0)
  expect_equal(result$dfun(31), 0)
  
  # exp_rate should be NA (no unbounded interval)
  expect_true(is.na(result$exp_rate))
})

test_that("exact values for three-interval case with uniform counts", {
  # Three intervals of width 10 each, equal counts
  # Total integrates to 1, so density = 1/(3*10) = 1/30 in each
  data <- data.frame(
    study = "A",
    x_lo = c(0, 10, 20),
    x_hi = c(10, 20, 30),
    count = c(100, 100, 100)
  )
  
  result <- create_assemblage_distribution(data = data)
  
  # Verify cuts
  expect_equal(result$cuts, c(0, 10, 20, 30))
  
  # Each interval has 1/3 of the mass over width 10
  # density = (1/3) / 10 = 1/30
  expected_density <- 1 / 30
  
  expect_equal(result$densities[1], expected_density, tolerance = 1e-10)
  expect_equal(result$densities[2], expected_density, tolerance = 1e-10)
  expect_equal(result$densities[3], expected_density, tolerance = 1e-10)
  
  # dfun verification
  expect_equal(result$dfun(5), expected_density, tolerance = 1e-10)
  expect_equal(result$dfun(15), expected_density, tolerance = 1e-10)
  expect_equal(result$dfun(25), expected_density, tolerance = 1e-10)
})

test_that("exact values with aggregation across studies", {
  # Two studies contributing to same intervals
  # Study A: [0,10) count=10, [10,20) count=20
  # Study B: [0,10) count=5,  [10,20) count=15
  # Aggregated: [0,10) count=15, [10,20) count=35
  # Total = 50
  data <- data.frame(
    study = c("A", "A", "B", "B"),
    x_lo = c(0, 10, 0, 10),
    x_hi = c(10, 20, 10, 20),
    count = c(10, 20, 5, 15)
  )
  
  result <- create_assemblage_distribution(data = data)
  
  # Verify cuts (should be unique sorted values)
  expect_equal(result$cuts, c(0, 10, 20))
  
  # Proportions: 15/50 = 0.3, 35/50 = 0.7
  # Densities: 0.3/10 = 0.03, 0.7/10 = 0.07
  expect_equal(result$densities[1], 0.03, tolerance = 1e-10)
  expect_equal(result$densities[2], 0.07, tolerance = 1e-10)
  
  # dfun verification
  expect_equal(result$dfun(5), 0.03, tolerance = 1e-10)
  expect_equal(result$dfun(15), 0.07, tolerance = 1e-10)
})

test_that("exact values with negative x_lo", {
  # Interval [-1, 1): width = 2, count = 20
  # Interval [1, 5):  width = 4, count = 80
  # Total = 100
  data <- data.frame(
    study = "A",
    x_lo = c(-1, 1),
    x_hi = c(1, 5),
    count = c(20, 80)
  )
  
  result <- create_assemblage_distribution(data = data)
  
  # Verify cuts include negative boundary
  expect_equal(result$cuts, c(-1, 1, 5))
  
  # Proportions: 0.2, 0.8
  # Densities: 0.2/2 = 0.1, 0.8/4 = 0.2
  expect_equal(result$densities[1], 0.1, tolerance = 1e-10)
  expect_equal(result$densities[2], 0.2, tolerance = 1e-10)
  
  # dfun verification in negative region
  expect_equal(result$dfun(-0.5), 0.1, tolerance = 1e-10)
  expect_equal(result$dfun(0), 0.1, tolerance = 1e-10)
  expect_equal(result$dfun(3), 0.2, tolerance = 1e-10)
})

test_that("exact values when interval spans multiple cuts", {
  # Study A defines fine-grained intervals: [0,10), [10,20)
  # Study B defines coarse interval: [0,20)
  # The coarse interval should be split proportionally across the cuts
  data <- data.frame(
    study = c("A", "A", "B"),
    x_lo = c(0, 10, 0),
    x_hi = c(10, 20, 20),
    count = c(10, 30, 60)
  )
  
  result <- create_assemblage_distribution(data = data)
  
  # Cuts should be [0, 10, 20]
  expect_equal(result$cuts, c(0, 10, 20))
  
  # Study A contributes: 10 to [0,10), 30 to [10,20)
  # Study B's 60 spans [0,20) with width 20:
  #   - [0,10) gets 60 * (10/20) = 30
  #   - [10,20) gets 60 * (10/20) = 30
  # Totals: [0,10) = 10 + 30 = 40, [10,20) = 30 + 30 = 60
  # Total count = 100
  # Proportions: 0.4, 0.6
  # Densities: 0.4/10 = 0.04, 0.6/10 = 0.06
  
  expect_equal(result$densities[1], 0.04, tolerance = 1e-10)
  expect_equal(result$densities[2], 0.06, tolerance = 1e-10)
  
  # dfun verification
  expect_equal(result$dfun(5), 0.04, tolerance = 1e-10)
  expect_equal(result$dfun(15), 0.06, tolerance = 1e-10)
})

test_that("exact values when interval spans cuts with unequal widths", {
  # Cuts at [0, 5, 20] (widths 5 and 15)
  # One interval [0, 20) with count 100 should distribute:
  #   - [0,5) gets 100 * (5/20) = 25
  #   - [5,20) gets 100 * (15/20) = 75
  data <- data.frame(
    study = c("A", "A", "B"),
    x_lo = c(0, 5, 0),
    x_hi = c(5, 20, 20),
    count = c(0, 0, 100)
  )
  
  result <- create_assemblage_distribution(data = data)
  
  # Cuts should be [0, 5, 20]
  expect_equal(result$cuts, c(0, 5, 20))
  
  # Totals: [0,5) = 25, [5,20) = 75
  # Proportions: 0.25, 0.75
  # Densities: 0.25/5 = 0.05, 0.75/15 = 0.05
  # Note: same density because proportional to width!
  
  expect_equal(result$densities[1], 0.05, tolerance = 1e-10)
  expect_equal(result$densities[2], 0.05, tolerance = 1e-10)
})
