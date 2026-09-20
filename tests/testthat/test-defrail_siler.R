# ============================================================================
# Tests for defrail_siler()
# ============================================================================
#
# HOW THIS FILE IS ORGANIZED
# --------------------------
# testthat files live in tests/testthat/ and are named test-<something>.R.
# Each test_that() block is one independent test: it should test ONE idea,
# have a clear description (this shows up in test output when it fails),
# and use expect_*() calls to make assertions. If an expect_*() fails,
# testthat reports it and moves to the next test_that() block -- one
# failure doesn't stop the whole file from running, which is why you want
# many small test_that() blocks rather than one giant one.
#
# A FIXTURE is just test data you set up once and reuse across multiple
# tests. Below, `test_regime` is a fixture: a made-up (but realistic-shaped)
# Siler mortality regime, small and easy to reason about by hand, so that
# when a test fails you can actually check the numbers yourself rather than
# trusting a huge real dataset.
#
# A NOTE ON TOLERANCE: almost none of these tests use exact equality
# (expect_identical / expect_equal with no tolerance), because defrail_siler
# involves numerical integration (RK4). Two numerically-integrated values
# that are "the same" mathematically will differ in their 15th decimal
# place, and sometimes more, depending on step size. expect_equal()'s
# `tolerance` argument lets you say "close enough" instead of "bit-for-bit
# identical" -- get in the habit of choosing tolerance deliberately (and
# writing a comment explaining why you chose it), rather than copy-pasting
# a default.
# ============================================================================

test_regime <- data.frame(
  a1 = 0.4,     # infant mortality component
  b1 = 1.2,     # rate of infant mortality decline
  a2 = 0.01,    # age-independent ("background") hazard
  a3 = 0.0002,  # senescent mortality component
  b3 = 0.08     # rate of senescent mortality increase
)
# a3 < 1 here on purpose, so compute_siler_risk() takes the "traditional
# Siler parameters" branch rather than the "robust Siler parameters" branch.
# If you also want to test the a3 > 1 (robust parameterization) branch,
# duplicate the fixture with a3 > 1 and re-run the relevant tests against it.


# ----------------------------------------------------------------------------
# 1. STRUCTURAL TESTS
# ----------------------------------------------------------------------------
# These don't check that the *numbers* are right -- they check that the
# function's contract (what it promises to return) holds. Structural tests
# are cheap to write and catch a surprising number of real bugs (wrong
# column names, off-by-one row counts, etc.) that would otherwise only
# surface later as a cryptic error somewhere downstream in the pipeline.

test_that("defrail_siler returns a data frame with the documented structure", {
  result <- defrail_siler(test_regime, frailty_variance = 0.01, max_age = 50)
  
  expect_s3_class(result, "data.frame")
  expect_named(result, c("age", "mu0"))
  expect_type(result$age, "integer")
  expect_type(result$mu0, "double")
  # ages should be a clean, gapless integer sequence starting at 0 --
  # anything else would break match(pop$age, mu0_lookup$age) calls
  # downstream in compute_cause_specific_hazards().
  expect_identical(result$age, 0:min(50, 110))
})

test_that("defrail_siler respects max_age as a hard cap on returned rows", {
  result <- defrail_siler(test_regime, frailty_variance = 0.2,
                          max_age = 80, max_survivable_age = 110)
  expect_equal(max(result$age), 80)
  expect_equal(nrow(result), 81)  # ages 0:80 inclusive
})


test_that("the guard evaluates survivorship at max_survivable_age regardless of max_age", {
  flat_regime <- data.frame(a1 = 0, b1 = 1, a2 = 0.0001, a3 = 0, b3 = 0.001)
  
  # Same implausible scenario, checked at three different max_age values --
  # all three should behave identically, since max_age should never affect
  # what the guard actually checks.
  for (ma in c(50, 110, 200)) {
    expect_error(
      defrail_siler(flat_regime, frailty_variance = 50,
                    max_age = ma, max_survivable_age = 110),
      "survive to age 110"
    )
  }
})

test_that("mu0 is non-negative and finite at every age", {
  # A hazard can never legitimately be negative, and RK4 blowing up
  # (e.g. due to a bad step size or a pathological mortality_regime)
  # would show up as Inf/NaN before it showed up as anything else.
  result <- defrail_siler(test_regime, frailty_variance = 1, max_age = 100)
  expect_true(all(is.finite(result$mu0)))
  expect_true(all(result$mu0 >= 0))
})


test_that("the ODE's mu0 at the exact point a=0 equals the raw Siler hazard at a=0", {
  # At a = 0, H0 = 0 by the initial condition -- no selection has happened
  # yet, so the correction factor is exactly 1 and mu0(0) should equal
  # mu_bar(0) exactly, with no approximation. This is the true boundary
  # condition; note this is NOT the same as defrail_siler()'s age == 0
  # output row, which reflects the midpoint-offset value at age 0.5, not
  # this exact point.
  ode <- solve_defrail_ode(test_regime, s2 = 0.2,
                           max_survivable_age = 110, step = 0.01)
  expect_equal(ode$mu0_fine[1], compute_siler_risk(0, test_regime), tolerance = 0)
})


# ----------------------------------------------------------------------------
# 2. INPUT VALIDATION TESTS
# ----------------------------------------------------------------------------
# These confirm the function fails LOUDLY and CLEARLY on bad input, rather
# than silently returning nonsense. This matters a lot for a function whose
# output feeds into a stochastic simulation -- a silent garbage-in/
# garbage-out failure here could masquerade as a "surprising" simulation
# result for a long time before anyone suspects the lookup table itself.

test_that("negative frailty_variance is rejected with an informative error", {
  expect_error(
    defrail_siler(test_regime, frailty_variance = -0.1),
    "non-negative"
  )
})

test_that("frailty_variance = 0 and frailty_variance = NULL behave identically", {
  # The docstring says NULL should be treated the same as 0. This is easy
  # to accidentally break if someone edits the `if` condition later
  # (e.g. changes `||` to `&&`, or reorders the null check) -- an explicit
  # test locks the documented behavior in place.
  result_zero <- defrail_siler(test_regime, frailty_variance = 0, max_age = 50)
  result_null <- defrail_siler(test_regime, frailty_variance = NULL, max_age = 50)
  expect_equal(result_zero, result_null)
})

test_that("an implausibly high frailty_variance triggers the survivorship guard", {
  # defrail_siler has a built-in sanity check: if the defrailed hazard
  # would let more than 0.1% of the population survive past
  # max_survivable_age, it stops rather than silently producing a lookup
  # table that would populate your ABM with Methuselah-aged agents. This
  # test deliberately tries to trigger that guard using an extremely low
  # mortality regime + an enormous frailty variance, and checks that the
  # function actually stops as documented rather than, say, producing
  # NaN or an infinite loop.
  flat_regime <- data.frame(a1 = 0, b1 = 1, a2 = 0.0001, a3 = 0, b3 = 0.001)
  expect_error(
    defrail_siler(flat_regime, frailty_variance = 50),
    "survive to age"
  )
})


# ----------------------------------------------------------------------------
# 3. THE REGRESSION TEST -- the one that would have caught the actual bug
# ----------------------------------------------------------------------------
# This is the most important test in this file, so it gets extra
# explanation. The bug you just found was: the frailty_variance == 0
# branch and the frailty_variance > 0 (ODE) branch used two SEPARATE code
# paths, and those paths drifted out of sync (one applied a midpoint age
# offset, the other didn't). Nothing about either branch was "wrong" in
# isolation -- the bug only existed in the SEAM between them.
#
# The mathematical fact that exposes this seam: as frailty_variance -> 0,
# the ODE solution must smoothly approach the raw (undefrailed) hazard.
# There's no discontinuity in the underlying math at s2 = 0 -- the
# correction term exp(s2 * Lambda(a)) approaches exp(0) = 1 continuously.
# So defrail_siler(regime, 0) and defrail_siler(regime, 1e-6) should be
# ALMOST indistinguishable. If they're not, one of the two branches is
# using a different convention than the other -- which is exactly what
# happened here.
#
# This test is deliberately written to not assume anything about which
# age convention (raw age vs. midpoint-offset age) is "correct" -- it only
# asserts that both branches agree with EACH OTHER in the boundary limit.
# That makes it a durable regression test: it stays valid even if you
# later change the offset convention, as long as you remember to apply
# the change consistently to both branches.

test_that("defrail_siler is continuous at the frailty_variance = 0 boundary", {
  result_exact <- defrail_siler(test_regime, frailty_variance = 0, max_age = 60)
  result_tiny  <- defrail_siler(test_regime, frailty_variance = 1e-6, max_age = 60)
  
  # tolerance = 1e-3 here is deliberately loose-ish: frailty_variance = 1e-6
  # is not IDENTICAL to 0, so a tiny, genuine difference is expected and
  # fine. What this test rules out is a LARGE, systematic difference (like
  # the ~40% gap you found) that indicates the two branches disagree about
  # something structural (like an age offset), not just about numerical
  # precision.
  expect_equal(result_tiny$mu0, result_exact$mu0, tolerance = 1e-3)
})


# ----------------------------------------------------------------------------
# 4. CORRECTNESS AGAINST THE CLOSED-FORM SOLUTION
# ----------------------------------------------------------------------------
# For gamma frailty, there's an exact analytic solution to the defrailing
# ODE (this is the same relationship the docstring's RK4 approach is
# numerically approximating):
#
#     mu0(a) = mu_bar(a) * exp(s2 * Lambda(a))
#
# where Lambda(a) = integral of mu_bar from 0 to a (the target/observed
# cumulative hazard). This test computes that closed form independently
# (via a simple numerical integral, NOT via defrail_siler's own RK4 loop --
# it would be circular to check the function against itself), and confirms
# the RK4 output tracks it closely. This is the test that actually verifies
# defrail_siler's core numerical method is correct, as opposed to just
# internally consistent.
#

# closed form of cumulative hazard
reference_mu0 <- function(ages, mortality_regime, s2, integration_step = 1e-4) {
  vapply(ages, function(a) {
    if (a <= 0) return(compute_siler_risk(0, mortality_regime))
    fine   <- seq(0, a, by = integration_step)
    haz    <- compute_siler_risk(fine, mortality_regime)
    Lambda <- sum(haz[-length(haz)]) * integration_step  # left-Riemann approx of the integral
    mu_bar <- compute_siler_risk(a, mortality_regime)
    mu_bar * exp(s2 * Lambda)
  }, numeric(1))
}

test_that("defrail_siler matches the closed-form gamma-frailty solution", {
  s2 <- 0.4
  mu0_tab <- defrail_siler(test_regime, frailty_variance = s2, max_age = 40)
  
  # Restrict to a handful of ages rather than all 40, so the test runs fast
  # -- the reference function does its own numerical integration per age,
  # which is O(n) per age and would be slow across every single age.
  check_ages <- c(0, 5, 10, 20, 30, 40)
  ages_for_reference <- check_ages + 0.5  # + 0.5 because the lookup table calculates the hazard for the midpoint between annual ages
  
  expected <- reference_mu0(ages_for_reference, test_regime, s2)
  actual   <- mu0_tab$mu0[mu0_tab$age %in% check_ages]
  
  # tolerance = 1e-2 (1%) is fairly loose because the reference integral
  # itself is only a rough left-Riemann-sum approximation, not a
  # high-precision integrator. This test is meant to catch a genuinely
  # wrong implementation (wrong exponent, wrong sign, swapped shape/rate,
  # etc.), not to certify RK4's numerical precision to many decimal places.
  expect_equal(actual, expected, tolerance = 1e-2)
})


# ----------------------------------------------------------------------------
# 5. DIRECTIONAL / MONOTONICITY TESTS
# ----------------------------------------------------------------------------
# These encode qualitative facts you know must be true from the theory,
# even without computing exact numbers. They're valuable because they're
# robust to small implementation changes (e.g. switching integrators) --
# they'll only fail if something is fundamentally backwards.

test_that("higher frailty_variance produces a higher defrailed hazard at later ages", {
  # Intuition: more variance in frailty means stronger selection over time
  # (the low-frailty survivors increasingly dominate the observed
  # population), so the TRUE individual-level hazard needed to reproduce
  # the same population-level Siler curve must be pushed higher as age
  # (and variance) increase, to compensate for that selection.
  low  <- defrail_siler(test_regime, frailty_variance = 0.2, max_age = 60)
  high <- defrail_siler(test_regime, frailty_variance = 1.0, max_age = 60)
  
  # Check this at an age well past 0 (selection needs time to accumulate;
  # right at birth there's essentially no difference -- see the boundary
  # test above).
  older_ages <- low$age >= 20
  expect_true(all(high$mu0[older_ages] > low$mu0[older_ages]))
})

test_that("defrailed hazard is non-decreasing in its own correction relative to raw Siler", {
  # For any frailty_variance > 0, the defrailed hazard should be >= the
  # raw (undefrailed) Siler hazard at every age past 0 -- the correction
  # factor exp(s2 * Lambda(a)) is always >= 1 for s2, Lambda(a) >= 0.
  raw       <- compute_siler_risk(0:60 + 0.5, test_regime)
  defrailed <- defrail_siler(test_regime, frailty_variance = 0.5, max_age = 60)$mu0
  expect_true(all(defrailed >= raw - 1e-8))  # small epsilon for floating-point comparison
})

# Quick plot to check the relationship between raw and defrailed Siler values.
# library(ggplot2)
# age <- 0:60
# raw_vs_defrailed <- data.frame(age = age, raw = raw, defrailed = defrailed)
# ggplot(raw_vs_defrailed, aes(x = age, y = raw)) +
#   geom_line() +
#   geom_line(aes(x = age, y = defrailed, color = "red"))



# ----------------------------------------------------------------------------
# 6. NUMERICAL STABILITY TEST
# ----------------------------------------------------------------------------
# This checks that the RK4 integration has actually converged -- i.e. that
# step = 0.01 (the default) is fine enough that halving it further doesn't
# meaningfully change the answer. If this test ever starts failing (e.g.
# after changing max_age, or using a much more extreme mortality_regime),
# it's a sign the default step size needs revisiting for that use case.

test_that("results are stable across reasonable step size choices", {
  coarse <- defrail_siler(test_regime, frailty_variance = 0.5, max_age = 60, step = 0.01)
  fine   <- defrail_siler(test_regime, frailty_variance = 0.5, max_age = 60, step = 0.001)
  
  expect_equal(coarse$mu0, fine$mu0, tolerance = 1e-3)
})
