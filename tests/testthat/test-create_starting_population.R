# ============================================================================
# Tests for create_pop()
# ============================================================================
#
# SCOPING NOTE: create_pop() branches heavily on age_structured. When TRUE,
# it delegates most of the real work to create_pop_stable_age(), which we
# don't have visibility into here. So these tests focus mostly on the
# age_structured = FALSE (cohort) path, where create_pop()'s own logic is
# fully self-contained and directly testable -- which conveniently is also
# the path your frailty-variance experiments actually use. The
# age_structured = TRUE tests below are deliberately light: they check
# create_pop()'s own input-validation behavior, not create_pop_stable_age()'s
# internals (that belongs in create_pop_stable_age()'s own test file).
#
# A NOTE ON pop_config: this is a plain list, not a validated object, so
# create_pop() will happily accept a pop_config missing a field it expects
# (treating it as NULL) rather than erroring. That's convenient for testing
# -- each fixture below only sets the fields relevant to what it's testing
# -- but also means create_pop() can't tell the difference between "this
# field was deliberately left unset" and "the caller forgot it." Keep that
# in mind if you ever see create_pop() behave as if a feature is "off" when
# you meant to turn it on -- check the pop_config key name and spelling
# first (this is literally what happened upstream in Simulate_Cemetery's
# own pop_config construction -- see the conversation this file came from).
# ============================================================================

# A minimal, fully-specified pop_config for the cohort (age_structured =
# FALSE) case. Individual tests below override just the field(s) they care
# about, using modifyList(), so each test's intent is visible at a glance
# without repeating the whole fixture.
pop_config_minimal <- list(
  model_lesions           = FALSE,
  annual_exposure         = NULL,
  lesion_formation_rate   = NULL,
  lesion_formation_window = c(0, 0),
  frailty_variance        = NULL
)


# ----------------------------------------------------------------------------
# 1. STRUCTURAL TESTS -- the basic contract, independent of any optional
#    features (lesions, frailty, exposure) being turned on
# ----------------------------------------------------------------------------

test_that("a cohort population has one row per agent, all starting at age 0", {
  pop <- create_pop(pop0_size = 100, age_structured = FALSE,
                    pop_config = pop_config_minimal)
  
  expect_equal(nrow(pop), 100)
  expect_true(all(pop$age == 0))
  # agent_id should uniquely identify every agent -- downstream code (e.g.
  # decedent record-keeping) likely relies on this being a clean 1:n key.
  expect_equal(pop$agent_id, 1:100)
  expect_equal(length(unique(pop$agent_id)), 100)
})

test_that("with no optional features enabled, only agent_id and age are created", {
  # This is a useful baseline: it documents exactly which columns are
  # "extra" (lesion, frailty, acquired_frailty, n_stress_events) versus
  # which are always present -- so if a future change accidentally adds
  # a column unconditionally, this test will catch it.
  pop <- create_pop(pop0_size = 10, age_structured = FALSE,
                    pop_config = pop_config_minimal)
  expect_named(pop, c("agent_id", "age"))
})


# ----------------------------------------------------------------------------
# 2. LESION COLUMN TESTS
# ----------------------------------------------------------------------------

test_that("model_lesions = TRUE adds a lesion column positioned right after age", {
  cfg <- modifyList(pop_config_minimal, list(
    model_lesions           = TRUE,
    lesion_formation_window = c(0, 5)
  ))
  pop <- create_pop(pop0_size = 10, age_structured = FALSE, pop_config = cfg)
  
  expect_true("lesion" %in% names(pop))
  # relocate() is supposed to place lesion immediately after age -- check
  # the actual column ORDER, not just presence, since downstream code or
  # documentation might assume this position.
  expect_equal(which(names(pop) == "lesion"), which(names(pop) == "age") + 1)
})

test_that("agents within the lesion formation window get lesion = 0 (not yet lesioned, but eligible)", {
  # In a cohort, everyone starts at age 0. If the window includes age 0,
  # everyone should start with lesion = 0 (eligible, not-yet-formed) --
  # NOT NA, which is reserved for "outside the window."
  cfg <- modifyList(pop_config_minimal, list(
    model_lesions           = TRUE,
    lesion_formation_window = c(0, 5)   # includes age 0
  ))
  pop <- create_pop(pop0_size = 10, age_structured = FALSE, pop_config = cfg)
  expect_true(all(pop$lesion == 0))
})

test_that("agents outside the lesion formation window get lesion = NA", {
  # Same cohort (everyone at age 0), but now the window doesn't open until
  # age 1 -- so nobody is currently eligible, and lesion should be NA for
  # everyone rather than 0 (0 would incorrectly imply "eligible, hasn't
  # formed one," which isn't true at age 0 if the window hasn't opened).
  cfg <- modifyList(pop_config_minimal, list(
    model_lesions           = TRUE,
    lesion_formation_window = c(1, 5)   # excludes age 0
  ))
  pop <- create_pop(pop0_size = 10, age_structured = FALSE, pop_config = cfg)
  expect_true(all(is.na(pop$lesion)))
})

test_that("model_lesions = FALSE never creates a lesion column, regardless of window", {
  # Guards against a future edit accidentally hoisting lesion-column
  # creation out from under the model_lesions check.
  cfg <- modifyList(pop_config_minimal, list(
    model_lesions           = FALSE,
    lesion_formation_window = c(0, 5)
  ))
  pop <- create_pop(pop0_size = 10, age_structured = FALSE, pop_config = cfg)
  expect_false("lesion" %in% names(pop))
})


# ----------------------------------------------------------------------------
# 3. FRAILTY COLUMN TESTS
# ----------------------------------------------------------------------------
# These are the most important tests in this file, since frailty
# initialization is exactly where the bug in this conversation lived.

test_that("frailty_variance = NULL creates neither frailty nor acquired_frailty", {
  cfg <- modifyList(pop_config_minimal, list(frailty_variance = NULL))
  pop <- create_pop(pop0_size = 10, age_structured = FALSE, pop_config = cfg)
  
  expect_false("frailty" %in% names(pop))
  expect_false("acquired_frailty" %in% names(pop))
})

test_that("frailty_variance = 0 gives every agent frailty = 1 exactly", {
  cfg <- modifyList(pop_config_minimal, list(frailty_variance = 0))
  pop <- create_pop(pop0_size = 10, age_structured = FALSE, pop_config = cfg)
  
  expect_true("frailty" %in% names(pop))
  expect_true(all(pop$frailty == 1))
})

test_that("frailty_variance > 0 draws from a Gamma distribution with mean 1 and the specified variance", {
  # Statistical test: can't check individual draws exactly (they're
  # random), but with a large enough n, the SAMPLE mean and variance
  # should land close to the theoretical mean (1) and variance (s2).
  # Setting a seed makes this test deterministic -- without it, this test
  # would occasionally fail by chance alone (a "flaky" test), which
  # erodes trust in your test suite over time. A generous tolerance
  # (rather than a tiny one) is intentional: it should be loose enough
  # that ordinary sampling variability at this n doesn't cause spurious
  # failures, but tight enough to catch a real error (e.g. a swapped
  # shape/rate parameterization, which would give a wildly different
  # mean or variance, not just a slightly-off one).
  set.seed(42)
  s2 <- 0.3
  cfg <- modifyList(pop_config_minimal, list(frailty_variance = s2))
  pop <- create_pop(pop0_size = 100000, age_structured = FALSE, pop_config = cfg)
  
  expect_equal(mean(pop$frailty), 1,  tolerance = 0.02)
  expect_equal(var(pop$frailty),  s2, tolerance = 0.02)
})

test_that("frailty is never negative, for any variance", {
  # Sanity check on the distributional family itself: Gamma is only
  # defined on non-negative values, so this should always hold by
  # construction -- but it's cheap insurance against, e.g., someone later
  # swapping rgamma() for a different distribution that doesn't share
  # this property, without updating this assumption elsewhere in the
  # codebase.
  set.seed(1)
  cfg <- modifyList(pop_config_minimal, list(frailty_variance = 2))
  pop <- create_pop(pop0_size = 10000, age_structured = FALSE, pop_config = cfg)
  expect_true(all(pop$frailty >= 0))
})

test_that("REGRESSION: acquired_frailty is not created without exposure/lesion mechanism active, even if frailty is non-null", {
  cfg <- modifyList(pop_config_minimal, list(
    frailty_variance      = 0,     # birth frailty variance set...
    annual_exposure       = NULL,  # ...but no exposure...
    lesion_formation_rate = NULL   # ...and no lesion mechanism either
  ))
  pop <- create_pop(pop0_size = 10, age_structured = FALSE, pop_config = cfg)
  expect_false("acquired_frailty" %in% names(pop))
})


# ----------------------------------------------------------------------------
# 4. n_stress_events COLUMN TESTS
# ----------------------------------------------------------------------------

test_that("n_stress_events is created when annual_exposure is set", {
  cfg <- modifyList(pop_config_minimal, list(annual_exposure = 0.1))
  pop <- create_pop(pop0_size = 10, age_structured = FALSE, pop_config = cfg)
  
  expect_true("n_stress_events" %in% names(pop))
  expect_true(all(pop$n_stress_events == 0L))
  expect_type(pop$n_stress_events, "integer")
})

test_that("n_stress_events is created when lesion_formation_rate is set", {
  # Tested with lesion_formation_rate specifically because -- as flagged
  # in the surrounding conversation -- Simulate_Cemetery's own pop_config
  # construction currently never sets this key at all, meaning this path
  # is effectively dead code when called through Simulate_Cemetery, even
  # though create_pop()'s own logic correctly handles it. This test
  # confirms create_pop() itself is fine; the caller-side gap belongs in
  # a Simulate_Cemetery test file.
  cfg <- modifyList(pop_config_minimal, list(lesion_formation_rate = 0.05))
  pop <- create_pop(pop0_size = 10, age_structured = FALSE, pop_config = cfg)
  
  expect_true("n_stress_events" %in% names(pop))
  expect_true(all(pop$n_stress_events == 0L))
})

test_that("n_stress_events is absent when neither exposure mechanism is set", {
  pop <- create_pop(pop0_size = 10, age_structured = FALSE,
                    pop_config = pop_config_minimal)
  expect_false("n_stress_events" %in% names(pop))
})


# ----------------------------------------------------------------------------
# 5. age_structured = TRUE: input validation only
# ----------------------------------------------------------------------------
# These deliberately don't check create_pop_stable_age()'s actual output --
# that function needs its own test file. Here we're only checking
# create_pop()'s own responsibility: handling (or failing to handle) a
# missing mortality_regime.

test_that("age_structured = TRUE with mortality_regime = NULL does not silently succeed", {
  expect_error(
    create_pop(pop0_size = 10, age_structured = TRUE,
               pop_config = pop_config_minimal, mortality_regime = NULL)
  )
})
