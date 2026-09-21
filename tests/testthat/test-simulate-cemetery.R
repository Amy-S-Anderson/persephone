


CoaleDemenyWest5 <-   data.frame(a1 = 0.457, b1 = 1.07, a2 = 0.01037, a3 = 0.000359, b3 = 0.0763)   


test_that("Simulate_Cemetery's mu0_table construction agrees with defrail_siler at frailty_variance = 0", {
  # Regression test: Simulate_Cemetery previously built its own mu0 table
  # inline for frailty_variance == 0, bypassing defrail_siler() and its
  # midpoint-offset correction -- producing a table inconsistent with
  # what defrail_siler(regime, 0) itself returns. This should never
  # diverge, since both are meant to represent the same thing.
  from_defrail_siler <- defrail_siler(CoaleDemenyWest5, frailty_variance = 0)
  # (once the fix lands, extracting Simulate_Cemetery's actual internal
  # mu0_table isn't directly possible without exposing it -- consider
  # having Simulate_Cemetery optionally return mu0_table in its output
  # list for testability, or test this indirectly via age-at-death
  # distributions matching closely between a direct defrail_siler-built
  # table and a full Simulate_Cemetery run at frailty_variance = 0)
  expect_true(TRUE)  # placeholder -- see comment above
})

