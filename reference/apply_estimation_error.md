# Apply Estimation Error to Agent Population

Applies age estimation error to agents in the archaeological sample.
Simulates the measurement error inherent in osteological age estimation
methods, where uncertainty increases with age and systematic bias causes
underestimation of older individuals.

## Usage

``` r
apply_estimation_error(
  state,
  error_model = "bespoke_increasing_sd",
  error_param = c(0.1375, 2, 60, 25)
)
```

## Arguments

- state:

  A data.frame representing the current population state. Must contain
  the following columns:

  agent_id

  :   Integer. Unique identifier for each agent.

  age

  :   Numeric. True age at death for dead agents.

  lesion

  :   Logical. TRUE if the agent has a lesion.

  dead

  :   Logical. TRUE if the agent is dead, FALSE if alive.

  was_deposited

  :   Logical. TRUE if the agent was deposited.

  in_sample

  :   Logical. TRUE if the agent is in the sample.

- error_model:

  Character. The estimation error model to use. Currently only
  `"bespoke_increasing_sd"` is supported.

- error_param:

  Numeric vector containing the error parameters. For the
  `"bespoke_increasing_sd"` model, this must be a length-4 vector:

  1 sd_per_year

  :   Rate at which SD increases per year after age 20. Default
      calibration: 0.1375 gives SD ≈ 7.5 at age 60.

  2 sd_at_20

  :   SD of random error at age 20. Default: 2.

  3 bias_start

  :   Age at which systematic bias begins. Default: 60.

  4 bias_at_90

  :   Mean systematic error at age 90. Default: 25.

## Value

A data.frame with the same structure as `state`, plus a new column:

- `estimated_age`: Numeric. The estimated age with measurement error
  applied. NA for agents not in sample.

All original columns are preserved unchanged.

## Details

This function implements the final stage of archaeological filtration,
adding realistic age estimation error to recovered skeletal remains. The
error model has two components:

1.  **Random error**: Standard deviation increases with age. SD = 0 at
    age 0, rises linearly to `sd_at_20` at age 20, then continues
    increasing at `sd_per_year` for ages beyond 20.

2.  **Systematic bias**: Older individuals are systematically
    underestimated. Bias = 0 up to `bias_start`, then increases
    linearly, reaching `bias_at_90` at age 90.

Only agents who are in the sample (`in_sample = TRUE`) receive estimated
ages. Agents not in the sample have `estimated_age = NA`.

The true age (`age`) is preserved; the estimated age is stored in a new
column `estimated_age`.

## See also

[`apply_preservation`](https://amy-s-anderson.github.io/persephone/reference/apply_preservation.md)
for the preservation filtering module that precedes estimation error.

## Examples

``` r
# Create a test population of sampled individuals
state <- data.frame(
  agent_id = 1:10,
  age = c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90),
  lesion = rep(FALSE, 10),
  dead = rep(TRUE, 10),
  was_deposited = rep(TRUE, 10),
  in_sample = rep(TRUE, 10)
)

# Default error parameters
error_param <- c(0.1375, 2, 60, 25)

# Apply estimation error
set.seed(42)
result <- apply_estimation_error(
  state = state,
  error_model = "bespoke_increasing_sd",
  error_param = error_param
)
```
