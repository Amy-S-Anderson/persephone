# Apply Deposition Filtering to Agent Population

Applies deposition filtering to a population of dead agents. Determines
which deceased individuals are deposited into the archaeological record
based on age cutoff criteria.

## Usage

``` r
apply_deposition(state, deposition_model, deposition_param, dx)
```

## Arguments

- state:

  A data.frame representing the current population state. Must contain
  the following columns:

  agent_id

  :   Integer. Unique identifier for each agent.

  age

  :   Numeric. Age at death for dead agents.

  lesion

  :   Logical. TRUE if the agent has a lesion.

  dead

  :   Logical. TRUE if the agent is dead, FALSE if alive.

  was_deposited

  :   Logical. Will be updated by this function.

  in_sample

  :   Logical. TRUE if the agent is in the sample.

- deposition_model:

  Character. The deposition model to use. Currently only `"cutoff"` is
  supported.

- deposition_param:

  Numeric vector containing the deposition parameters. For the
  `"cutoff"` model, this must be length 1:

  1 cutoff_age

  :   The minimum age for deposition. Dead agents with age \>=
      cutoff_age are deposited; those below are not.

- dx:

  Numeric. The time step size in years. Not used in the current
  implementation but included for interface consistency with other
  modules.

## Value

A data.frame with the same structure as `state`, with updated values:

- `was_deposited`: Set to TRUE for dead agents at or above the cutoff
  age; FALSE otherwise.

- All other columns: Unchanged.

## Details

This function implements the first stage of archaeological filtration.
Dead agents whose age at death meets the deposition criteria are marked
as deposited (`was_deposited = TRUE`).

Currently, only the `"cutoff"` model is supported, which deposits all
dead agents at or above a minimum age. This models scenarios where very
young individuals (e.g., neonates) are less likely to enter the
archaeological record due to burial practices or taphonomic factors.

The deposition decision is deterministic: all dead agents at or above
the cutoff age are deposited; all below are not.

Note: This module handles deposition only. Living agents are not
processed. Preservation filtering (whether deposited remains are
recoverable) is handled by a separate module.

## See also

[`apply_preservation`](https://amy-s-anderson.github.io/persephone/reference/apply_preservation.md)
for the preservation filtering module that follows deposition.

## Examples

``` r
# Create a test population of deceased individuals
state <- data.frame(
  agent_id = 1:10,
  age = c(0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7),
  lesion = rep(FALSE, 10),
  dead = rep(TRUE, 10),
  was_deposited = rep(FALSE, 10),
  in_sample = rep(TRUE, 10)
)

# Apply deposition with cutoff at age 2.5 (exclude infants/toddlers)
result <- apply_deposition(
  state = state,
  deposition_model = "cutoff",
  deposition_param = c(2.5),
  dx = 1
)
```
