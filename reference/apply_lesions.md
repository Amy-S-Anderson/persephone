# Apply Lesion Acquisition to Agent Population

Applies a single time step of lesion acquisition to a population of
agents. Agents without lesions may acquire one based on a constant
transition rate, optionally restricted to an age window.

## Usage

``` r
apply_lesions(state, lesion_model, lesion_param, dx)
```

## Arguments

- state:

  A data.frame representing the current population state. Must contain
  the following columns:

  agent_id

  :   Integer. Unique identifier for each agent.

  age

  :   Numeric. Current age of each agent.

  lesion

  :   Logical. TRUE if the agent has a lesion, FALSE otherwise.

  dead

  :   Logical. TRUE if the agent is dead, FALSE if alive.

  in_sample

  :   Logical. TRUE if the agent is in the sample.

- lesion_model:

  Character. The lesion acquisition model to use. One of:

  `"constant"`

  :   Constant rate, always possible (no age restriction)

  `"constant_to"`

  :   Constant rate up to a specified age

  `"constant_from"`

  :   Constant rate after a specified age

  `"constant_interval"`

  :   Constant rate within an age interval

- lesion_param:

  Numeric vector containing the lesion parameters. The required length
  depends on the model:

  `"constant"`

  :   Length 1: `c(k1)`

  `"constant_to"`

  :   Length 2: `c(k1, age_end)`

  `"constant_from"`

  :   Length 2: `c(k1, age_start)`

  `"constant_interval"`

  :   Length 3: `c(k1, age_start, age_end)`

  where `k1` is the transition rate (hazard coefficient, not
  probability).

- dx:

  Numeric. The time step size in years. Determines the lesion
  acquisition probability (`k1 * dx`).

## Value

A data.frame with the same structure as `state`, with updated values:

- `lesion`: Updated to TRUE for agents who acquired a lesion.

- All other columns: Unchanged (age is not modified).

## Details

This function implements a discrete-time approximation for lesion
acquisition. Within each time step of size `dx`, each eligible agent
faces a probability of acquiring a lesion equal to `k1 * dx`, where `k1`
is the transition rate (hazard coefficient).

An agent is eligible for lesion acquisition if:

- The agent is alive (`dead == FALSE`)

- The agent does not already have a lesion (`lesion == FALSE`)

- The agent's age falls within the model's age window

This approximation is valid when `k1 * dx << 1`. For accuracy, use small
values of `dx` (e.g., 0.01 to 1).

Note: This module handles lesion acquisition only. Mortality and age
increments are handled by separate modules.

## See also

`apply_mortality` for the mortality module.

## Examples

``` r
# Create a small test population
state <- data.frame(
  agent_id = 1:10,
  age = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
  lesion = rep(FALSE, 10),
  dead = rep(FALSE, 10),
  in_sample = rep(TRUE, 10)
)

# Apply lesion acquisition with window [0, 6)
result <- apply_lesions(
  state = state,
  lesion_model = "constant_to",
  lesion_param = c(0.05, 6),  # k1 = 0.05, window ends at age 6
  dx = 1
)
```
