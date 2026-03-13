# Apply Mortality to Agent Population

Applies a single time step of mortality to a population of agents using
the Usher3 illness-death model. Agents with lesions experience elevated
mortality risk compared to agents without lesions.

## Usage

``` r
apply_mortality_usher3(state, mortality_model, mortality_param, dx)
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

- mortality_model:

  Character. The mortality model to use. Currently only `"usher3"` is
  supported. Any other value will raise an error.

- mortality_param:

  Numeric vector containing the mortality parameters. For the `"usher3"`
  model, this must be a length-6 vector:

  `[1]` k2

  :   The mortality multiplier for agents with lesions. A value of 1.0
      means no excess mortality; values \> 1 indicate elevated risk.
      Must be non-negative.

  `[2:6]` b_siler

  :   The five Siler hazard parameters using the demohaz
      parameterization: `c(b1, b2, b3, b4, b5)`. See
      [`hsiler`](https://rdrr.io/pkg/demohaz/man/hsiler.html) for
      details.

- dx:

  Numeric. The time step size in years. Determines both the mortality
  probability (hazard \* dx) and the age increment for survivors.
  Smaller values give more accurate results but require more iterations.

## Value

A data.frame with the same structure as `state`, with updated values:

- `dead`: Updated to TRUE for agents who died this time step.

- `age`: Incremented by `dx` for surviving agents; unchanged for agents
  who died (preserving age at death).

- `in_sample`: Unchanged (passed through from input).

## Details

This function implements a discrete-time approximation of the Usher3
illness-death model for mortality. Within each time step of size `dx`,
each living agent faces a probability of death equal to `hazard * dx`,
where the hazard depends on the agent's lesion status:

- Agents without lesions: `p_die = hsiler(age, b_siler) * dx`

- Agents with lesions: `p_die = k2 * hsiler(age, b_siler) * dx`

This approximation is valid when `hazard * dx << 1`. For accuracy, use
small values of `dx` (e.g., 0.01 to 1).

Agents who die retain their age at death. Surviving agents have their
age incremented by `dx`.

Note: This module handles mortality only. Lesion acquisition and age
filtration (archaeological sampling bias) are handled by separate
modules.

## See also

[`hsiler`](https://rdrr.io/pkg/demohaz/man/hsiler.html) for the Siler
hazard function.

## Examples

``` r
# Create a small test population
state <- data.frame(
  agent_id = 1:10,
  age = c(0, 5, 10, 20, 30, 40, 50, 60, 70, 80),
  lesion = c(FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE),
  dead = rep(FALSE, 10),
  in_sample = rep(TRUE, 10)
)

# Siler parameters (Gage & Dyke 1986, Table 2, Level 15)
a_siler <- c(0.175, 1.40, 0.00368, 0.000075, 0.0917)
b_siler <- demohaz::trad_to_demohaz_siler_param(a_siler)
mortality_param <- c(1.2, b_siler)  # k2 = 1.2

# Apply one year of mortality
result <- apply_mortality_usher3(
  state = state,
  mortality_model = "usher3",
  mortality_param = mortality_param,
  dx = 1
)
```
