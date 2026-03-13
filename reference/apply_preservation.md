# Apply Preservation Filtering to Agent Population

Applies preservation filtering to deposited agents. Determines which
deposited skeletal remains are preserved and recoverable based on an
age-dependent Siler hazard model.

## Usage

``` r
apply_preservation(state, preservation_model, preservation_param, dx)
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

  :   Logical. TRUE if the agent was deposited.

  in_sample

  :   Logical. Will be updated by this function.

- preservation_model:

  Character. The preservation model to use. Currently only `"siler"` is
  supported.

- preservation_param:

  Numeric vector containing the preservation parameters. For the
  `"siler"` model, this must be a length-5 vector containing the Siler
  hazard parameters in the demohaz parameterization:
  `c(b1, b2, b3, b4, b5)`. See
  [`hsiler`](https://rdrr.io/pkg/demohaz/man/hsiler.html) for details.

- dx:

  Numeric. The time step size in years. Determines the filtering
  probability (`hsiler(age) * dx`).

## Value

A data.frame with the same structure as `state`, with updated values:

- `in_sample`: Set to FALSE for non-deposited agents and for deposited
  agents who fail the Siler filter; unchanged for deposited agents who
  pass the filter.

- All other columns: Unchanged.

## Details

This function implements the second stage of archaeological filtration.
Deposited agents are subject to a stochastic filtering process based on
the Siler hazard function, which can model age-dependent preservation
biases in the archaeological record.

The filtering uses a discrete-time approximation: for each deposited
agent, the probability of being filtered out (not preserved) is
`hsiler(age) * dx`. This approximation is valid when `hazard * dx << 1`.

Agents who are not deposited (`was_deposited = FALSE`) automatically
have `in_sample` set to FALSE, as they cannot be recovered if they were
never deposited.

Agents who are already `in_sample = FALSE` remain FALSE (no
re-sampling).

Note: This module handles preservation filtering only. Deposition
(whether an agent enters the archaeological record) is handled by a
separate module.

## See also

[`apply_deposition`](https://amy-s-anderson.github.io/persephone/reference/apply_deposition.md)
for the deposition filtering module that precedes preservation,
[`hsiler`](https://rdrr.io/pkg/demohaz/man/hsiler.html) for the Siler
hazard function.

## Examples

``` r
# Create a test population of deposited individuals
state <- data.frame(
  agent_id = 1:10,
  age = c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90),
  lesion = rep(FALSE, 10),
  dead = rep(TRUE, 10),
  was_deposited = rep(TRUE, 10),
  in_sample = rep(TRUE, 10)
)

# Siler parameters (Gage & Dyke 1986, Table 2, Level 15)
a_siler <- c(0.175, 1.40, 0.00368, 0.000075, 0.0917)
b_siler <- demohaz::trad_to_demohaz_siler_param(a_siler)

# Apply preservation filtering
result <- apply_preservation(
  state = state,
  preservation_model = "siler",
  preservation_param = b_siler,
  dx = 1
)
```
