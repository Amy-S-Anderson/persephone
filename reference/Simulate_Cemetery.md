# Simulate a cemetery using the Persephone ABM

Runs the Persephone agent-based model: a birth cohort ages through time,
facing annual risks of skeletal lesion formation and Siler-model
mortality. Individuals with lesions may experience modified mortality
risk. The simulation ends when fewer than 10 individuals remain alive.

## Usage

``` r
Simulate_Cemetery(
  cohort_size,
  lesion_formation_rate,
  formation_window_opens = 0,
  formation_window_closes,
  mortality_risk_type = "proportional",
  relative_mortality_risk = 1,
  mortality_regime,
  deposition_param = 0,
  loss_strength = "none",
  age_noise = FALSE
)
```

## Arguments

- cohort_size:

  Integer. Number of individuals in the starting cohort.

- lesion_formation_rate:

  Numeric. Annual probability of developing a lesion (between 0 and 1).

- formation_window_opens:

  Numeric. Age at which lesions can start forming. Default 0.

- formation_window_closes:

  Numeric. Age at which new lesions stop forming.

- mortality_risk_type:

  Character. How lesions modify mortality: "proportional",
  "time_decreasing", or "time_increasing".

- relative_mortality_risk:

  Numeric. Mortality multiplier for individuals with lesions. 1 = no
  effect, 2 = double risk. Default 1.

- mortality_regime:

  Data frame with Siler parameters (a1, b1, a2, a3, b3).

## Value

A list with two elements:

- individual_outcomes:

  Data frame of all individuals with age at death and lesion status.

- survivors:

  Data frame of survivor counts and lesion prevalence at each age.

## Examples

``` r
result <- Simulate_Cemetery(
  cohort_size = 500,
  lesion_formation_rate = 0.10,
  formation_window_closes = 5,
  mortality_regime = CoaleDemenyWestF5
)
```
