# Roll for death for a single agent (v1 implementation)

Roll for death for a single agent (v1 implementation)

## Usage

``` r
apply_mortality_v1(
  cohort,
  i,
  age_based_risk,
  mortality_risk_type,
  relative_mortality_risk
)
```

## Arguments

- cohort:

  Population data frame

- i:

  Index of the agent

- age_based_risk:

  Baseline Siler mortality risk at current age

- mortality_risk_type:

  One of "proportional", "time_decreasing", "time_increasing"

- relative_mortality_risk:

  Multiplier for lesion-bearing individuals

## Value

Updated cohort data frame
