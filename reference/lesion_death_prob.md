# Compute lesion-modified Siler death probability

Compute lesion-modified Siler death probability

## Usage

``` r
lesion_death_prob(age, regime, risk_type, rmr)
```

## Arguments

- age:

  Numeric. Age in years.

- regime:

  Data frame with Siler parameters (a1, b1, a2, a3, b3).

- risk_type:

  Character. One of "proportional", "time_decreasing",
  "time_increasing".

- rmr:

  Numeric. Relative mortality risk multiplier.

## Value

Numeric death probability.
