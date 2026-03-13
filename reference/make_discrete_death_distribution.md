# Compute discrete age-at-death distribution from a probability function

Compute discrete age-at-death distribution from a probability function

## Usage

``` r
make_discrete_death_distribution(prob_fun, age_max = 100)
```

## Arguments

- prob_fun:

  Function taking age and returning death probability.

- age_max:

  Integer. Maximum age to compute. Default 100.

## Value

Data frame with columns age, qx, px, S_start, dx.
