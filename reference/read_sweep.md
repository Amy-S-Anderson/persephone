# Read and combine saved sweep results

Read and combine saved sweep results

## Usage

``` r
read_sweep(mortality_regime, lesion_rates, reps, rmr)
```

## Arguments

- mortality_regime:

  Data frame with Siler parameters

- lesion_rates:

  Vector of lesion formation rates that were swept

- reps:

  Number of replicates (unused but kept for API consistency)

- rmr:

  Relative mortality risk value

## Value

Combined data frame of all sweep results
