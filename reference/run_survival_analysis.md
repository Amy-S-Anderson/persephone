# Run survival analysis on sweep results

Fits Kaplan-Meier curves and log-rank tests comparing survival between
individuals with and without lesions, for each combination of mortality
regime, lesion formation rate, and replicate.

## Usage

``` r
run_survival_analysis(sweep_data, parallel = TRUE, workers = NULL)
```

## Arguments

- sweep_data:

  Data frame of combined sweep results

- parallel:

  Logical. Use parallel processing? Default TRUE.

- workers:

  Number of parallel workers. Default: detectCores() - 1.

## Value

A list with survival_data (for plotting) and logrank_results
