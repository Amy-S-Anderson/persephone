# Sweep across values of a target parameter

Sweep across values of a target parameter

## Usage

``` r
run_model_sweep(
  base_params,
  root_output_directory,
  target_param,
  target_param_values,
  numreps = 10
)
```

## Arguments

- base_params:

  Named list of base simulation parameters

- root_output_directory:

  Root directory for saving sweep outputs

- target_param:

  Character name of the parameter to sweep

- target_param_values:

  Vector of values to sweep over

- numreps:

  Number of replicates per parameter value

## Value

Combined data frame of all simulation results
