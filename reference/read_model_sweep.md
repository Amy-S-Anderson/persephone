# Read saved sweep results from disk

Read saved sweep results from disk

## Usage

``` r
read_model_sweep(
  root_output_directory,
  target_param,
  target_param_values,
  data_file = "sim_cemetery.csv"
)
```

## Arguments

- root_output_directory:

  Root directory containing sweep outputs

- target_param:

  Character name of the swept parameter

- target_param_values:

  Vector of parameter values to read

- data_file:

  Filename to read from each replicate directory

## Value

Nested list of data frames
