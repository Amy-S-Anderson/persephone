# Roll for lesion formation for a single agent

Roll for lesion formation for a single agent

## Usage

``` r
form_lesion(
  cohort,
  i,
  formation_window_opens,
  formation_window_closes,
  lesion_formation_rate
)
```

## Arguments

- cohort:

  Population data frame

- i:

  Index of the agent

- formation_window_opens:

  Age at which lesions can start forming

- formation_window_closes:

  Age at which lesions stop forming

- lesion_formation_rate:

  Annual probability of lesion formation

## Value

Updated cohort data frame
