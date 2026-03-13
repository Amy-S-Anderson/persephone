# Create Empirical Distribution from Assemblage Data

Creates an empirical probability distribution from aggregated assemblage
age-at-death data. The distribution uses piecewise constant density
within finite age intervals and an exponential tail for the final
unbounded interval (if present).

## Usage

``` r
create_assemblage_distribution(data)
```

## Arguments

- data:

  A data.frame containing assemblage age-at-death data with the
  following required columns:

  study

  :   Character. Identifier for the assemblage/study.

  x_lo

  :   Numeric. Lower bound of the age interval. Can be negative (e.g.,
      -0.5 for pre-birth deaths).

  x_hi

  :   Numeric. Upper bound of the age interval. Can be `Inf` for the
      final open-ended age class.

  count

  :   Numeric. Number of individuals in this interval for this study.
      Must be non-negative.

## Value

A list containing:

- dfun:

  Function. The density function. Takes a numeric vector of ages and
  returns corresponding density values.

- rfun:

  Function. The sampling function. Takes an integer n and returns n
  random samples from the distribution.

- cuts:

  Numeric vector. The interval boundaries (cut points) in ascending
  order. Does not include Inf.

- densities:

  Numeric vector. The piecewise constant density values for each finite
  interval. Length is `length(cuts) - 1`.

- exp_rate:

  Numeric. The rate parameter (lambda) for the exponential tail
  distribution. `NA` if no unbounded interval exists.

- tail_density:

  Numeric. The density at the start of the exponential tail (for
  continuity). `NA` if no unbounded interval exists.

## Details

This function aggregates interval-censored age-at-death data across
multiple archaeological studies/assemblages and constructs a normalized
probability distribution. The approach is:

1.  Aggregate counts across studies for overlapping intervals

2.  Identify unique cut points from interval boundaries

3.  Assign piecewise constant density within each finite interval,
    proportional to the count and inversely proportional to width

4.  For unbounded final intervals (x_hi = Inf), use an exponential tail
    distribution with rate derived from expected remaining life

5.  Normalize so the total probability integrates to 1

**Note on the exponential tail:** The rate parameter for the exponential
tail is derived from the expected remaining life, which is a notional
approximation. This assumes that individuals in the final age class have
an average remaining lifespan equal to the midpoint of the preceding
finite interval's width. This is a simplification and should be
interpreted cautiously.

**Limitations:** This approach does not account for uncertainty in age
estimation. It treats interval-censored data as if uniformly distributed
within each interval. This is a notional/exploratory tool.

## Examples

``` r
# Simple assemblage with two intervals
data <- data.frame(
  study = c("site_A", "site_A"),
  x_lo = c(0, 20),
  x_hi = c(20, Inf),
  count = c(30, 70)
)
dist <- create_assemblage_distribution(data)

# Evaluate density at various ages
dist$dfun(c(10, 25, 50))
#> [1] 0.015000000 0.027258027 0.007809556

# Sample 100 ages from the distribution
samples <- dist$rfun(100)
```
