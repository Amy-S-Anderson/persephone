#' @title Create Empirical Distribution from Assemblage Data
#'
#' @description Creates an empirical probability distribution from aggregated
#'   assemblage age-at-death data. The distribution uses piecewise constant
#'   density within finite age intervals and an exponential tail for the
#'   final unbounded interval (if present).
#'
#' @details This function aggregates interval-censored age-at-death data across
#'   multiple archaeological studies/assemblages and constructs a normalized
#'   probability distribution. The approach is:
#'
#'   \enumerate{
#'     \item Aggregate counts across studies for overlapping intervals
#'     \item Identify unique cut points from interval boundaries
#'     \item Assign piecewise constant density within each finite interval,
#'           proportional to the count and inversely proportional to width
#'     \item For unbounded final intervals (x_hi = Inf), use an exponential
#'           tail distribution with rate derived from expected remaining life
#'     \item Normalize so the total probability integrates to 1
#'   }
#'
#'   \strong{Note on the exponential tail:} The rate parameter for the
#'   exponential tail is derived from the expected remaining life, which is
#'   a notional approximation. This assumes that individuals in the final
#'   age class have an average remaining lifespan equal to the midpoint of
#'   the preceding finite interval's width. This is a simplification and
#'   should be interpreted cautiously.
#'
#'   \strong{Limitations:} This approach does not account for uncertainty in
#'   age estimation. It treats interval-censored data as if uniformly
#'   distributed within each interval. This is a notional/exploratory tool.
#'
#' @param data A data.frame containing assemblage age-at-death data with the
#'   following required columns:
#'   \describe{
#'     \item{study}{Character. Identifier for the assemblage/study.}
#'     \item{x_lo}{Numeric. Lower bound of the age interval. Can be negative
#'       (e.g., -0.5 for pre-birth deaths).}
#'     \item{x_hi}{Numeric. Upper bound of the age interval. Can be \code{Inf}
#'       for the final open-ended age class.}
#'     \item{count}{Numeric. Number of individuals in this interval for this
#'       study. Must be non-negative.}
#'   }
#'
#' @return A list containing:
#'   \describe{
#'     \item{dfun}{Function. The density function. Takes a numeric vector of
#'       ages and returns corresponding density values.}
#'     \item{rfun}{Function. The sampling function. Takes an integer n and
#'       returns n random samples from the distribution.}
#'     \item{cuts}{Numeric vector. The interval boundaries (cut points) in
#'       ascending order. Does not include Inf.}
#'     \item{densities}{Numeric vector. The piecewise constant density values
#'       for each finite interval. Length is \code{length(cuts) - 1}.}
#'     \item{exp_rate}{Numeric. The rate parameter (lambda) for the exponential
#'       tail distribution. \code{NA} if no unbounded interval exists.}
#'     \item{tail_density}{Numeric. The density at the start of the exponential
#'       tail (for continuity). \code{NA} if no unbounded interval exists.}
#'   }
#'
#' @examples
#' # Simple assemblage with two intervals
#' data <- data.frame(
#'   study = c("site_A", "site_A"),
#'   x_lo = c(0, 20),
#'   x_hi = c(20, Inf),
#'   count = c(30, 70)
#' )
#' dist <- create_assemblage_distribution(data)
#'
#' # Evaluate density at various ages
#' dist$dfun(c(10, 25, 50))
#'
#' # Sample 100 ages from the distribution
#' samples <- dist$rfun(100)
#'
#' @export
create_assemblage_distribution <- function(data) {
  
  # ---------------------------------------------------------------------------
 # Input validation
  # ---------------------------------------------------------------------------
  
  # Check that data is a data.frame
 if (!is.data.frame(data)) {
    stop("data must be a data.frame")
  }
  
  # Check for required columns
  required_cols <- c("study", "x_lo", "x_hi", "count")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Check that x_lo is numeric
  if (!is.numeric(data$x_lo)) {
    stop("x_lo must be numeric")
  }
  
  # Check that x_hi is numeric
  if (!is.numeric(data$x_hi)) {
    stop("x_hi must be numeric")
  }
  
  # Check that count is numeric
  if (!is.numeric(data$count)) {
    stop("count must be numeric")
  }
  
  # Check for negative counts
  if (any(data$count < 0)) {
    stop("count cannot contain negative values")
  }
  
  # Check that x_lo < x_hi for all rows
  if (any(data$x_lo >= data$x_hi)) {
    stop("x_lo must be less than x_hi for all intervals")
  }
  
  # Check that total count is non-zero
  if (sum(data$count) == 0) {
    stop("Total count cannot be zero")
  }
  
  # ---------------------------------------------------------------------------
  # Aggregate counts and identify intervals
  # ---------------------------------------------------------------------------
  
  # Determine if there's an unbounded (Inf) interval
  has_unbounded <- any(is.infinite(data$x_hi))
  
  # Get all unique finite cut points
  all_boundaries <- c(data$x_lo, data$x_hi)
  finite_boundaries <- all_boundaries[is.finite(all_boundaries)]
  cuts <- sort(unique(finite_boundaries))
  
  # Number of finite intervals
  n_intervals <- length(cuts) - 1
  
  # Aggregate counts for each interval defined by consecutive cuts.
  # Each row's count is distributed proportionally across the cut intervals
  # it spans, based on relative widths (assuming uniform distribution within
  # the original interval).
  #
  # Special case: rows with x_hi = Inf have all their count assigned to the
  # unbounded tail (we cannot assume uniform distribution over infinity).
  interval_counts <- numeric(n_intervals)
  unbounded_count <- 0
  
  for (j in seq_len(nrow(data))) {
    row_lo <- data$x_lo[j]
    row_hi <- data$x_hi[j]
    row_count <- data$count[j]
    
    # If x_hi is Inf, all count goes to the unbounded tail
    if (is.infinite(row_hi)) {
      unbounded_count <- unbounded_count + row_count
      next
    }
    
    # For finite intervals, distribute count proportionally across cut intervals
    row_width <- row_hi - row_lo
    
    for (i in seq_len(n_intervals)) {
      cut_lo <- cuts[i]
      cut_hi <- cuts[i + 1]
      
      # Calculate overlap between [row_lo, row_hi) and [cut_lo, cut_hi)
      overlap_lo <- max(row_lo, cut_lo)
      overlap_hi <- min(row_hi, cut_hi)
      overlap_width <- max(0, overlap_hi - overlap_lo)
      
      if (overlap_width > 0) {
        # Proportion of row's count that falls in this cut interval
        proportion <- overlap_width / row_width
        interval_counts[i] <- interval_counts[i] + row_count * proportion
      }
    }
  }
  
  # Total count for normalization
  total_count <- sum(interval_counts) + unbounded_count
  
  # ---------------------------------------------------------------------------
  # Calculate densities
  # ---------------------------------------------------------------------------
  
  # Proportions for each interval
  interval_props <- interval_counts / total_count
  unbounded_prop <- unbounded_count / total_count
  
  # Interval widths
  widths <- diff(cuts)
  
  # Piecewise constant densities (proportion / width)
  densities <- interval_props / widths
  
  # ---------------------------------------------------------------------------
  # Exponential tail (if unbounded interval exists)
  # ---------------------------------------------------------------------------
  
  if (has_unbounded && unbounded_prop > 0) {
    # Derive exponential rate from expected remaining life
    # This is a notional approximation: we use the average width of finite
    # intervals as a rough estimate of expected remaining life
    # 
    # For an exponential distribution: E[X - x0 | X > x0] = 1/lambda
    # So lambda = 1 / expected_remaining_life
    #
    # We estimate expected_remaining_life as the mean width of finite intervals
    # This is crude but serves our notional purposes
    expected_remaining_life <- mean(widths)
    exp_rate <- 1 / expected_remaining_life
    
    # The tail density at the start must be chosen so that:
    # (sum of finite densities * widths) + (tail_density / exp_rate) = 1
    # We already know unbounded_prop should integrate to unbounded_prop
    # So: tail_density / exp_rate = unbounded_prop
    # Therefore: tail_density = unbounded_prop * exp_rate
    tail_density <- unbounded_prop * exp_rate
  } else {
    exp_rate <- NA
    tail_density <- NA
  }
  
  # ---------------------------------------------------------------------------
  # Define density function (dfun)
  # ---------------------------------------------------------------------------
  
  dfun <- function(x) {
    # Vectorized density evaluation
    result <- numeric(length(x))
    
    for (i in seq_along(x)) {
      xi <- x[i]
      
      # Check if outside the defined range
      if (xi < cuts[1]) {
        result[i] <- 0
        next
      }
      
      # Check finite intervals
      in_finite <- FALSE
      for (j in seq_len(n_intervals)) {
        if (xi >= cuts[j] && xi < cuts[j + 1]) {
          result[i] <- densities[j]
          in_finite <- TRUE
          break
        }
      }
      
      if (in_finite) next
      
      # Check exponential tail
      if (has_unbounded && !is.na(exp_rate) && xi >= cuts[length(cuts)]) {
        # Exponential decay from the last cut point
        x0 <- cuts[length(cuts)]
        result[i] <- tail_density * exp(-exp_rate * (xi - x0))
      } else {
        # Beyond range with no tail
        result[i] <- 0
      }
    }
    
    return(result)
  }
  
  # ---------------------------------------------------------------------------
  # Define sampling function (rfun)
  # ---------------------------------------------------------------------------
  
  rfun <- function(n) {
    # Sample from the piecewise constant + exponential tail distribution
    samples <- numeric(n)
    
    # Probability weights for each interval (including tail)
    if (has_unbounded && unbounded_prop > 0) {
      probs <- c(interval_props, unbounded_prop)
      n_regions <- n_intervals + 1
    } else {
      probs <- interval_props
      n_regions <- n_intervals
    }
    
    # Sample which interval each draw comes from
    interval_idx <- sample(seq_len(n_regions), n, replace = TRUE, prob = probs)
    
    for (i in seq_len(n)) {
      idx <- interval_idx[i]
      
      if (idx <= n_intervals) {
        # Sample uniformly within finite interval
        samples[i] <- runif(1, min = cuts[idx], max = cuts[idx + 1])
      } else {
        # Sample from exponential tail
        x0 <- cuts[length(cuts)]
        samples[i] <- x0 + rexp(1, rate = exp_rate)
      }
    }
    
    return(samples)
  }
  
  # ---------------------------------------------------------------------------
  # Return result
  # ---------------------------------------------------------------------------
  
  list(
    dfun = dfun,
    rfun = rfun,
    cuts = cuts,
    densities = densities,
    exp_rate = exp_rate,
    tail_density = tail_density
  )
}

