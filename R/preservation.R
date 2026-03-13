#' @title Apply Preservation Filtering to Agent Population
#'
#' @description Applies preservation filtering to deposited agents. Determines
#'   which deposited skeletal remains are preserved and recoverable based on
#'   an age-dependent Siler hazard model.
#'
#' @details This function implements the second stage of archaeological
#'   filtration. Deposited agents are subject to a stochastic filtering
#'   process based on the Siler hazard function, which can model age-dependent
#'   preservation biases in the archaeological record.
#'
#'   The filtering uses a discrete-time approximation: for each deposited
#'   agent, the probability of being filtered out (not preserved) is
#'   \code{hsiler(age) * dx}. This approximation is valid when
#'   \code{hazard * dx << 1}.
#'
#'   Agents who are not deposited (\code{was_deposited = FALSE}) automatically
#'   have \code{in_sample} set to FALSE, as they cannot be recovered if they
#'   were never deposited.
#'
#'   Agents who are already \code{in_sample = FALSE} remain FALSE (no
#'   re-sampling).
#'
#'   Note: This module handles preservation filtering only. Deposition
#'   (whether an agent enters the archaeological record) is handled by a
#'   separate module.
#'
#' @param state A data.frame representing the current population state. Must
#'   contain the following columns:
#'   \describe{
#'     \item{agent_id}{Integer. Unique identifier for each agent.}
#'     \item{age}{Numeric. Age at death for dead agents.}
#'     \item{lesion}{Logical. TRUE if the agent has a lesion.}
#'     \item{dead}{Logical. TRUE if the agent is dead, FALSE if alive.}
#'     \item{was_deposited}{Logical. TRUE if the agent was deposited.}
#'     \item{in_sample}{Logical. Will be updated by this function.}
#'   }
#'
#' @param preservation_model Character. The preservation model to use.
#'   Currently only \code{"siler"} is supported.
#'
#' @param preservation_param Numeric vector containing the preservation
#'   parameters. For the \code{"siler"} model, this must be a length-5 vector
#'   containing the Siler hazard parameters in the demohaz parameterization:
#'   \code{c(b1, b2, b3, b4, b5)}. See \code{\link[demohaz]{hsiler}} for
#'   details.
#'
#' @param dx Numeric. The time step size in years. Determines the filtering
#'   probability (\code{hsiler(age) * dx}).
#'
#' @return A data.frame with the same structure as \code{state}, with updated
#'   values:
#'   \itemize{
#'     \item \code{in_sample}: Set to FALSE for non-deposited agents and for
#'       deposited agents who fail the Siler filter; unchanged for deposited
#'       agents who pass the filter.
#'     \item All other columns: Unchanged.
#'   }
#'
#' @examples
#' # Create a test population of deposited individuals
#' state <- data.frame(
#'   agent_id = 1:10,
#'   age = c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90),
#'   lesion = rep(FALSE, 10),
#'   dead = rep(TRUE, 10),
#'   was_deposited = rep(TRUE, 10),
#'   in_sample = rep(TRUE, 10)
#' )
#'
#' # Siler parameters (Gage & Dyke 1986, Table 2, Level 15)
#' a_siler <- c(0.175, 1.40, 0.00368, 0.000075, 0.0917)
#' b_siler <- demohaz::trad_to_demohaz_siler_param(a_siler)
#'
#' # Apply preservation filtering
#' result <- apply_preservation(
#'   state = state,
#'   preservation_model = "siler",
#'   preservation_param = b_siler,
#'   dx = 1
#' )
#'
#' @seealso \code{\link{apply_deposition}} for the deposition filtering
#'   module that precedes preservation, \code{\link[demohaz]{hsiler}} for
#'   the Siler hazard function.
#'
#' @export
apply_preservation <- function(state,
                               preservation_model,
                               preservation_param,
                               dx) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------

  if (is.null(preservation_model)) {
    stop("preservation_model must be specified")
  }

  # ---------------------------------------------------------------------------
  # Extract parameters (model-specific)
  # ---------------------------------------------------------------------------

  if (preservation_model == "siler") {
    if (length(preservation_param) != 5) {
      stop("For siler model, preservation_param must be length 5: c(b1, b2, b3, b4, b5)")
    }
    b_siler <- preservation_param

  } else {
    stop("preservation_model must be 'siler'. Other models are not yet supported.")
  }

  # ---------------------------------------------------------------------------
  # Copy state to avoid modifying input
  # ---------------------------------------------------------------------------

  result <- state

  # ---------------------------------------------------------------------------
  # Handle non-deposited agents
  # ---------------------------------------------------------------------------

  # Non-deposited agents cannot be in sample
  not_deposited_idx <- which(!state$was_deposited)
  result$in_sample[not_deposited_idx] <- FALSE

  # ---------------------------------------------------------------------------
  # Identify deposited agents eligible for preservation filtering
  # ---------------------------------------------------------------------------

  # Only deposited agents who are currently in_sample are filtered
  # (already FALSE stays FALSE)
  eligible_idx <- which(state$was_deposited & state$in_sample)

  # If no one is eligible, return early
  if (length(eligible_idx) == 0) {
    return(result)
  }

  # ---------------------------------------------------------------------------
  # Apply preservation filtering (model-specific)
  # ---------------------------------------------------------------------------

  if (preservation_model == "siler") {
    # Get ages of eligible agents
    ages <- state$age[eligible_idx]

    # Calculate Siler hazard at each age
    h <- demohaz::hsiler(ages, b_siler)

    # Probability of being filtered out (not preserved)
    p_filtered <- h * dx

    # Stochastic filtering
    u <- runif(length(eligible_idx))
    filtered_out <- u < p_filtered
    filtered_idx <- eligible_idx[filtered_out]

    # Mark filtered agents as not in sample
    result$in_sample[filtered_idx] <- FALSE
  }

  # ---------------------------------------------------------------------------
  # Return updated state
  # ---------------------------------------------------------------------------

  return(result)
}

