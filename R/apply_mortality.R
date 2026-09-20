


#' @title Apply Mortality to Agent Population
#'
#' @description Applies a single time step of mortality to a population of
#'   agents, based on age-dependent mortality hazards (defined by Siler
#'   function parameters) and, optionally, other sources of mortality risk
#'   (e.g., individual frailty, skeletal lesion presence, etc.). Optionally
#'   splits mortality into competing "illness" and "trauma" causes, so that
#'   risk factors (e.g., lesions) can be specified to affect illness-related
#'   mortality without affecting trauma-related mortality.
#'
#'   This function is an orchestrator: it calls three sub-functions in
#'   sequence, each independently documented and testable:
#'   \enumerate{
#'     \item \code{compute_hazard_multiplier()} -- risk factors + exposure
#'       -> a single per-agent multiplier
#'     \item \code{compute_cause_specific_hazards()} -- age-based risk +
#'       multiplier + trauma/illness split -> illness_hazard, trauma_hazard
#'     \item \code{resolve_deaths()} -- hazards -> death dice, cause
#'       assignment, survivor/decedent split
#'   }
#'   See each sub-function's documentation for the logic it owns.
#'
#' @param pop Population data frame
#' @param current_time Integer, the current time step in the model
#' @param mu0_lookup data.frame, created by defrail_siler at the top of
#'   Simulate_Cemetery, optionally merged with a p_trauma column via
#'   build_p_trauma_lookup(). If mu0_lookup has no p_trauma column, or
#'   p_trauma is 0 at every age, all mortality is attributed to illness and
#'   no cause_of_death column is added to decedents (fully backward
#'   compatible behavior).
#' @param mortality_risk_type One of "proportional", "time_decreasing",
#'   "time_increasing"
#' @param force_death Logical. Defaults to FALSE. If TRUE, every agent dies
#'   this step regardless of hazard (used e.g. to close out a simulation
#'   when the population has dwindled to a handful of survivors); cause of
#'   death is still assigned from real age-specific hazards. See
#'   \code{resolve_deaths()}.
#' @param risk_factors Named list of risk factor columns. For continuous
#'   columns (e.g. frailty, acquired_frailty) whose values ARE the
#'   proportional hazard, set the value to NULL. For binary (0/1) columns
#'   (e.g. lesion), supply the proportional hazard as a numeric scalar
#'   (e.g. list(lesion = 2.1)). NOTE: when cause-of-death splitting is
#'   active (i.e. mu0_lookup$p_trauma has nonzero values), all risk_factors
#'   are currently applied to illness-related mortality only; trauma-related
#'   mortality reflects the age-based baseline hazard alone. This matches
#'   the intended use case (e.g. skeletal lesions raising risk of dying from
#'   illness/infection but not from accidents), but means risk_factors
#'   intended to modify trauma risk are not yet supported.
#' @param exposure_hazard_multiplier Numeric. Multiplier applied to exposed
#'   agents when the exposure hazard is transient.
#'
#' @return A list with:
#'   \item{pop}{Surviving agents.}
#'   \item{decedents}{Agents who died this step, with year_died (and
#'     cause_of_death, if cause-of-death splitting is active) added.}
#' @keywords internal
#' @export
apply_mortality <- function(pop,
                            mu0_lookup,
                            mortality_risk_type         = "proportional",
                            force_death                 = FALSE,
                            risk_factors                = list(),
                            current_time,
                            exposure_hazard_multiplier  = 1) { # If this argument isn't specified in the ABM loop, then it defaults to 1 (no multiplying effect of exposure)
  
  # updated acquired_frailty, if relevant:
  if ("acquired_frailty" %in% names(pop)) {
    pop <- acquire_frailty(pop, exposure_hazard_multiplier)
  }
  
  # --- Step 1: risk factors + exposure -> hazard multiplier --------------
  hazard_multiplier <- compute_hazard_multiplier(
    pop                        = pop,
    risk_factors               = risk_factors
  )

  # --- Step 2: age-based risk + multiplier + trauma/illness split --------
  cause_hazards <- compute_cause_specific_hazards(
    pop                  = pop,
    mu0_lookup           = mu0_lookup,
    hazard_multiplier    = hazard_multiplier,
    mortality_risk_type  = mortality_risk_type
  )

  # --- Step 3: death dice, cause assignment, survivor/decedent split -----
  result <- resolve_deaths(
    pop                      = pop,
    illness_hazard           = cause_hazards$illness_hazard,
    trauma_hazard            = cause_hazards$trauma_hazard,
    split_active             = cause_hazards$split_active,
    current_time             = current_time,
    force_death              = force_death
  )

  result
}


