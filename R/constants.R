# =====================================================================
# constants.R
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   Defines shared values for time conversion, rate-function arguments,
#   numerical calculations and state labels. These are internal settings,
#   not arguments users normally supply.
#
# PLACE IN THE PIPELINE
#   Read from everywhere; depends on nothing. It lets the rest of the
#   package be read without guessing what a literal such as 365.25 or 1e-7
#   means.
# =====================================================================


# ---------------------------------------------------------------------
# Time
# ---------------------------------------------------------------------

#' Days per year, for the day-count-to-years conversion
#'
#' Rate functions are defined on fractional years, so day counts must be
#' converted to years. This factor does that conversion: it builds the
#' `age`, `calTime` and `duration` passed to a rate function, and the age
#' reported in the output (`age_days / DAYS_PER_YEAR`, rounded to two
#' decimals). Calendar dates are never scaled by it; they use the exact
#' calendar in [getInDays()].
#'
#' The value is the mean length of a Gregorian year.
#'
#' @keywords internal
#' @noRd
DAYS_PER_YEAR <- 365.25


#' The margin by which an individual stops being at risk before `max_age`
#'
#' Ordinary waiting-time calculations end at
#' `max_age - MAX_AGE_MARGIN`, about 3.65 days before the maximum age.
#'
#' This matters for rate tables indexed by whole age using
#' `as.integer(age)`. When `max_age` is a whole number, a table covering
#' ages 0 to `max_age - 1` therefore needs no entry for age `max_age`
#' during these calculations.
#'
#' Episode tables end observation at the same boundary; see
#' [convertToLongFormat()].
#'
#' @keywords internal
#' @noRd
MAX_AGE_MARGIN <- 0.01


# ---------------------------------------------------------------------
# Rate functions
# ---------------------------------------------------------------------

#' The only argument names a transition-rate function may declare
#'
#' Rates may use any subset of these arguments, in any order, including
#' no arguments for a constant rate. Other argument names are rejected.
#' The declared names tell the registry which time scales a rate uses; see
#' `new_rate_registry()`.
#'
#' @keywords internal
#' @noRd
RATE_ARGS <- c("age", "calTime", "duration")


# ---------------------------------------------------------------------
# Numerical tolerances
# ---------------------------------------------------------------------

#' Relative tolerance when retrying numerical integration
#'
#' After an integration error, the solver retries with a looser
#' tolerance, allowing a less precise calculation. It returns
#' `reliable = FALSE` to record the retry.
#'
#' @keywords internal
#' @noRd
INTEGRATE_FALLBACK_REL_TOL <- 1e-2

#' Tolerance for the waiting-time search
#'
#' Controls how closely [stats::uniroot()] locates the waiting time
#' at which accumulated hazard reaches the random threshold. The
#' tolerance is in years and is much smaller than one day. 
#'
#' @keywords internal
#' @noRd
UNIROOT_TOL <- .Machine$double.eps^0.5

#' Tolerance for checks that a set of probabilities sums to one
#'
#' Allows small numerical rounding differences when checking that
#' probabilities in newborn rules sum to one. Passed to
#' [base::all.equal()].
#'
#' @keywords internal
#' @noRd
PROB_TOL <- 1e-7


# ---------------------------------------------------------------------
# State labels
# ---------------------------------------------------------------------

#' Separator between state-variable values
#'
#' Joins values into labels such as `"f/DE/1"`. Individual values must
#' not contain `"/"`, because it separates the label into its parts.
#'
#' The model-object interface specifies the meaning of states through
#' explicit declarations and rules. The classic interface retains
#' conventions for `"dead"`, `"rest"`, `"f"`, `"no"` and `"low"`;
#' see `legacy_roles()` in `R/compat.R`.
#'
#' @keywords internal
#' @noRd
STATE_SEP <- "/"

#' Allowed descriptions of exits into absorbing states
#'
#' Used by the `exits` argument of [micsim_model()] to identify death,
#' emigration or another exit. An unspecified kind is recorded as `NA`.
#' The core engine records these descriptions but handles all absorbing
#' states alike. A subsystem whose response depends on why someone left
#' the population can look up the exit kind in the model.
#'
#' @keywords internal
#' @noRd
EXIT_KINDS <- c("death", "emigration", "other")
