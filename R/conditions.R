# =====================================================================
# conditions.R
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   Builds the package's explicit error and warning objects. These carry
#   a class as well as a message, so a test can check that one specific
#   rule fired instead of matching English text.
#
# PLACE IN THE PIPELINE
#   Used by every layer: registry.R, solver.R, the model validator, and
#   the engine.
#
# WHAT YOU NEED TO KNOW
#   An R error is an object, not just a string. When an error carries the
#   class c("my_error", "error", "condition"), calling code can catch that
#   exact error with tryCatch(..., my_error = function(e) ...), and a test
#   can write expect_error(..., class = "my_error").
#
# GUARANTEES
#   * Errors built here have the supplied micsim_ subclass. Errors from
#     base R or user-written callbacks may retain their original classes.
#   * Errors are raised with stop(), never printed with cat(), so they can
#     be caught, suppressed and tested.
# =====================================================================


#' Build a MicSim error object with a specific class
#'
#' @param subclass Character vector of class names, most specific first,
#'   e.g. `"micsim_rate_error"`. `"error"` and `"condition"` are appended
#'   automatically.
#' @param message The message shown to the user. One or more lines.
#' @param call The call to report, or `NULL` to suppress it. `NULL` is the
#'   default because MicSim's errors are about the *model*, not about the
#'   internal function that happened to notice the problem, and printing
#'   the internal call misleads the reader.
#' @param ... Further named fields stored on the condition object. Tests
#'   and calling code may read them; for example the model validator
#'   attaches `problems`, the full character vector of faults.
#'
#' @return A condition object. This function does not raise it; use
#'   `micsim_abort()` for that.
#' @keywords internal
#' @noRd
micsim_condition <- function(subclass, message, call = NULL, ...) {
  structure(
    class = c(subclass, "error", "condition"),
    list(message = message, call = call, ...)
  )
}

#' Raise a classed MicSim error
#'
#' Calls `stop()` with an error object from `micsim_condition()`.
#' Internal checks use this helper to report a problem with its MicSim
#' error class and message.
#'
#' @inheritParams micsim_condition
#' @return Never returns.
#' @keywords internal
#' @noRd

micsim_abort <- function(subclass, message, call = NULL, ...) {
  stop(micsim_condition(subclass, message, call = call, ...))
}


#' Build a classed MicSim warning
#'
#' The warning counterpart of `micsim_condition()`, so that a caller can
#' muffle one specific warning with `withCallingHandlers()` and a test can
#' write `expect_warning(..., class = )`.
#'
#' @param subclass,message,call,... As in `micsim_condition()`.
#' @return A condition object of class `c(subclass, "warning", "condition")`.
#' @keywords internal
#' @noRd
micsim_warning <- function(subclass, message, call = NULL, ...) {
  structure(
    class = c(subclass, "warning", "condition"),
    list(message = message, call = call, ...)
  )
}

#' Raise the "this rate function is unusable" error
#'
#' Used by rate checks when a supplied function has unsupported argument
#' names, does not return one value per input, returns `NA`, or returns
#' a negative rate. The calling check supplies the explanation; this
#' helper raises the corresponding error.
#'
#' @param rate_label Human-readable transition label, e.g. `"f/DE/0 -> f/DE/1"`.
#' @param problem One sentence describing what is wrong.
#' @param ... Further fields stored on the condition (e.g. `argnames`).
#' @return Never returns.
#' @keywords internal
#' @noRd
abort_rate <- function(rate_label, problem, ...) {
  micsim_abort(
    "micsim_rate_error",
    paste0("Transition rate for `", rate_label, "` is unusable.\n  ", problem),
    ...
  )
}

#' Raise the "these arguments do not describe a model" error
#'
#' Used by the model validator, which collects *all* faults and reports
#' them in one go. `problems` is kept as a field so that a test can assert
#' on an individual rule without parsing the assembled message.
#'
#' @param problems A data frame with one row per finding,
#' as built in R/validate.R, or a character vector with one fault per element.
#' @return Never returns.
#' @keywords internal
#' @noRd
abort_invalid_model <- function(problems) {
  # A character vector is one finding per element; validate() passes a
  # data frame with one row per finding (see R/validate.R).
  if (is.character(problems)) {
    problems <- data.frame(transition = NA_character_, check = NA_character_,
                           severity = "error", where = NA_character_,
                           n_points = NA_integer_, detail = problems,
                           stringsAsFactors = FALSE)
  }
  # The condition also carries the classes of the kinds of fault found, so
  # a caller catching e.g. `micsim_rate_error` still catches a model whose
  # rates are faulty.
  kinds <- c(
    if (any(problems$check %in% c("evaluation", "type", "length", "NA",
                                  "negative"))) "micsim_rate_error",
    if (any(problems$check %in% "newborn")) "micsim_newborn_rule_error"
  )
  micsim_abort(
    c("micsim_invalid_model", kinds),
    format_problems(problems, "problem(s) found in the simulation model"),
    problems = problems
  )
}

#' Raise the "the numerical solve failed" error
#'
#' Used when the numerical search for a waiting time fails. For example,
#' the search can fail to find a time at which accumulated risk reaches
#' the drawn value, even though the endpoint calculations indicate one
#' should exist. The simulation stops rather than returning that failed
#' calculation as a waiting time.
#'
#' @param problem One sentence describing what failed.
#' @param ... Further fields stored on the condition.
#' @return Never returns.
#' @keywords internal
#' @noRd
abort_solver <- function(problem, ...) {
  micsim_abort("micsim_solver_error", problem, ...)
}
