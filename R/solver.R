# =====================================================================
# solver.R
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   Draws how long an individual waits in their current state before one
#   specific transition happens, ignoring every other transition.
#
# PLACE IN THE PIPELINE
#   registry.R --(hazard closures)--> solver.R --(waiting time)--> engine.R
#
#   The engine calls this file once per possible destination and keeps the
#   earliest waiting time, which is how competing risks are resolved. The
#   solver itself knows nothing about competition.
#
# WHAT YOU NEED TO KNOW  (concepts in ?MicSim-concepts)
#   hazard rate, integrated hazard, inverse-transform sampling, forced
#   transition.
#
# GUARANTEES
#   * The solver reads and writes nothing outside its own arguments. The
#     engine hands it the random budget (`budget`); called without one, it
#     draws exactly one through `draw_u`, never by calling runif() directly.
#   * It returns Inf whenever no event occurs inside the remaining
#     horizon. Inf is not an error here, but a normal answer.
#   * It never hides a degraded numerical integration. Every return value
#     carries `reliable`, and the engine counts the FALSE ones.
#
# UNITS
#   Everything in this file is in fractional YEARS. There are no day
#   counts and no dates here. The engine does the converting
#   on the way in and out.
#
# NOTE ON RANDOM NUMBERS
#   Every call uses exactly one budget, whichever path it takes, so the
#   draws of other destinations and later schedulings never move. The
#   engine draws the budgets of a whole scheduling at once from the
#   individual's own stream (rng_streams.R) and passes each in `budget`.
# =====================================================================


#' Draw a standard-exponential variate
#'
#' The default random-number source for `draw_waiting_time()`. The
#' expression is `-log(1 - runif(1))` rather than the equivalent-in-law
#' `-log(runif(1))`. Both give the same distribution, but use the uniform draw
#' differently.
#'
#' @return A single positive numeric value.
#' @keywords internal
#' @noRd
default_draw_u <- function() {
  -log(1 - stats::runif(1L))
}


#' Integrate a hazard, reporting degraded accuracy
#'
#' @description
#' Computes the integrated hazard \eqn{\int_{lower}^{upper} \lambda(s)\,
#' \mathrm{d}s} by adaptive quadrature, and says so when it had to
#' compromise.
#'
#' @details
#' [stats::integrate()] fails on integrands it cannot resolve, typically
#' narrow spikes. For example, a fertility schedule that is zero everywhere 
#' except a 30-year window can defeat the default tolerance when the integration
#' range is a whole lifetime. Retrying at a looser tolerance is better than
#' not integrating at all, but a looser tolerance means a less accurate answer 
#' in the main numerical step of the whole method.
#' So the retry still happens, but the fact that it happened is returned to the 
#' caller and ends up in the run diagnostics.
#'
#' @param lambda Vectorised hazard, a function of years elapsed.
#' @param lower,upper Integration limits, in years.
#'
#' @return A list with
#'   \describe{
#'     \item{`value`}{the integrated hazard, a non-negative number}
#'     \item{`reliable`}{`FALSE` if the default tolerance failed and the
#'       looser tolerance was used}
#'   }
#' @seealso [MicSim-concepts]
#' @keywords internal
#' @noRd
integrate_hazard <- function(lambda, lower, upper, label = "this transition") {
  attempt <- tryCatch(
    stats::integrate(lambda, lower, upper),
    # Calling the rate can fail. It's error already names the transition, so
    # it is passed on. Retrying at a looser tolerance wouldn't fix this.
    micsim_rate_error = function(e) stop(e),
    error = function(e) {
      # An Inf the forced finder did not see makes integrate() fail. Say
      # so in terms of the rate (checked by value, not by the wording of
      # integrate()'s message, which is translated).
      probe <- tryCatch(lambda(seq(lower, upper, length.out = 1001L)),
                        error = function(e2) NULL)
      if (is.numeric(probe) && any(is.infinite(probe))) {
        abort_inf_between(label)
      }
      NULL
    }
  )
  if (!is.null(attempt)) {
    return(list(value = attempt$value, reliable = TRUE))
  }

  # Fallback. stop.on.error = FALSE makes integrate() return its best
  # estimate instead of failing. The reliability flag is a RETURN VALUE,
  # not a superassignment into the enclosing frame, we don't use
  # `<<-`.
  degraded <- stats::integrate(
    lambda, lower, upper,
    rel.tol       = INTEGRATE_FALLBACK_REL_TOL,
    stop.on.error = FALSE
  )
  list(value = degraded$value, reliable = FALSE)
}


#' Draw a waiting time by inverting the integrated hazard
#'
#' @description
#' Draws how long an individual waits in their current state before
#' one specific transition occurs, as if that transition were the only one
#' available. The engine calls this once per possible destination and
#' keeps the earliest draw; see [MicSim-concepts], section
#' "Competing risks".
#'
#' @details
#' Write \eqn{\lambda(x)} for the transition rate \eqn{x} years from now,
#' and \eqn{\Lambda(\tau) = \int_0^{\tau} \lambda(s)\,\mathrm{d}s} for the
#' integrated hazard. If \eqn{T} denotes the waiting time, then
#' \eqn{\Lambda(T)} follows a standard exponential distribution. 
#' Then, drawing \eqn{u} from that distribution and solving
#' \deqn{\int_0^{\tau} \lambda(s)\,\mathrm{d}s = u}
#' for \eqn{\tau} yields a draw of \eqn{T}.
#'
#' In words: we accumulate risk forward in time until the accumulated risk
#' reaches a randomly drawn budget. The date at which the budget runs out
#' is the event date. If the budget is never exhausted before the
#' observation window closes, no event happens and the answer is `Inf`.
#'
#' Three things can shorten the procedure:
#'
#' 1. A rate may be `Inf`, which encodes an event on a *known* date rather
#'    than a random one (school enrolment). That is checked first, by the
#'    `forced_time()` finder `bind_hazard()` builds. The budget,
#'    which is still drawn, is discarded.
#' 2. A rate may be identically zero over the whole window (fertility for
#'    a 70-year-old), in which case no event can occur.
#' 3. The accumulated risk over the whole window may fall short of the
#'    drawn budget, again meaning no event.
#'
#' @section Assumptions:
#' * `lambda` is vectorised, non-negative and finite. Vectorisation is a
#'   hard requirement of the method and not a convenience, because the 
#'   quadrature routine evaluates it at many values of `x` in one call.
#' * The all-zero shortcut (step 2) tests the rate at whole years from
#'   now, at the end of the window, and just after every point at which
#'   age, calendar time or duration reaches a whole number (`crossings`).
#'   That is exact for any rate that is piecewise constant on whole ages,
#'   years and durations -- every rate table -- because the path through
#'   the window enters each cell of such a table at one of those points.
#'   The end of a window has to be included. Otherwise, a fertility table
#'   keyed on whole ages reads as zero for a woman whose window ends a 
#'   few months after her 15th birthday. Without `crossings`, 
#'   a monthly grid is used instead.
#' * A rate that is `Inf` somewhere in the window but nowhere the forced
#'   finder looks (`bind_hazard()`) is an error.
#' * The waiting times drawn for different destinations use independent
#'   random budgets. The independence is a computational device, not a
#'   modelling assumption, because the joint law of waiting
#'   time and destination depends only on the cause-specific rates.
#'
#' @param lambda Vectorised hazard: a function of `x`, the number of
#'   **years** elapsed since now, returning a non-negative numeric vector
#'   the same length as `x`. Built by `bind_hazard()`, which has already
#'   bound the individual's current age, calendar time and duration.
#' @param horizon_years Remaining observation window in **years**: the
#'   time to whichever comes first, the end of the simulation or the
#'   individual reaching `maxAge`. A non-positive value returns `Inf` and
#'   is not an error, because the engine legitimately asks about
#'   individuals who have run out of window. The budget is still used.
#' @param draw_u Function of no arguments returning one standard
#'   exponential variate, used when `budget` is `NULL`. Defaults to
#'   `default_draw_u()`. It is an argument for two reasons: tests can supply
#'   a fixed value.
#' @param forced_time Function of the remaining window returning the
#'   earliest forced point in it, or `NA`; built by `bind_hazard()`. `NULL`
#'   skips the check. The registry passes `NULL` when a rate is declared
#'   `may_be_infinite = FALSE`, which saves tens of thousands of rate 
#'   evaluations per scheduling call.
#' @param crossings Function of the remaining window returning the points
#'   at which a clock the rate uses reaches a whole number; built by
#'   `bind_hazard()`. `NULL` falls back to a monthly grid for the all-zero
#'   test.
#' @param label Human-readable name of the transition being solved, used
#'   only in error messages. The engine passes the edge's label, e.g.
#'   `"f/DE/0 -> f/DE/1"`, so that a bad rate names itself.
#' @param budget The standard exponential budget for this destination,
#'   drawn by the engine from the individual's stream. `NULL` draws one with
#'   `draw_u`. Every call uses exactly one budget, whichever path it takes.
#'
#' @return A list with
#'   \describe{
#'     \item{`time`}{waiting time in years, or `Inf` if no event occurs
#'       within `horizon_years`}
#'     \item{`u`}{the random budget that was drawn, or `NA_real_` for a forced
#'       transition}
#'     \item{`reliable`}{`FALSE` if any quadrature step was degraded,
#'       including those inside the root search}
#'     \item{`forced`}{`TRUE` if `time` came from an `Inf` rate rather than
#'       from the inversion}
#'   }
#'
#' @seealso `bind_hazard()`, which builds `lambda`, and [MicSim-concepts].
#'
#' @examples
#' # A constant hazard of 0.1 per year. The waiting time is then
#' # Exponential(0.1), so tau = u / 0.1 exactly.
#' constant <- function(x) rep(0.1, length(x))
#' draw_waiting_time(constant, horizon_years = 100, draw_u = function() 1)$time
#' # -> 10
#'
#' @keywords internal
#' @noRd
draw_waiting_time <- function(lambda,
                              horizon_years,
                              draw_u     = default_draw_u,
                              forced_time = NULL,
                              crossings   = NULL,
                              label       = "this transition",
                              budget      = NULL) {

  # ------------------------------------------------------------------
  # The random risk budget: a standard exponential variate. The engine
  # supplies it; otherwise it is drawn here, once, before anything else,
  # so that every call uses exactly one whatever path it then takes.
  # ------------------------------------------------------------------
  u <- if (is.null(budget)) draw_u() else budget

  # ------------------------------------------------------------------
  # Step (0). No window left.
  # The individual has reached maxAge, or the simulation horizon closes
  # today (e.g. an immigrant arriving on the last day). Nothing can happen;
  # the budget is spent and discarded.
  # ------------------------------------------------------------------
  if (horizon_years <= 0) {
    return(no_event(u = u, reliable = TRUE))
  }

  # ------------------------------------------------------------------
  # Step (1). Forced transitions.
  # A rate of Inf means "this event happens on a known date". The
  # canonical case is elementary-school enrolment: the user writes a rate
  # that returns Inf at age 7 and 0 elsewhere. Where the finder looks is
  # documented on bind_hazard().
  #
  # The budget is not needed on this path and is discarded.
  # ------------------------------------------------------------------
  if (!is.null(forced_time)) {
    at <- forced_time(horizon_years)
    if (!is.na(at)) {
      return(list(
        time     = at,
        u        = NA_real_,
        reliable = TRUE,
        forced   = TRUE
      ))
    }
  }

  # ------------------------------------------------------------------
  # Step (2). All-zero shortcut.
  # Sampling a fertility waiting time for a 70-year-old means integrating
  # a function that is zero everywhere. Quadrature would do that
  # correctly but at real cost, and this is by far the most common case
  # in a full life-course model. Tested at whole years only; see the
  # Assumptions section above.
  # ------------------------------------------------------------------
  checkpoints <- zero_test_points(horizon_years, crossings)
  rate_at_checkpoints <- lambda(checkpoints)
  assert_rate_values(rate_at_checkpoints, length(checkpoints), label)
  if (any(is.infinite(rate_at_checkpoints))) abort_inf_between(label)
  if (sum(rate_at_checkpoints) == 0) {
    return(no_event(u = u, reliable = TRUE))
  }

  # From here on, the rate is integrated. An Inf the forced finder did not
  # see makes integrate() fail; integrate_hazard() turns that into an
  # error that names the rate.
  integrand <- lambda

  # ------------------------------------------------------------------
  # Step (3). Total risk accumulated over the whole remaining window.
  # If the drawn budget exceeds it, the budget is never exhausted and no
  # event occurs before the window closes.
  # ------------------------------------------------------------------
  total <- integrate_hazard(integrand, 0, horizon_years, label)
  if (total$value <= 0 || u > total$value) {
    return(no_event(u = u, reliable = total$reliable))
  }

  # ------------------------------------------------------------------
  # Step (4). Invert: find tau with integral from 0 to tau of lambda = u.
  #
  # Define f(x) = Lambda(x) - u. Then f(0) = -u < 0 and
  # f(horizon) = total - u > 0, so a root is bracketed and uniroot()
  # cannot fail for want of a sign change. We hand those two values to
  # uniroot() explicitly via f.lower / f.upper, which saves two full
  # quadratures per event; each one costs as much as everything above.
  # ------------------------------------------------------------------
  # Each step of the root search integrates again; whether any of those
  # integrations had to fall back is recorded in `inner`, so that the
  # reliability flag covers the whole inversion, not only the total.
  inner <- new.env(parent = emptyenv())
  inner$reliable <- TRUE
  accumulated_minus_budget <- function(x) {
    step <- integrate_hazard(integrand, 0, x, label)
    if (!step$reliable) inner$reliable <- FALSE
    step$value - u
  }

  root <- tryCatch(
    stats::uniroot(
      accumulated_minus_budget,
      lower   = 0,
      upper   = horizon_years,
      f.lower = -u,
      f.upper = total$value - u,
      tol     = UNIROOT_TOL
    ),
    error = function(e) e
  )

  if (inherits(root, "micsim_rate_error")) stop(root)
  if (inherits(root, "condition")) {
    # A bracketed root that uniroot() still cannot find means the
    # integrated hazard is not behaving monotonically under quadrature,
    # i.e. the numbers are untrustworthy. Producing a plausible-looking
    # wrong waiting time here would be worse than stopping.
    abort_solver(
      paste0(
        "Could not invert the integrated hazard.\n",
        "  The root search failed on a bracket that contains a root, which\n",
        "  points to a numerically pathological transition rate.\n",
        "  Underlying message: ", conditionMessage(root)
      ),
      horizon_years = horizon_years,
      u = u,
      integrated_hazard = total$value
    )
  }

  list(
    time     = root$root,
    u        = u,
    reliable = total$reliable && inner$reliable,
    forced   = FALSE
  )
}


#' Where the all-zero shortcut samples a rate
#'
#' @param horizon_years The remaining window.
#' @param crossings `NULL`, or a function returning the whole-number
#'   crossings of the rate's clocks in the window.
#' @return Sorted, unique points in `[0, horizon_years]`.
#' @keywords internal
#' @noRd
zero_test_points <- function(horizon_years, crossings) {
  H <- horizon_years
  inside <- if (is.null(crossings)) {
    seq(0, H, by = 1 / 12)
  } else {
    # Just AFTER each crossing, so that a table indexed by as.integer(age)
    # is read in the new cell rather than, by rounding, the old one.
    pmin(crossings(H) + 1e-9, H)
  }
  sort(unique(c(0, seq_len(floor(H)), inside, H)))
}


#' Raise the error for an Inf rate the forced finder cannot place
#'
#' @param label The transition's label.
#' @return Never returns.
#' @keywords internal
#' @noRd
abort_inf_between <- function(label) {
  abort_rate(label, paste0(
    "It returns Inf at a point the simulation cannot place as a forced\n",
    "  transition. A forced (Inf) rate must be Inf at a whole age, calendar\n",
    "  year or duration (e.g. ifelse(age == 7, Inf, 0)), or on an interval\n",
    "  of at least a day (e.g. ifelse(age >= 7, Inf, 0))."
  ))
}


#' Result when no event occurs
#'
#' Four places in `draw_waiting_time()` return the same thing. Naming
#' it makes each of them read as a statement about the model rather than
#' as list assembly.
#'
#' @param u The budget that was drawn, or `NA_real_` if none was.
#' @param reliable Whether the quadrature behind the decision was reliable.
#' @return A waiting-time result list with `time = Inf`.
#' @keywords internal
#' @noRd
no_event <- function(u, reliable) {
  list(time = Inf, u = u, reliable = reliable, forced = FALSE)
}


#' Check the values a rate returned
#'
#' @description
#' Guards the two places where the solver evaluates a user-supplied rate
#' and would otherwise fail later with an opaque message from
#' [stats::integrate()].
#'
#' @details
#' The registry checks the shape of every rate once at construction, and
#' the model validator checks values over the model's full age and
#' calendar grid. Neither can cover every combination of arguments the engine
#' actually reaches, so these two checks run on values the solver has
#' computed anyway, at no extra cost.
#'
#' Note that `Inf` is permitted, because it encodes a forced transition. 
#' `NA` and negative values are not, because both make the integrated hazard 
#' meaningless.
#'
#' @param values The numeric vector the rate returned.
#' @param expected_length How many values were asked for.
#' @param label Human-readable transition name, used in the message.
#' @return `invisible(TRUE)`, or raises `micsim_rate_error`.
#' @keywords internal
#' @noRd
assert_rate_values <- function(values, expected_length, label) {
  if (!is.numeric(values)) {
    abort_rate(label, paste0("It returned ", class(values)[1L],
                             " instead of a numeric vector."))
  }
  if (length(values) != expected_length) {
    abort_rate(label, paste0(
      "It is not vectorised: asked for ", expected_length,
      " value(s), got ", length(values), ".\n",
      "  Rate functions must return one value per element of their input.\n",
      "  A constant rate must therefore be written as, for example,\n",
      "  function(age, calTime) rep(0.1, length(age))."
    ))
  }
  if (anyNA(values)) {
    abort_rate(label, "It returned NA. The integrated hazard is then undefined.")
  }
  if (any(values < 0)) {
    abort_rate(label, "It returned a negative value. Rates cannot be negative.")
  }
  invisible(TRUE)
}
