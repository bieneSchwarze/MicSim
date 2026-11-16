# =====================================================================
# registry.R
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   Turns the user's description of "which transitions are possible and at
#   what rate" into a lookup table the engine can use without thinking.
#
# PLACE IN THE PIPELINE
#   codec.R --(state indices)--> registry.R --(hazard closures)--> solver.R
#
#   The registry is built once, before the simulation starts. During the
#   run the engine asks it two questions and nothing else:
#     feasible(i)                 -> which transitions leave state i?
#     bind_hazard(edge, ...)      -> give me that transition's rate as a
#                                    one-argument function of elapsed years
#
# WHAT YOU NEED TO KNOW  (concepts in ?MicSim-concepts)
#   transition rate, the three time scales (age / calTime / duration),
#   composite state, absorbing state.
#
#   Rates may depend on age, calendar time (calTime) and time since a state
#   variable last changed (duration). bind_hazard() fixes their current
#   values and returns a function of elapsed years for the solver. The solver
#   then does not need to know anything about which time scales the underlying
#   rate actually depends on.
#
# GUARANTEES
#   * Every rate is a resolved function object, never a name looked up at
#     run time.
#   * Every edge knows, before the run starts, which time scales its rate
#     uses, which state variable it changes, and whether it can be forced.
#   * Nothing here is mutable. The registry is read-only for the engine.
# =====================================================================


#' Build the transition-rate registry
#'
#' @description
#' Expands transition specifications into individual origin-destination pairs
#' and groups them by origin state for quick lookup during a simulation.
#'
#' @details
#' A transition is a triple: origin state, destination state, and the
#' function giving the rate at which the move happens. For each triple the
#' registry records, once and for all:
#'
#' \describe{
#'   \item{the rate closure}{the function object itself, so nothing is ever
#'     looked up by name during the run}
#'   \item{its time-scale dependence}{derived from the function's argument
#'     names via [base::formals()]. A rate declaring `function(age,
#'     calTime)` depends on age and calendar time and not on duration.}
#'   \item{the affected state variable}{which one of the composite state's
#'     components this transition changes. If a rate uses `duration`, it
#'     uses the time since this component last changed. For a smoking
#'     relapse rate, the relevant variable could be `smoking_status`.}
#'   \item{whether the rate can be infinite}{an `Inf` rate encodes a forced
#'     event on a known date. Detecting one costs a rate evaluation at
#'     every day of the remaining horizon, so the engine wants to know in
#'     advance which edges can skip that check.}
#' }
#'
#' @section Assumptions:
#' * A rate function can take no arguments or any combination of `age`,
#'   `calTime` and `duration`. For a constant rate, use a function that
#'   takes no arguments and returns a single number, such as
#'   `function() 0.05`. The engine rejects any other argument names.
#' * If a transition changes two or more state variables at once, the
#'   notion of "time since the last change of the affected variable" is
#'   ill-defined. Here: if the rate uses `duration` the model is rejected;
#'   if it does not, the situation is harmless and is allowed through.
#' * Transitions into absorbing states have no single affected variable
#'   either: the codec knows which states are absorbing, whatever they are
#'   called. They receive `duration = Inf`.
#'
#' @param transitions A list of transitions. Each element is a list with
#'   components `from`, `to` and `rate`, or `move` (a [move()] selector)
#'   and `rate`. `from` may be `NA` as shorthand for "every transient
#'   state except to", which is how mortality is normally declared. `rate`
#'   must be a function, not a name.
#' @param codec A `state_codec`, from `new_state_codec()`. It supplies the
#'   enumeration of states and the label-to-index mapping.
#' @param probe Optional grid from `micsim_rate_probe()`, with components
#'   `ages`, `cal_years` and `durations`. `NULL` keeps every edge marked as
#'   possibly infinite.
#'
#' @return An object of class `rate_registry`: a list with
#'   \describe{
#'     \item{`feasible`}{a function taking a state index and returning its
#'       outgoing edges, or `NULL` if there are none}
#'     \item{`edges`}{the flat list of all edges}
#'     \item{`codec`}{the supplied state codec}
#'   }
#'
#' @seealso `bind_hazard()`, `micsim_rate_probe()`, [MicSim-concepts]
#'
#' @examples
#' codec <- new_state_codec(
#'   state_space = list(sex = c("m", "f"), fert = c("0", "1")),
#'   absorbing   = "dead"
#' )
#' mortality <- function(age, calTime) 0.00003 * exp(0.1 * age)
#' fertility <- function(age, calTime) ifelse(age > 15 & age < 45, 0.1, 0)
#' registry <- new_rate_registry(
#'   transitions = list(
#'     list(from = NA,      to = "dead",  rate = mortality),
#'     list(from = "f/0",   to = "f/1",   rate = fertility)
#'   ),
#'   codec = codec
#' )
#' registry
#'
#' @keywords internal
#' @noRd
new_rate_registry <- function(transitions, codec, probe = NULL) {

  # --- guard the inputs ----------------------------------------------
  if (!inherits(codec, "state_codec")) {
    micsim_abort(
      "micsim_registry_error",
      "`codec` must be a state_codec, as returned by new_state_codec()."
    )
  }
  if (!is.list(transitions) || length(transitions) == 0L) {
    micsim_abort(
      "micsim_registry_error",
      "`transitions` must be a non-empty list of transition specifications."
    )
  }

  # --- expand and resolve every transition into one or more edges -----
  # from = NA expands to one edge per transient origin; move() expands to
  # one per matching origin. Both can produce many edges from one declaration.
  edges <- list()
  # Cache probe results by rate function, since many edges can share a rate.
  cache <- new.env(parent = emptyenv())
  cache$rates <- list()
  cache$answers <- logical(0)
  infinite_for <- function(rate) {
    if (is.null(probe) || !is.function(rate)) return(TRUE)
    for (k in seq_along(cache$rates)) {
      if (identical(cache$rates[[k]], rate)) return(cache$answers[k])
    }
    argnames <- names(formals(rate))
    if (is.null(argnames)) argnames <- character(0)
    answer <- if (all(argnames %in% RATE_ARGS)) {
      rate_may_be_infinite(rate, argnames, probe)
    } else {
      TRUE                               # new_rate_edge() reports the error
    }
    cache$rates[[length(cache$rates) + 1L]] <- rate
    cache$answers <- c(cache$answers, answer)
    answer
  }
  for (d in seq_along(transitions)) {
    tr <- transitions[[d]]
    pairs <- expand_transition(tr, codec)
    infinite <- if (nrow(pairs) > 0L) infinite_for(tr$rate) else TRUE
    for (k in seq_len(nrow(pairs))) {
      e <- new_rate_edge(unname(pairs[k, 1L]), unname(pairs[k, 2L]), tr$rate,
                         codec, may_be_infinite = infinite)
      # Retain the original declaration for reporting and reestimation.
      e$decl <- d
      edges[[length(edges) + 1L]] <- e
    }
  }

  # Give each edge its list position as an ID for subsystem selections.
  for (k in seq_along(edges)) edges[[k]]$id <- k

  # Repeated pairs would be scheduled as separate competing risks, adding
  # their hazards. Reject them to avoid counting the same transition twice.
  pair_key <- vapply(edges, function(e) paste(e$from_idx, e$to_idx),
                     character(1L))
  twice <- unique(pair_key[duplicated(pair_key)])
  if (length(twice) > 0L) {
    labels <- vapply(edges[match(twice, pair_key)], `[[`, character(1L),
                     "label")
    micsim_abort(
      "micsim_registry_error",
      paste0(
        "Transition(s) declared more than once: ",
        toString(utils_head(labels, 5L)),
        if (length(labels) > 5L) paste0(" and ", length(labels) - 5L, " more") else "",
        ".\n",
        "  A pair declared twice would be drawn as two independent risks,\n",
        "  doubling its hazard. Declare each transition once (note that\n",
        "  from = NA already covers every transient state)."
      )
    )
  }

  # --- index the edges by origin state -------------------------------
  # Allocate each origin's edge list at its final size, then fill it.
  # The engine can then retrieve outgoing edges directly by state index.
  origin_of <- vapply(edges, function(e) e$from_idx, integer(1L))
  counts <- tabulate(origin_of, nbins = codec$n)

  by_origin <- vector("list", codec$n)
  for (i in which(counts > 0L)) {
    by_origin[[i]] <- vector("list", counts[i])
  }
  filled <- integer(codec$n)
  for (k in seq_along(edges)) {
    i <- origin_of[k]
    filled[i] <- filled[i] + 1L
    by_origin[[i]][[filled[i]]] <- edges[[k]]
  }

  # feasible(i): the outgoing edges of state index i, or NULL if none.
  # States without outgoing edges, including absorbing states, return NULL.
  feasible <- function(state_idx) by_origin[[state_idx]]

  structure(
    list(feasible = feasible, edges = edges, codec = codec,
         probed = !is.null(probe)),
    class = "rate_registry"
  )
}


#' Resolve one transition specification into (origin, destination) pairs
#'
#' Accepts a [move()] selector, explicit origin labels, or `from = NA`
#' for every transient state. Checks that whatever the user wrote names
#' states that exist.
#'
#' @param tr One transition specification: a list with `from`, `to` and
#'   `rate`, or with `edge` and `rate`.
#' @param codec A `state_codec`.
#' @return A two-column character matrix of origin and destination labels.
#' @keywords internal
#' @noRd
expand_transition <- function(tr, codec) {
  if (is.null(tr$rate)) {
    micsim_abort(
      "micsim_registry_error",
      "Every transition must have `from`, `to` and `rate` components."
    )
  }
  if (!is.null(tr$move)) {
    if (!inherits(tr$move, "micsim_move")) {
      micsim_abort(
        "micsim_registry_error",
        "A transition can be declared with move(), but not with moves_of()."
      )
    }
    return(resolve_move(tr$move, codec))
  }
  if (is.null(tr$to) || length(tr$to) != 1L || is.na(tr$to)) {
    micsim_abort(
      "micsim_registry_error",
      "Every transition must have `from`, `to` and `rate` components, with\n  one destination state."
    )
  }
  to <- as.character(tr$to)
  origins <- expand_origins(tr, codec)
  # A transition must change the state. With from = NA ("every transient
  # state"), the destination itself is left out of the origins; listing
  # it explicitly is an error. Otherwise x -> x would be simulated as an
  # event that changes nothing but resets x's duration clock.
  if (length(tr$from) == 1L && is.na(tr$from)) {
    origins <- setdiff(origins, to)
    if (length(origins) == 0L) {
      micsim_abort(
        "micsim_registry_error",
        paste0("A transition from every state to \"", to, "\" has no origin:\n",
               "  \"", to, "\" is the only transient state, and a transition must\n",
               "  change the state.")
      )
    }
  } else if (to %in% origins) {
    micsim_abort(
      "micsim_registry_error",
      paste0("A transition must change the state, but \"", to, "\" is both an\n",
             "  origin and the destination. Remove it from `from`.")
    )
  }
  cbind(from = origins, to = rep(to, length(origins)))
}


#' Resolve the origin states of one transition specification
#'
#' Handles the `from = NA` shorthand and checks that whatever the user
#' wrote names states that actually exist.
#'
#' @param tr One transition specification: a list with `from`, `to`, `rate`.
#' @param codec A `state_codec`.
#' @return Character vector of origin state labels.
#' @keywords internal
#' @noRd
expand_origins <- function(tr, codec) {
  if (is.null(tr$to) || is.null(tr$rate)) {
    micsim_abort(
      "micsim_registry_error",
      "Every transition must have `from`, `to` and `rate` components."
    )
  }
  if (length(tr$from) == 1L && is.na(tr$from)) {
    return(codec$transient)
  }
  if (anyNA(tr$from)) {
    micsim_abort(
      "micsim_registry_error",
      "`from` is either NA, meaning every transient state, or state labels; not a mix."
    )
  }
  unknown <- setdiff(tr$from, codec$transient)
  if (length(unknown) > 0L) {
    micsim_abort(
      "micsim_registry_error",
      paste0("Unknown origin state(s): ", toString(unknown), ".",
             if (any(unknown %in% codec$absorbing))
               "\n  Absorbing states cannot be left, so they cannot be origins." else "")
    )
  }
  as.character(tr$from)
}


#' Build one resolved transition edge
#'
#' @description
#' An edge is everything the engine needs to know about a single
#' (origin, destination) pair, computed once so that nothing has to be
#' recomputed during the run.
#'
#' @param from_label,to_label State labels, e.g. `"f/DE/0"` and `"dead"`.
#' @param rate The user's rate function.
#' @param codec A `state_codec`.
#' @param may_be_infinite Whether the rate can return `Inf`, as probed by
#'   `new_rate_registry()`; `TRUE`, the safe answer, by default.
#'
#' @return A list with the fields documented in `new_rate_registry()`
#' @keywords internal
#' @noRd
new_rate_edge <- function(from_label, to_label, rate, codec,
                          may_be_infinite = TRUE) {

  label <- paste(from_label, "->", to_label)

  if (!is.function(rate)) {
    abort_rate(label, paste0(
      "It is a ", class(rate)[1L], ", not a function.\n",
      "  Pass the rate function itself, not its name as a string."
    ))
  }

  argnames <- names(formals(rate))
  if (is.null(argnames)) argnames <- character(0)
  illegal <- setdiff(argnames, RATE_ARGS)
  if (length(illegal) > 0L) {
    abort_rate(label, paste0(
      "It declares argument(s) ", toString(illegal), ".\n",
      "  Rate functions may take only ", toString(RATE_ARGS),
      ", in any subset and any order."
    ))
  }

  uses_duration <- "duration" %in% argnames
  from_idx <- codec$encode(from_label)
  to_idx   <- codec$encode(to_label)

  # Absorbing destinations are recognised by the codec, not by their
  # label. They have no affected variable, so their duration is Inf. 
  to_absorbing <- codec$is_absorbing(to_idx)
  changed <- if (to_absorbing) integer(0) else
    changed_state_variables(from_label, to_label)
  affected <- if (length(changed) == 1L) changed else NA_integer_

  # A duration-dependent rate needs an unambiguous clock to measure
  # against. If this transition changes several state variables at once
  # there is no such clock, so a rate that depends on duration here is
  # rejected rather than silently given duration = Inf.
  if (uses_duration && is.na(affected) && !to_absorbing) {
    abort_rate(label, paste0(
      "It depends on `duration`, but this transition changes more than one\n",
      "  state variable at once, so 'time since the last change of the\n",
      "  affected variable' is not defined. Split it into single-variable\n",
      "  transitions, or drop the duration dependence."
    ))
  }

  assert_rate_shape(rate, argnames, label)

  list(
    from_label    = from_label,
    to_label      = to_label,
    label         = label,
    from_idx      = from_idx,
    to_idx        = to_idx,
    rate          = rate,
    argnames      = argnames,
    uses_age      = "age" %in% argnames,
    uses_caltime  = "calTime" %in% argnames,
    uses_duration = uses_duration,
    affected_var  = affected,
    changed_vars  = changed,
    may_be_infinite = isTRUE(may_be_infinite)
  )
}


#' Find all state variables changed by a transition
#'
#' @description
#' Compares the origin and destination labels component by component and
#' returns the positions of the components that differ.
#'
#' @details
#' Compares slash-separated labels component by component. The returned
#' positions identify the duration clocks to reset after a transition.
#' Returns `integer(0)` when the labels have different numbers of
#' components, which is the case for a transition into an absorbing state
#' such as `"dead"`.
#'
#' @param from_label,to_label Composite state labels.
#' @return Integer vector of positions, possibly empty.
#' @keywords internal
#' @noRd
changed_state_variables <- function(from_label, to_label) {
  from_parts <- strsplit(from_label, "/", fixed = TRUE)[[1L]]
  to_parts   <- strsplit(to_label,   "/", fixed = TRUE)[[1L]]
  if (length(from_parts) != length(to_parts)) {
    return(integer(0))
  }
  as.integer(which(from_parts != to_parts))
}


#' Check that a rate function has the right shape
#'
#' @description
#' Calls the rate once on a short vector and checks that it returns one
#' number per input. This catches the single most common user error, a
#' rate written as `function(age) 0.1` instead of
#' `function(age) rep(0.1, length(age))`
#'
#' @details
#' Only the *shape* is checked here, not the values. Checking values requires
#' knowing the model's age range and simulation horizon, which the registry
#' does not have; `validate()` does that on the model's grid
#' (see `check_one_rate()` in R/validate.R).
#'
#' @section Assumptions:
#' The registry does not know which ages the model can reach, so it cannot
#' know which ages a rate has to work at. A perfectly good rate may be
#' defined only from age 17 -- a table of adult transition rates, say,
#' whose first row is age 17 and which is indexed as `vec[age - 16]`. Ask
#' it about age 0 and it fails, correctly.
#'
#' So the probe tries several sets of ages and uses the first that
#' evaluates. If none of them does, the check is SKIPPED rather than
#' failed. 
#'
#' @param rate The rate function.
#' @param argnames Its argument names.
#' @param label Human-readable transition label, for the error message.
#' @return `invisible(TRUE)`, or raises `micsim_rate_error`.
#' @keywords internal
#' @noRd
assert_rate_shape <- function(rate, argnames, label) {

  if (length(argnames) == 0L) {
    # A rate declaring no arguments is a constant. It is allowed, and
    # bind_hazard() recycles it, so only the type matters here.
    out <- tryCatch(rate(), error = function(e) e)
    if (inherits(out, "condition")) {
      abort_rate(label, paste0("Evaluating it raised an error: ",
                               conditionMessage(out)))
    }
    if (!is.numeric(out) || length(out) != 1L) {
      abort_rate(label, "A rate taking no arguments must return one number.")
    }
    return(invisible(TRUE))
  }

  # Candidate probes, tried in order. Working ages first, because a rate
  # table that does not start at zero is common and one that does not
  # cover working ages is not.
  candidates <- list(
    list(age = c(20, 30, 40), calTime = c(2000, 2010, 2020),
         duration = c(0, 1, 2)),
    list(age = c(0, 30, 60),  calTime = c(2000, 2010, 2020),
         duration = c(0, 1, 2)),
    list(age = c(50, 60, 70), calTime = c(1990, 2000, 2010),
         duration = c(0, 5, 10)),
    list(age = c(1, 2, 3),    calTime = c(2015, 2016, 2017),
         duration = c(0, 1, 2))
  )

  for (probe in candidates) {
    args <- probe[intersect(names(probe), argnames)]
    out <- tryCatch(do.call(rate, args), error = function(e) e)
    if (inherits(out, "condition")) next          # not defined here; try again

    if (!is.numeric(out) || length(out) != 3L) {
      abort_rate(label, paste0(
        "It is not vectorised: given 3 value(s) per argument it returned ",
        length(out), ".\n",
        "  Rate functions must return one value per element of their input,\n",
        "  because the simulation integrates them numerically.\n",
        "  A constant rate is written function(age, calTime) rep(0.1, length(age))."
      ))
    }
    return(invisible(TRUE))
  }

  # The rate could not be evaluated at any of the probe points. That is
  # not evidence of a fault: see the Assumptions section.
  invisible(TRUE)
}


#' Prepare a grid for checking rate functions
#'
#' @description
#' Builds the argument grid used by `rate_may_be_infinite()` and
#' `validate()` (via `check_one_rate()`). Kept as its own object so
#' that the model layer, which knows `maxAge` and the simulation
#' horizon, can hand the same grid to both.
#'
#' @param ages Integer ages to test, typically `0:maxAge`.
#' @param cal_years Integer calendar years to test, typically the years
#'   spanned by the simulation horizon.
#' @param durations Integer durations to test, typically `0:maxAge`.
#' @param max_points Cap on the total number of rate evaluations. Rates
#'   using all three time scales would otherwise need
#'   `length(ages) * length(cal_years) * length(durations)` evaluations,
#'   which for a 120-year, 100-year model is 1.5 million. When the cap
#'   would be exceeded the duration grid is thinned, because duration is
#'   the scale on which rates are smoothest in practice.
#'
#' @return A list with components `ages`, `cal_years`, `durations`.
#' @seealso `new_rate_registry()`
#' @keywords internal
#' @noRd
micsim_rate_probe <- function(ages, cal_years, durations,
                              max_points = 1e6) {
  total <- length(ages) * length(cal_years) * length(durations)
  if (total > max_points && length(durations) > 1L) {
    keep <- max(2L, floor(max_points / (length(ages) * length(cal_years))))
    idx <- unique(round(seq(1, length(durations), length.out = keep)))
    durations <- durations[idx]
  }
  list(ages = ages, cal_years = cal_years, durations = durations)
}


#' Check whether a rate may return `Inf`
#'
#' @description
#' An `Inf` rate marks a forced transition. Looking for one can involve
#' checking the rate at every day in the remaining simulation window. A grid
#' check lets the solver skip that search for rates with no infinite values
#' on the grid.
#'
#' @details
#' The check evaluates all combinations of grid values for the arguments
#' the rate uses. A rate with no arguments is evaluated once. When the
#' check cannot tell whether a rate may be infinite, the solver still
#' searches for a forced transition.
#'
#' @param rate The rate function.
#' @param argnames Its argument names.
#' @param probe A grid from `micsim_rate_probe()`, or `NULL`.
#' @return `TRUE` or `FALSE`.
#' @keywords internal
#' @noRd
rate_may_be_infinite <- function(rate, argnames, probe = NULL) {
  if (is.null(probe)) {
    return(TRUE)
  }
  # A rate that cannot be evaluated on the whole lattice -- a table that
  # starts at age 17, say -- may still be Inf somewhere it can be
  # evaluated, so it keeps the safe answer rather than failing here.
  values <- tryCatch(evaluate_rate_on_grid(rate, argnames, probe),
                     error = function(e) NULL)
  if (is.null(values) || !is.numeric(values)) return(TRUE)
  any(is.infinite(values))
}


#' Build the model's whole-year grid for detecting infinite rates
#'
#' @description
#' Includes every whole age and duration from zero up to, but below,
#' `max_age`. It also includes every calendar year from the start to the
#' end of the simulation, including both endpoint years. Ages start at
#' zero because a fertility subsystem can introduce newborns during the run.
#'
#' The grid keeps every whole-year value so the check can detect a rate
#' that is infinite at a single value, such as when `duration == 5`.
#'
#' @param max_age Sharp maximum age, in years.
#' @param horizon Length-2 vector of `yyyymmdd` dates.
#' @return A grid from `micsim_rate_probe()`.
#' @keywords internal
#' @noRd
forced_probe_grid <- function(max_age, horizon) {
  top   <- max(0L, as.integer(ceiling(max_age)) - 1L)
  years <- getYear(getInDays(horizon))
  if (anyNA(years)) return(NULL)

  micsim_rate_probe(ages = 0:top, cal_years = seq.int(years[1L], years[2L]),
                    durations = 0:top, max_points = Inf)
}


#' Evaluate a rate over a probe grid
#'
#' Calls the rate once with all combinations of the grid values for its
#' arguments. For example, a rate that takes `age` and `calTime` receives
#' vectors containing every pair of ages and calendar years. A constant
#' rate with no arguments, such as `function() 0.05`, is called directly.
#'
#' @param rate The rate function.
#' @param argnames Its argument names.
#' @param probe A grid from `micsim_rate_probe()`.
#' @return Numeric vector of rate values.
#' @keywords internal
#' @noRd
evaluate_rate_on_grid <- function(rate, argnames, probe) {
  if (length(argnames) == 0L) {
    return(as.numeric(rate()))
  }
  pieces <- list()
  if ("age"      %in% argnames) pieces$age      <- probe$ages
  if ("calTime"  %in% argnames) pieces$calTime  <- probe$cal_years
  if ("duration" %in% argnames) pieces$duration <- probe$durations

  lattice <- expand.grid(pieces, KEEP.OUT.ATTRS = FALSE)
  do.call(rate, as.list(lattice))
}


#' Describe the inputs at one point on the rate grid
#'
#' Rebuilds the grid `evaluate_rate_on_grid()` used and names the
#' arguments at one position, so that a failing rate says WHERE it failed
#' rather than only that it did.
#'
#' @param argnames The rate's argument names.
#' @param probe A grid from `micsim_rate_probe()`.
#' @param index Position in the evaluated vector.
#' @return A single string, e.g. `"age 0, calTime 2014"`.
#' @keywords internal
#' @noRd
describe_grid_point <- function(argnames, probe, index) {
  if (length(argnames) == 0L) return("a constant rate")
  pieces <- list()
  if ("age"      %in% argnames) pieces$age      <- probe$ages
  if ("calTime"  %in% argnames) pieces$calTime  <- probe$cal_years
  if ("duration" %in% argnames) pieces$duration <- probe$durations
  lattice <- expand.grid(pieces, KEEP.OUT.ATTRS = FALSE)
  row <- lattice[index, , drop = FALSE]
  paste(names(row), unlist(row), collapse = ", ")
}


#' Express a transition rate as a function of time from now
#'
#' @description
#' Builds a function that gives the transition rate at a specified number
#' of years from now. The solver uses this function to calculate a waiting
#' time for the transition.
#'
#' @details
#' All three time scales advance together as the individual waits. If the
#' person is 34.2 years old today, then `x` years from now they are
#' `34.2 + x` years old, the calendar reads `caltime0_years + x`, and they
#' have been in the affected state for `duration0_years + x` years. So
#' binding the three starting values is enough to reduce the rate to
#' \eqn{\lambda(x)}.
#'
#' Two closures are returned:
#'
#' \describe{
#'   \item{`lambda`}{the exact hazard, used for integration and inversion}
#'   \item{`forced_time`}{a function of the remaining window `H` (years)
#'     that returns the earliest `x` in `[0, H]` at which the rate is
#'     `Inf` -- a forced transition -- or `NA` if there is none. `NULL`
#'     when the registry has established that the rate can never be
#'     infinite.}
#' }
#'
#' **Where a forced transition is looked for.** Users write forced rates
#' as `ifelse(age == 7, Inf, 0)`: an exact equality on a whole number,
#' which a scan of fractional ages essentially never hits. `forced_time()`
#' therefore evaluates the rate at three kinds of point, with the time
#' scales at their exact values:
#'
#' \describe{
#'   \item{now, in completed years}{`x = 0`, with `trunc(age)`,
#'     `trunc(calTime)` and `trunc(duration)`. A child aged 7.4 still in
#'     "no school" with an `age == 7` rate is enrolled straight away}
#'   \item{every whole-number crossing}{each `x` at which age, calendar
#'     time or duration (whichever the rate uses) reaches a whole number,
#'     with that time scale set to exactly that number. A newborn with an
#'     `age == 7` rate is forced at exactly age 7.}
#'   \item{every day}{`x = k / 365.25`, at exact values, which catches a
#'     rate that is `Inf` over an interval such as `age > 7`.}
#' }
#'
#' @param edge One edge from a `rate_registry`.
#' @param age0_years The individual's age now, in fractional years.
#' @param caltime0_years Calendar time now, in fractional years
#'   (e.g. 2014.0 for 1 January 2014).
#' @param duration0_years Years since the last change of the state variable
#'   this edge changes. `Inf` when that is undefined, which is the case for
#'   transitions into absorbing states and for multi-variable transitions.
#'
#' @return A list with components `lambda`, `forced_time` and
#'   `crossings` (the points in a window at which a clock the rate uses
#'   reaches a whole number), all suitable for passing to
#'   `draw_waiting_time()`.
#'
#' @seealso `draw_waiting_time()`, [MicSim-concepts]
#'
#' @examples
#' fert <- function(age, calTime) ifelse(age > 15 & age < 45, 0.1, 0)
#' edge <- list(rate = fert, argnames = c("age", "calTime"),
#'              uses_age = TRUE, uses_caltime = TRUE, uses_duration = FALSE,
#'              may_be_infinite = FALSE)
#' h <- bind_hazard(edge, age0_years = 30, caltime0_years = 2014,
#'                  duration0_years = Inf)
#' h$lambda(c(0, 5, 20))   # rate now, in 5 years, in 20 years
#'
#' @keywords internal
#' @noRd
bind_hazard <- function(edge,
                        age0_years,
                        caltime0_years,
                        duration0_years) {

  rate <- edge$rate

  # Evaluate a rate taking no arguments once, then return one value per x
  # as required by numerical integration.
  if (length(edge$argnames) == 0L) {
    constant <- as.numeric(rate())
    lambda <- function(x) rep(constant, length(x))
    forced_time <- if (isTRUE(edge$may_be_infinite)) {
      function(horizon_years) if (is.infinite(constant)) 0 else NA_real_
    } else {
      NULL
    }
    return(list(lambda = lambda, forced_time = forced_time,
                crossings = function(horizon_years) numeric(0)))
  }

  # The rate at explicitly given values of the three clocks.
  rate_at <- function(age, cal, dur) {
    args <- list()
    if (edge$uses_age)      args$age      <- age
    if (edge$uses_caltime)  args$calTime  <- cal
    if (edge$uses_duration) args$duration <- dur
    do.call(rate, args)
  }

  # The hazard proper, called many times per integration. Written out per
  # combination of clocks
  a0 <- age0_years; c0 <- caltime0_years; d0 <- duration0_years
  lambda <- switch(
    paste0(edge$uses_age, edge$uses_caltime, edge$uses_duration),
    "TRUEFALSEFALSE" = function(x) rate(age = a0 + x),
    "FALSETRUEFALSE" = function(x) rate(calTime = c0 + x),
    "FALSEFALSETRUE" = function(x) rate(duration = d0 + x),
    "TRUETRUEFALSE"  = function(x) rate(age = a0 + x, calTime = c0 + x),
    "TRUEFALSETRUE"  = function(x) rate(age = a0 + x, duration = d0 + x),
    "FALSETRUETRUE"  = function(x) rate(calTime = c0 + x, duration = d0 + x),
    function(x) rate(age = a0 + x, calTime = c0 + x, duration = d0 + x)
  )

  forced_time <- if (isTRUE(edge$may_be_infinite)) {
    function(horizon_years) {
      forced_probe(rate_at, edge, age0_years, caltime0_years,
                   duration0_years, horizon_years, edge$label)
    }
  } else {
    NULL
  }

  # Where a rate keyed on whole ages, years or durations can change: the
  # points at which one of the clocks it uses reaches a whole number. The
  # solver's all-zero test samples just after each of them.
  crossings <- function(horizon_years) {
    out <- numeric(0)
    for (v0 in c(if (edge$uses_age) age0_years,
                 if (edge$uses_caltime) caltime0_years,
                 if (edge$uses_duration) duration0_years)) {
      if (!is.finite(v0)) next
      first <- floor(v0) + 1
      last  <- floor(v0 + horizon_years)
      if (last >= first) out <- c(out, seq(first, last) - v0)
    }
    out
  }

  list(lambda = lambda, forced_time = forced_time, crossings = crossings)
}


#' Find the earliest sampled time with an infinite rate
#'
#' @description
#' Evaluates the rate, in one vectorised call, at the points described in
#' `bind_hazard()`: now in completed years, every whole-number crossing of
#' the time scales the rate uses, and every day of the window. Returns the
#' earliest point at which it is `Inf`.
#'
#' @param rate_at Function of `(age, cal, dur)` vectors returning rates.
#' @param edge The edge, for which time scales the rate uses.
#' @param age0,cal0,dur0 The clocks now, in years.
#' @param horizon_years The remaining window `H`.
#' @param label The transition's label, for error messages.
#' @return The earliest `x` in `[0, H]` at which the rate is `Inf`, or
#'   `NA_real_`.
#' @keywords internal
#' @noRd
forced_probe <- function(rate_at, edge, age0, cal0, dur0, horizon_years,
                         label) {

  H <- horizon_years

  # (1) Check now using completed years. trunc(Inf) preserves an undefined
  #     duration as Inf.
  x   <- 0
  age <- trunc(age0)
  cal <- trunc(cal0)
  dur <- trunc(dur0)

  # (2) Every day of the window, at exact values.
  days <- seq.int(0L, trunc(H * DAYS_PER_YEAR)) / DAYS_PER_YEAR
  x   <- c(x, days)
  age <- c(age, age0 + days)
  cal <- c(cal, cal0 + days)
  dur <- c(dur, dur0 + days)

  # (3) Every whole-number crossing of a time scale the rate uses, with
  #     that scale set to exactly the whole number.
  crossings <- function(v0) {
    if (!is.finite(v0)) return(numeric(0))
    first <- floor(v0) + 1
    last  <- floor(v0 + H)
    if (last < first) numeric(0) else seq(first, last)
  }
  if (edge$uses_age) {
    w <- crossings(age0)
    x <- c(x, w - age0); age <- c(age, w); cal <- c(cal, cal0 + w - age0)
    dur <- c(dur, dur0 + w - age0)
  }
  if (edge$uses_caltime) {
    w <- crossings(cal0)
    x <- c(x, w - cal0); age <- c(age, age0 + w - cal0); cal <- c(cal, w)
    dur <- c(dur, dur0 + w - cal0)
  }
  if (edge$uses_duration) {
    w <- crossings(dur0)
    x <- c(x, w - dur0); age <- c(age, age0 + w - dur0)
    cal <- c(cal, cal0 + w - dur0); dur <- c(dur, w)
  }

  values <- rate_at(age, cal, dur)
  assert_rate_values(values, length(x), label)

  hit <- is.infinite(values) & x <= H
  if (!any(hit)) return(NA_real_)
  min(x[hit])
}


#' Print a rate registry
#'
#' Shows the number of states and resolved transitions. For each listed
#' transition, it shows the arguments its rate uses. A rate with no
#' arguments is labelled `constant`.
#'
#' Whether an edge can be forced (`Inf` rate) is reported only when a probe
#' grid was supplied at construction. Without one every edge is
#' conservatively treated as possibly forced.
#'
#' @param x A `rate_registry`.
#' @param max_lines Maximum number of transitions to list.
#' @param ... Not used; present for compatibility with the [print()] generic.
#' @return `x`, invisibly. Called for the side effect of printing.
#' @exportS3Method print rate_registry
#' @noRd
print.rate_registry <- function(x, max_lines = 20L, ...) {
  cat("<rate_registry>\n")
  cat("  ", length(x$edges), " resolved transition(s) over ",
      x$codec$n, " states\n", sep = "")

  if (!isTRUE(x$probed)) {
    cat("  forced-transition probe not run: every edge is treated as",
        " possibly forced\n", sep = "")
  }

  shown <- x$edges[seq_len(min(max_lines, length(x$edges)))]
  for (edge in shown) {
    scales <- if (length(edge$argnames)) {
      paste(edge$argnames, collapse = ", ")
    } else {
      "constant"
    }
    forced <- if (isTRUE(x$probed) && isTRUE(edge$may_be_infinite)) {
      "  [forced transitions possible]"
    } else {
      ""
    }
    cat("    ", edge$label, "   [depends on: ", scales, "]", forced,
        "\n", sep = "")
  }
  if (length(x$edges) > max_lines) {
    cat("    ... and ", length(x$edges) - max_lines, " more\n", sep = "")
  }
  invisible(x)
}
