# =====================================================================
# class_model.R
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   Provides the user-facing helpers for describing a simulation as one
#   reusable model object. A model keeps the inputs together so users can
#   inspect, save and run the same specification with different seeds.
#
# PLACE IN THE PIPELINE
#   user-defined values --> state_space() ------------------+
#   user-written rates --> transition() --> transitions() -+
#   population, dates and optional inputs -----------------+
#       --> micsim_model() [this file] --> model object
#       --> validate() [optional] --> simulate() [simulate.R]
#
# WHAT YOU NEED TO KNOW
#   The user defines the state variables and their possible values, supplies
#   population data and dates, and writes functions giving transition rates.
#   Here the functions themselves are passed to transition().
#
#   Building a model prepares the inputs but does not simulate life courses.
#   validate() adds sampled rate checks and returns the model, so a typical
#   workflow reads:
#
#       micsim_model(...) |>
#         add_subsystem(fertility(...)) |>
#         validate() |> simulate(seed = 234) |> as_long()
#
# NO RESERVED LABELS
#   Every label in a model is the user's own. A state is absorbing because
#   transitions lead into it and none lead out, not because it is called
#   "dead"; who can give birth, which transition is a school enrolment and
#   so on are told to the subsystem that needs to know. See selectors.R.
#
# GUARANTEES
#   * A model that has been constructed is structurally valid: its state
#     space, transitions, rates and newborn rule have all been checked
#     against each other. Nothing that could be caught here is deferred to
#     run time.
#   * Building a model runs no simulation and consumes no random numbers.
# =====================================================================


#' Describe a state space
#'
#' @description
#' Names the model's state variables and lists the values each may take.
#' Supply this object to [micsim_model()].
#' 
#' @details
#' The state space contains all combinations of the supplied values. For
#' example, two sexes, three nationalities and two parity values give twelve
#' states. A label such as `"f/DE/1"` follows the order of the arguments.
#'
#' Values must be distinct, both within and across variables. Avoid `"/"`
#' inside a value because that character separates label components.
#'
#' @param ... Named character vectors, one per state variable.
#'
#' @return An object of class `micsim_state_space`: the named list of 
#'   state-variable values.
#'
#' @seealso [micsim_model()], [MicSim-concepts]
#'
#' @examples
#' state_space(sex = c("m", "f"), nat = c("DE", "AT"), fert = c("0", "1"))
#'
#' @export
state_space <- function(...) {
  vars <- list(...)

  if (length(vars) == 0L || is.null(names(vars)) || any(names(vars) == "")) {
    micsim_abort(
      "micsim_model_error",
      paste0(
        "state_space() takes one NAMED character vector per state\n",
        "  variable, e.g. state_space(sex = c(\"m\", \"f\"))."
      )
    )
  }

  vars <- lapply(vars, as.character)
  values <- unlist(vars, use.names = FALSE)
  if (anyNA(values) || any(!nzchar(values)) ||
        any(grepl(STATE_SEP, values, fixed = TRUE))) {
    micsim_abort(
      "micsim_model_error",
      paste0(
        "State values must be non-empty strings without \"", STATE_SEP, "\",\n",
        "  which separates the components of a state label."
      )
    )
  }
  if (any(lengths(vars) == 0L) || any(vapply(vars, anyDuplicated, 1L) > 0L)) {
    micsim_abort(
      "micsim_model_error",
      "Every state variable needs at least one value, each listed once."
    )
  }
  if (anyDuplicated(values) > 0L) {
    micsim_abort(
      "micsim_model_error",
      paste0(
        "State values must be unique across variables. Repeated: ",
        toString(unique(values[duplicated(values)])), ".\n",
        "  Transition patterns match on values rather than on positions,\n",
        "  so a value used twice could not be told apart."
      )
    )
  }

  structure(vars, class = "micsim_state_space")
}


#' Describe one allowed transition
#'
#' @description
#' Describes an allowed move from one state to another and supplies the
#' function giving its rate.
#'
#' @details
#' Three ways to say which states:
#'
#' \describe{
#'   \item{`transition("female/0", "female/1", rate)`}{full state labels.
#'     `from` may be a vector of labels.}
#'   \item{`transition(NA, "deceased", rate)`}{`NA` means every transient
#'     state, which is how mortality is normally written.}
#'   \item{`transition(move("edu", "none", "primary"), rate)`}{one state
#'     variable changes and everything else stays, for every state (or
#'     those satisfying the edge's `where`). See [move()].}
#' }
#'
#' A destination that is not a transient state is an absorbing state: the
#' model creates it, under whatever name you give it.
#'
#' @param from Origin state label(s), `NA` for "every transient state except
#'   to", or a [move()] selector.
#' @param to Destination state label. Omitted when `from` is a [move()],
#'   which names its own destinations.
#' @param rate A function that returns transition rates, supplied as a
#'   function object rather than its name. It may accept any subset of
#'   `age`, `calTime` and `duration`, each measured in fractional years.
#'   It must work with vector inputs, returning one numeric rate for each
#'   corresponding set of input values. For a constant rate, you may
#'   instead supply a function with no arguments that returns a single
#'   numeric value.
#'
#' @return A transition, for passing to [transitions()].
#'
#' @seealso [transitions()], [move()], [micsim_model()]
#'
#' @examples
#' states <- state_space(
#'   sex = c("female", "male"),
#'   country = c("DE", "AT"),
#'   parity = c("0", "1", "2")
#' )
#' states
#'
#' # Full state labels follow the order of the variables in `states`.
#' 
#' # A constant rate is used for all moves between countries.
#' # This move from Germany to Austria applies to women with no children.
#' move_to_at <- function() 0.02
#' transition(from = "female/DE/0", to = "female/AT/0", rate = move_to_at)
#'
#' # A move from Austria to Germany applies to everyone living in Austria.
#' # Only the country changes. All other state variables remain constant.
#' move_to_de <- function() 0.03
#' transition(move("country", from = "AT", to = "DE"), rate = move_to_de)
#'
#' # Two origin states share the same emigration rate, and emigration is
#' # an absorbing state, as "abroad" is not in the state space.
#' emigration_rate <- function() 0.01
#' transition(from = c("female/DE/0", "male/DE/0"),
#'            to = "abroad", rate = emigration_rate)
#'
#' # First and second births change parity for women in either country.
#' first_birth <- function(age) ifelse(age >= 18 & age < 45, 0.10, 0)
#' second_birth <- function(age) ifelse(age >= 18 & age < 45, 0.08, 0)
#' transition(move("parity", from = "0", to = "1", where = sex == "female"),
#'            rate = first_birth)
#' transition(move("parity", from = "1", to = "2", where = sex == "female"),
#'            rate = second_birth)
#'
#' # Mortality applies in every transient state.
#' mortality <- function(age) 0.00003 * exp(0.1 * age)
#' transition(from = NA, to = "deceased", rate = mortality)
#'
#' @export
transition <- function(from, to, rate) {
  # The name the rate was given, when it was given by name, for reports
  # and plots ("mortRates(age, calTime)"). Read before anything is
  # evaluated.
  rate_expr <- if (!missing(rate)) substitute(rate) else if (!missing(to)) substitute(to)
  rate_name <- if (is.name(rate_expr)) as.character(rate_expr) else NA_character_
  if (inherits(from, "micsim_moves_of")) {
    micsim_abort(
      "micsim_model_error",
      paste0(
        "moves_of() selects among transitions that are already declared, so\n",
        "  it cannot declare one. Use move(var, from, to) instead."
      )
    )
  }
  if (inherits(from, "micsim_move")) {
    if (missing(rate)) {
      if (missing(to)) {
        micsim_abort("micsim_model_error", "A transition needs a `rate`.")
      }
      rate <- to                       # transition(move(...), rate)
    } else if (!missing(to)) {
      micsim_abort(
        "micsim_model_error",
        "A move() names its own destinations; give `rate` but not `to`."
      )
    }
    return(structure(list(move = from, rate = rate, rate_name = rate_name),
                     class = "micsim_transition"))
  }
  if (missing(to) || missing(rate)) {
    micsim_abort("micsim_model_error",
                 "A transition needs `from`, `to` and `rate`.")
  }
  structure(list(from = from, to = to, rate = rate, rate_name = rate_name),
            class = "micsim_transition")
}


#' Collect the transitions of a model
#'
#' @description
#' Collects the individual [transition()] declarations into one object for
#' [micsim_model()]. Each declaration describes an origin, destination and
#' rate function; the collection describes which moves the model allows.
#'
#' @examples
#' # Set the rates for the allowed transitions 
#' first_birth <- function(age) ifelse(age >= 18 & age < 45, 0.10, 0)
#' second_birth <- function(age) ifelse(age >= 18 & age < 45, 0.08, 0)
#' mortality <- function(age) 0.00003 * exp(0.1 * age)
#'
#' # Collect all the allowed transitions in one object
#' trans <- transitions(
#'   transition(move("parity", from = "0", to = "1", where = sex == "female"),
#'              rate = first_birth),
#'   transition(move("parity", from = "1", to = "2", where = sex == "female"),
#'              rate = second_birth),
#'   transition(from = NA, to = "deceased", rate = mortality)
#' )
#' trans
#'
#' @details
#' Order matters only for readability. Each destination receives its random
#' numbers in the order of the state space, so declaring the same transitions
#' in a different order gives the same population under the same seed.
#'
#' @param ... Objects from [transition()].
#'
#' @return An object of class `micsim_transitions`.
#'
#' @seealso [transition()], [micsim_model()]
#'
#' @export
transitions <- function(...) {
  trs <- list(...)
  if (length(trs) == 0L) {
    micsim_abort("micsim_model_error",
                 "A model needs at least one transition.")
  }
  if (!all(vapply(trs, inherits, logical(1L), "micsim_transition"))) {
    micsim_abort("micsim_model_error",
                 "transitions() takes objects built with transition().")
  }
  structure(unname(trs), class = "micsim_transitions")
}


#' Define a simulation model
#'
#' @description
#' Builds the object that says everything about a model: its states, its
#' transitions, who starts in it and how long it runs. Mechanisms that
#' ordinary transitions cannot express, such as fertility or a school
#' calendar, are added afterwards with [add_subsystem()].
#'
#' @details
#' Building a model enumerates the state space, resolves transitions and
#' checks population formats, dates and identifiers. It also checks rate
#' argument names and, where short test evaluations succeed, the number
#' of values returned. Subsystem-specific checks occur in [add_subsystem()].
#' These checks cannot establish that all rates are suitable for every
#' input encountered during a run. [validate()] provides an additional
#' screen for missing and negative rate values on a grid.
#'
#' **Absorbing states** are the destinations of transitions that are not
#' transient states. Call them anything: `"dead"`, `"deceased"`, `"tot"`
#' and `"abroad"` are all fine, and a model may have any number of them,
#' including none (then `max_age` alone ends each life course). If
#' something must tell a death from another kind of exit, say so in
#' `exits`.
#'
#' Nothing here consumes random numbers. A model can be built, printed,
#' saved with [saveRDS()] and simulated repeatedly, and every run from the
#' same seed gives the same result. A model saved under one MicSim version
#' and simulated under another gives a warning; rebuild it with micsim_model()
#' if the run fails or looks wrong.
#'
#' @param state_space User-defined state variables, from [state_space()].
#' @param transitions Allowed transitions, from [transitions()].
#' @param population Data frame describing the starting population, with
#'   columns `ID`, `birthDate` and `initState`, in that order. Use unique
#'   whole-number IDs within R's integer range, birth dates in `yyyymmdd`
#'   format (numeric, character, `Date`), and full transient-state labels.
#' @param horizon Two numeric dates in `yyyymmdd` format, `c(start, end)`, 
#'   with start before end.
#' @param max_age Maximum age for scheduling transitions, in years; default
#'   `100`.
#' @param immigrants Data frame with columns `ID`, `immigrDate`,
#'   `birthDate` and `immigrInitState`, or `NULL`.
#' @param exits Optional named character vector saying what kind of exit
#'   each absorbing state is: `"death"`, `"emigration"` or `"other"`, e.g.
#'   `c(deceased = "death", abroad = "emigration")`. Absorbing states left
#'   out are of unspecified kind. Nothing in core MicSim depends on it, but
#'   subsystems that must know whether someone died do.
#' @param ... Not used.
#'
#' @return A `micsim_model` containing the supplied inputs and internal
#'   `codec` (state mapping), `registry` (transition-rate lookup) and
#'   `micsim_version` (the MicSim version that built the model).
#'
#' @seealso [state_space()], [transitions()], [add_subsystem()],
#'   [validate()], [simulate.micsim_model()], [MicSim-concepts]
#'
#' @examples
#' mortality   <- function(age) 0.00003 * exp(0.1 * age)
#' first_birth <- function(age) ifelse(age >= 18 & age < 45, 0.10, 0)
#'
#' # The birth is declared once and reused by the fertility subsystem
#' birth <- move("parity", from = "0", to = "1", where = sex == "female")
#'
#' trans <- transitions(
#'   transition(birth, rate = first_birth),
#'   transition(from = NA, to = "deceased", rate = mortality)
#' )
#'
#' model <- micsim_model(
#'   state_space = state_space(sex = c("female", "male"), parity = c("0", "1")),
#'   transitions = trans,
#'   population = data.frame(
#'     ID        = 1:20,
#'     birthDate = rep("19850101", 20),
#'     initState = rep("female/0", 20)
#'   ),
#'   horizon = c(20140101, 20241231),
#'   max_age = 100,
#'   exits   = c(deceased = "death")
#' ) |>
#'   add_subsystem(fertility(
#'     births  = birth,
#'     newborn = newborn_rule(sex = c(male = 0.515, female = 0.485), parity = "0")
#'   ))
#' model
#'
#' @export
micsim_model <- function(state_space,
                         transitions,
                         population,
                         horizon,
                         max_age    = 100,
                         immigrants = NULL,
                         exits      = NULL,
                         ...) {

  refuse_retired_arguments(list(...))

  if (!inherits(state_space, "micsim_state_space")) {
    micsim_abort("micsim_model_error",
                 "`state_space` must be built with state_space().")
  }
  if (!inherits(transitions, "micsim_transitions")) {
    micsim_abort("micsim_model_error",
                 "`transitions` must be built with transitions().")
  }

  absorbing <- absorbing_states_of(transitions, state_space)
  codec <- new_state_codec(unclass(state_space), absorbing = absorbing)

  # Building the registry resolves every rate, derives its time-scale
  # dependence and checks its shape, so a malformed transition is an error
  # here rather than three minutes into a run.
  registry <- new_rate_registry(unclass(transitions), codec,
                                probe = forced_probe_grid(max_age, horizon))

  model <- structure(
    list(
      state_space = state_space,
      transitions = transitions,
      population  = standardise_people(population),
      horizon     = unname(horizon),
      max_age     = max_age,
      immigrants  = standardise_people(immigrants),
      exits       = resolve_exits(exits, absorbing),
      subsystems  = list(),
      codec       = codec,
      registry    = registry,
      # The MicSim version that built the model. A model stores internal
      # objects (state codec, rate registry, compiled subsystems) whose
      # structure may change between versions; simulate() compares this
      # stamp with the installed version.
      micsim_version = as.character(getNamespaceVersion("MicSim"))
    ),
    class = "micsim_model"
  )

  # Building the engine applies its input checks without running it.
  invisible(micsim_engine_for(model, verbose = FALSE))

  model
}


#' Identify absorbing states from transition destinations
#'
#' @description
#' Identifies absorbing states as transition destinations that are not
#' among the model's transient states. No particular label, such as
#' `"dead"`, is required.
#'
#' @details
#' Each absorbing state must have a non-empty label that contains no
#' `"/"` and does not match an individual state-variable value. These
#' checks help catch incomplete or misspelled transient state labels.
#'
#' @param transitions A `micsim_transitions`.
#' @param state_space A `micsim_state_space`.
#' @return Character vector of absorbing labels, in the order they are
#'   first used.
#' @keywords internal
#' @noRd
absorbing_states_of <- function(transitions, state_space) {
  ss <- unclass(state_space)
  grid <- expand.grid(rev(ss), stringsAsFactors = FALSE)[, rev(seq_along(ss)),
                                                        drop = FALSE]
  transient <- do.call(paste, c(unname(as.list(grid)), sep = STATE_SEP))
  values <- unlist(ss, use.names = FALSE)

  to <- unlist(lapply(unclass(transitions), function(tr) {
    if (is.null(tr$move)) as.character(tr$to) else NULL
  }), use.names = FALSE)
  absorbing <- unique(setdiff(to, transient))

  for (label in absorbing) {
    if (is.na(label) || !nzchar(label)) {
      micsim_abort("micsim_model_error",
                   "A transition's destination must be a non-empty state label.")
    }
    if (grepl(STATE_SEP, label, fixed = TRUE)) {
      micsim_abort(
        "micsim_model_error",
        paste0(
          "The destination `", label, "` is not a state of the model. A\n",
          "  composite label must combine one value of every state variable\n",
          "  (", toString(names(ss)), "), in that order."
        )
      )
    }
    if (label %in% values) {
      micsim_abort(
        "micsim_model_error",
        paste0(
          "The destination `", label, "` is a value of a state variable, not\n",
          "  a state. Write the full state label, or declare the change with\n",
          "  move(); an absorbing state needs a name of its own."
        )
      )
    }
  }
  absorbing
}


#' Check and assign exit types for absorbing states
#'
#' @description
#' Checks that each entry in `exits` names an absorbing state and supplies
#' a supported exit type. Names must be non-empty and unique. Absorbing
#' states without a declared exit type are assigned `NA`.
#'
#' @param exits Named character vector whose names are absorbing state
#'   labels and whose values are exit types, for example
#'   `c(dead = "death")`. Use `NULL` to leave all exit types unspecified.
#' @param absorbing Character vector of the model's absorbing state labels.
#' @return A character vector over `absorbing`, named by it, holding each
#'   state's kind or `NA`.
#' @keywords internal
#' @noRd
resolve_exits <- function(exits, absorbing) {
  out <- stats::setNames(rep(NA_character_, length(absorbing)), absorbing)
  if (is.null(exits)) return(out)
  if (!is.character(exits) || is.null(names(exits)) ||
        any(!nzchar(names(exits))) || anyDuplicated(names(exits)) > 0L) {
    micsim_abort(
      "micsim_model_error",
      "`exits` must be a named character vector, e.g. c(dead = \"death\")."
    )
  }
  unknown <- setdiff(names(exits), absorbing)
  if (length(unknown) > 0L) {
    micsim_abort(
      "micsim_model_error",
      paste0(
        "`exits` names state(s) that are not absorbing states of the model: ",
        toString(unknown), ".\n",
        "  Absorbing states: ",
        if (length(absorbing)) toString(absorbing) else "none", "."
      )
    )
  }
  bad <- setdiff(exits, EXIT_KINDS)
  if (length(bad) > 0L) {
    micsim_abort(
      "micsim_model_error",
      paste0("Unknown exit kind(s) ", toString(bad), ". Use one of ",
             toString(EXIT_KINDS), ".")
    )
  }
  out[names(exits)] <- unname(exits)
  out
}


#' Reject unsupported arguments to micsim_model()
#'
#' @description
#' Stops with an error if any arguments were supplied through `...`.
#' For recognised arguments, the error includes guidance on how to
#' update the model call.
#'
#' @param dots List of arguments supplied through `...` in [micsim_model()].
#' @return Invisibly returns `TRUE` when `dots` is empty. Otherwise,
#'   raises a `micsim_model_error`.
#' @keywords internal
#' @noRd
refuse_retired_arguments <- function(dots) {
  if (length(dots) == 0L) return(invisible(TRUE))
  nms <- names(dots)
  if (is.null(nms)) nms <- rep("", length(dots))
  advice <- c(
    absorbing        = "absorbing states are now the destinations that are not transient states; there is nothing to declare (use `exits` to say what kind of exit each is)",
    fertility        = "use add_subsystem(fertility(births = ..., newborn = ...))",
    newborns         = "the newborn rule is now the `newborn` argument of fertility()",
    school_enrolment = "use add_subsystem(school_enrolment(move(...), month = ...))",
    hooks            = "subsystems are added with add_subsystem()"
  )
  known <- intersect(nms, names(advice))
  msg <- if (length(known) > 0L) {
    paste0(
      "micsim_model() no longer takes ", toString(paste0("`", known, "`")),
      ":\n", paste0("  * ", known, ": ", advice[known], collapse = "\n")
    )
  } else {
    paste0("Unused argument(s) to micsim_model(): ",
           toString(ifelse(nzchar(nms), nms, "<unnamed>")), ".")
  }
  micsim_abort("micsim_model_error", msg)
}


#' Build the engine a model describes
#'
#' The one place that translates a `micsim_model` into engine arguments,
#' so that [micsim_model()] and [simulate.micsim_model()] cannot drift
#' apart.
#'
#' @param model A `micsim_model`.
#' @param verbose Passed to `new_engine()`.
#' @return A `micsim_engine`.
#' @keywords internal
#' @noRd
micsim_engine_for <- function(model, verbose) {
  new_engine(
    codec      = model$codec,
    registry   = model$registry,
    init_pop   = model$population,
    horizon    = model$horizon,
    max_age    = model$max_age,
    subsystems = unname(model$subsystems),
    immigrants = model$immigrants,
    verbose    = verbose
  )
}



#' Checks rate values before simulating
#' Replace the population of a model
#'
#' @description
#' Returns a model with a different starting population, immigrant
#' population, or both. All other settings are retained, and the new
#' populations are checked exactly as [micsim_model()] checks them.
#'
#' This allows a model developed on a sample to be applied to the full
#' population without repeating the model specification.
#'
#' @details
#' Only `population` and `immigrants` can be changed with this method.
#' To change the state space, transitions, simulation period, maximum age,
#' exit types or subsystems, rebuild the model with [micsim_model()] and,
#' where needed, [add_subsystem()].
#'
#' Replace populations through `update()` so that the model is checked.
#' Direct assignment, such as `model$population <- ...`, bypasses these
#' checks and is not supported.
#'
#' @param object A `micsim_model`.
#' @param population A new starting population, in the format
#'   [micsim_model()] takes. Leave it out to keep the current one.
#' @param immigrants A new immigrant table, or `NULL` for none. Leave it
#'   out to keep the current one.
#' @param ... Must be empty.
#' @return A `micsim_model`.
#' @seealso [micsim_model()], [add_subsystem()]
#' @examples
#' mort <- function(age, calTime) 0.0001 * exp(0.085 * age)
#' pop <- data.frame(ID = 1:200, birthDate = "19700101",
#'                   initState = rep(c("m", "f"), 100))
#' sample_model <- micsim_model(state_space(sex = c("m", "f")),
#'                              transitions(transition(NA, "dead", mort)),
#'                              pop[1:20, ], c(20000101, 20101231))
#' full_model <- update(sample_model, population = pop)
#' nrow(full_model$population)
#' @export
update.micsim_model <- function(object, population, immigrants, ...) {
  extra <- list(...)
  if (length(extra) > 0L) {
    micsim_abort("micsim_model_error", paste0(
      "update() can replace `population` and `immigrants` only.\n",
      "  Not accepted: ", toString(if (is.null(names(extra))) "(unnamed)" else names(extra)), ".\n",
      "  To change anything else, build the model again with micsim_model()."
    ))
  }
  if (!missing(population)) object$population <- standardise_people(population)
  if (!missing(immigrants)) object$immigrants <- standardise_people(immigrants)
  # the same checks micsim_model() and add_subsystem() make
  invisible(micsim_engine_for(object, verbose = FALSE))
  object
}


#' Check whether any subsystem can add individuals
#' 
#' @param model A `micsim_model`.
#' @return `TRUE` or `FALSE`.
#' @keywords internal
#' @noRd
model_spawns <- function(model) {
  any(vapply(model$subsystems, function(s) isTRUE(s$spawns), logical(1L)))
}


#' Check an object thoroughly
#'
#' @description
#' For a [micsim_model()], runs every check and reports everything it finds at once.
#' Errors: a rate that raises an error, returns something other than one number
#' per input, `NA` or a negative value on a grid of whole ages, calendar years and
#' durations; a newborn table that has no state for a mother. Warnings: a rate that
#' is zero everywhere or exceeds 10 per person-year; a transition nobody can reach;
#' a state nobody can leave. `Inf` is allowed for forced transitions.
#'
#' @details
#' During construction the model already checked that each rate is a function
#' with admissible arguments and that it returns one value per input.
#'
#' Validation checks the values at the grid points. For example, a rate
#' that passes at whole ages can still fail at a fractional age. Assessing
#' whether the model is suitable for the research question is a separate step.
#'
#' Because validation returns the model, you can pass it directly to
#' `model |> validate() |> simulate(seed = 234) `
#'
#' @param x An object to check, here `micsim_model`.
#' @param strict If `TRUE`, warnings are treated as errors.
#' @param ... Passed to methods.
#'
#' @return `x`, invisibly, if no error is found. Otherwise an error of class
#' `micsim_invalid_model` listing every problem. Warnings are given as one
#' `micsim_model_warning`. Both carry `$problems`, a data frame with one row per finding.
#'
#' @seealso [micsim_model()], [MicSim-concepts]
#'
#' @export
validate <- function(x, ...) {
  UseMethod("validate")
}


#' @rdname validate
#' @export
validate.micsim_model <- function(x, strict = FALSE, ...) {
  # Every check runs and every finding is reported at once; see
  # R/validate.R.
  validate_model(x, strict = strict)
}


#' Construct the rate-validation grid 
#'
#' @description
#' Builds integer-age, calendar-year and duration sequences for `validate()`.
#' This grid is a screening tool, not an exact set of reachable inputs.
#'
#' @details
#' The lower age is the youngest entry age rounded down, bounded at zero,
#' or zero when fertility is enabled. The upper age is at least that lower
#' bound and otherwise `ceiling(max_age) - 1`. Calendar years run from the
#' starting year to the ending year of the simulation horizon, inclusive.
#' Durations run from zero to the upper age bound minus the lower age
#' bound. 
#'
#' @param model A `micsim_model`.
#' @return A grid from `micsim_rate_probe()`.
#' @keywords internal
#' @noRd
reachable_grid <- function(model) {

  horizon_days <- getInDays(model$horizon)

  entry_ages <- (horizon_days[1L] -
                   getInDays(model$population$birthDate)) / DAYS_PER_YEAR
  if (!is.null(model$immigrants)) {
    entry_ages <- c(
      entry_ages,
      (getInDays(model$immigrants$immigrDate) -
         getInDays(model$immigrants$birthDate)) / DAYS_PER_YEAR
    )
  }
  youngest <- if (model_spawns(model)) {
    0                                  # newborns start at age zero
  } else {
    min(entry_ages, na.rm = TRUE)
  }

  lowest  <- max(0L, as.integer(floor(youngest)))
  highest <- max(lowest, as.integer(ceiling(model$max_age)) - 1L)
  years   <- getYear(horizon_days)

  # Durations run up to the oldest age, not to the spread of ages: every
  # duration clock of the starting population starts at birth, so a
  # 94-year-old has a duration of 94 the moment the run begins.
  micsim_rate_probe(
    ages      = seq.int(lowest, highest),
    cal_years = seq.int(years[1L], years[2L]),
    durations = seq.int(0L, highest)
  )
}


#' Print a simulation model
#'
#' Provides an overview of the model, including the number of states and
#' transitions, the starting population, the simulation period, and any
#' optional mechanisms enabled.
#'
#' @param x A `micsim_model`.
#' @param ... Not used; present for compatibility with the [print()] generic.
#' @return `x`, invisibly. Called for the side effect of printing.
#' @export
print.micsim_model <- function(x, ...) {
  cat("<micsim_model>\n")
  cat("  state space : ",
      paste(names(x$state_space), collapse = " x "), " = ",
      x$codec$K, " transient state(s)\n", sep = "")
  absorbing <- x$codec$absorbing
  if (length(absorbing) > 0L) {
    kinds <- x$exits[absorbing]
    shown <- ifelse(is.na(kinds), absorbing, paste0(absorbing, " (", kinds, ")"))
    cat("  absorbing   : ", toString(shown), "\n", sep = "")
  } else {
    cat("  absorbing   : none (life courses end at max age)\n")
  }
  cat("  transitions : ", length(x$registry$edges),
      " (from ", length(x$transitions), " declaration(s))\n", sep = "")
  cat("  population  : ", nrow(x$population), " individual(s)\n", sep = "")
  if (!is.null(x$immigrants)) {
    cat("  immigrants  : ", nrow(x$immigrants), " arrival(s)\n", sep = "")
  }
  cat("  horizon     : ", x$horizon[1L], " to ", x$horizon[2L],
      ", max age ", x$max_age, "\n", sep = "")
  if (length(x$subsystems) == 0L) {
    cat("  subsystems  : none\n")
  } else {
    cat("  subsystems  :\n")
    for (sub in x$subsystems) {
      cat("    ", sub$name, ": ", sub$describe, "\n", sep = "")
    }
  }
  invisible(x)
}
