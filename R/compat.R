# =====================================================================
# compat.R
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   Prepares inputs from the classic micSim() interface for the new
#   engine
#
# PLACE IN THE PIPELINE
#   the user's micSim() call
#        |
#        v
#   compat.R  -- reconstructs --> codec.R, registry.R,
#             --                  class_newborn_rule.R
#        |
#        v
#   engine.R --> the `pop` data frame
#
# THE THREE TRANSLATIONS IT PERFORMS
#   1. A TRANSITION MATRIX becomes a state codec and a rate registry. The
#      matrix's row and column names are the state space; its non-"0" cells
#      are the transitions; the strings in those cells are the NAMES of
#      rate functions, resolved to the functions themselves once, here,
#      instead of on every event.
#
#   2. FERTILITY PATTERNS become explicit transition labels. The `fertTr`
#      argument accepts partial patterns such as "0->1"; the engine wants
#      the full list.
#
#   3. THREE LOOSE NEWBORN ARGUMENTS become one validated object.
#
#   4. monthSchoolEnrol becomes a school_enrolment() subsystem, and the
#      fertility arguments a fertility() subsystem, compiled exactly as
#      add_subsystem() compiles them for the object API.
#
# THE ONE PLACE LABELS STILL MEAN SOMETHING
#   The classic interface uses specific labels for fertility and school
#   enrolment. Fertility patterns match states containing "f". School
#   enrolment uses transitions from "no" to "low" in one state variable.
#   buildTransitionMatrix() places "dead" and "rest" first among absorbing
#   states when those labels are present. Those conventions are kept for this 
#   legacy interface. The object API has no reserved labels at all.
#
# GUARANTEES
#   * The state space follows the matrix's rows and columns. Its order
#     decides which of an individual's random numbers each destination
#     receives (rng_streams.R), so changing it changes every simulated
#     population without changing the model.
#   * micSim() takes no seed: it draws the seed of the run from the
#     caller's random stream (one draw), so `set.seed(234); micSim(...)` is
#     reproducible.
# =====================================================================


#' Run a continuous-time microsimulation
#'
#' @description
#' Simulates individual life courses as sequences of changes between states.
#' This is the classic interface, which specifies transitions and their
#' rate functions in a transition matrix.
#'
#' @details
#' This is the classic interface based on a transition matrix whose
#' cells name rate functions. [micsim_model()] provides an alternative
#' interface that keeps the model inputs together and accepts functions
#' directly. Both use the same simulation engine.
#'
#' **The state space** contains all possible combinations of the model's state
#' variables, with each state written as its values joined by forward
#' slashes: `"f/DE/1"` for a German woman with one child. The absorbing
#' states are the columns of `transitionMatrix` after the transient ones;
#' any label works.
#'
#' **Labels with a meaning.** This interface keeps three conventions of
#' MicSim 3.0.0, because it has no other way to be told: only states
#' containing `"f"` can give birth through `fertTr`; with
#' `monthSchoolEnrol`, an allowed transition counts as enrolment
#'  if it changes one state variable from `"no"` to `"low"` and leaves
#'  the other values as they were; and [buildTransitionMatrix()] puts
#' `"dead"` and `"rest"` first among the absorbing states. 
#'
#' **Rates** are functions of `age`, `calTime` and `duration`, in any
#' subset, all in fractional years. They are transition intensities, not
#' probabilities: a rate may exceed one. They must be vectorised, because
#' the simulation integrates them numerically.
#'
#' **Reproducibility.** `micSim()` takes no seed. It draws the seed of
#' the run from the random stream the caller has set up (one draw), so
#' `set.seed(234)` immediately before the call makes a run reproducible.
#' Every individual then draws from their own random-number stream, as
#' described in [simulate.micsim_model()].
#'
#' @section Reproducibility across versions:
#' The return value keeps the same columns, row semantics and ordering.
#' Some values differ from earlier releases, most importantly because dates
#' are now computed on the exact civil calendar. Random numbers are now drawn
#' per individual, so the same seed gives a different, equally valid
#' population than in MicSim 3.0.0.
#'
#' @param initPop Data frame describing the starting population, with
#'   columns `ID`, `birthDate` and `initState`, in that order. Birth dates
#'   are `yyyymmdd`.
#' @param immigrPop Data frame describing individuals who enter the
#'   population after the simulation starts, with columns `ID`,
#'   `immigrDate`, `birthDate` and `immigrInitState`. `NULL` when the
#'   model has no immigration.
#' @param transitionMatrix Character matrix as produced by
#'   [buildTransitionMatrix()]: rows are transient states, columns are all
#'   states, and a cell holds `"0"` or the name of a rate function.
#' @param absStates The absorbing states. Accepted for backward
#'   compatibility but no longer needed: they are read from the matrix's
#'   column names, which removes the risk of the two disagreeing.
#' @param fixInitStates Integer positions of the state variables a newborn
#'   inherits from its mother.
#' @param varInitStates Matrix of admissible newborn states, one row per
#'   candidate and one column per state variable.
#' @param initStatesProb Weights for the rows of `varInitStates`. They must
#'   sum to one within each combination of inherited attributes.
#' @param maxAge Sharp maximum age, in years.
#' @param simHorizon Length-2 vector `c(startDate, endDate)`, `yyyymmdd`.
#' @param fertTr Transitions that create a newborn, as
#'   `"origin->destination"` patterns. Empty when the model has no
#'   fertility.
#' @param monthSchoolEnrol Month, 1 to 12, on whose first day
#'   elementary-school enrolment takes place. Empty when the model does
#'   not track school enrolment.
#' @param verbose Report progress through [message()], which can be
#'   suppressed or captured. Pass `FALSE` for a silent run.
#'
#' @return A data frame with one row per transition, plus one row per
#'   individual who never experienced one, and the columns `ID`,
#'   `birthDate`, `initState`, `From`, `To`, `transitionTime` and
#'   `transitionAge`, with `motherID` added when the model has fertility.
#'
#' @references
#' Zinn, S. (2014). The MicSim Package of R: An Entry-Level Toolkit for
#' Continuous-Time Microsimulation. *International Journal of
#' Microsimulation*, 7(3), 3-32. \doi{10.34196/IJM.00105}
#'
#' @seealso [buildTransitionMatrix()], [MicSim-concepts]
#'
#' @examples
#' # A mortality-only model: five individuals, one absorbing state.
#' stateSpace <- c("m", "f")
#' attr(stateSpace, "name") <- "sex"
#'
#' mortRates <- function(age, calTime) {
#'   0.00003 * exp(ifelse(calTime <= 2020, 0.1, 0.097) * age)
#' }
#'
#' initPop <- data.frame(
#'   ID        = 1:5,
#'   birthDate = c("19301231", "19990403", "19561015", "19911111", "19650101"),
#'   initState = c("f", "m", "f", "m", "m"),
#'   stringsAsFactors = FALSE
#' )
#'
#' transitionMatrix <- buildTransitionMatrix(
#'   allTransitions = NULL,
#'   absTransitions = cbind("dead", "mortRates"),
#'   stateSpace     = stateSpace
#' )
#'
#' set.seed(234)
#' pop <- micSim(
#'   initPop          = initPop,
#'   transitionMatrix = transitionMatrix,
#'   absStates        = "dead",
#'   maxAge           = 120,
#'   simHorizon       = c(20000101, 21001231),
#'   verbose          = FALSE
#' )
#' head(pop)
#'
#' @section Deprecated:
#' `micSim()` is deprecated as of MicSim 4.0.0 and will be removed in a
#' future release. It still works, and shows a message once per R session
#' naming its replacement: build the model with [micsim_model()] (and
#' [add_subsystem()]) and run it with [simulate.micsim_model()]. Set
#' `options(MicSim.show_deprecation = FALSE)` to turn the message off.
#'
#' @keywords internal
#' @export
micSim <- function(initPop,
                   immigrPop        = NULL,
                   transitionMatrix,
                   absStates        = NULL,
                   fixInitStates    = c(),
                   varInitStates    = c(),
                   initStatesProb   = c(),
                   maxAge           = 99,
                   simHorizon,
                   fertTr           = c(),
                   monthSchoolEnrol = c(),
                   verbose          = TRUE) {
  deprecate_message("micSim")
  # Rate function names are resolved in the environment micSim() was
  # called from, which is where the user's script defines them.
  legacy_micsim(
    initPop          = initPop,
    immigrPop        = immigrPop,
    transitionMatrix = transitionMatrix,
    fixInitStates    = fixInitStates,
    varInitStates    = varInitStates,
    initStatesProb   = initStatesProb,
    maxAge           = maxAge,
    simHorizon       = simHorizon,
    fertTr           = fertTr,
    monthSchoolEnrol = monthSchoolEnrol,
    verbose          = verbose,
    caller           = parent.frame()
  )
}


#' Run the classic interface without the deprecation message
#'
#' @description
#' The body of [micSim()]. `micSim()` shows the deprecation message and
#' calls this; the workers of [micSimParallel()] call it directly, so the
#' message is not repeated on every core.
#'
#' @param initPop,immigrPop,transitionMatrix,fixInitStates,varInitStates,initStatesProb,maxAge,simHorizon,fertTr,monthSchoolEnrol,verbose
#'   As in [micSim()].
#' @param caller Environment in which to start looking up the rate
#'   functions named in `transitionMatrix`, normally where the user's
#'   script called `micSim()`.
#' @return The `pop` data frame, as [micSim()] documents it.
#' @keywords internal
#' @noRd
legacy_micsim <- function(initPop, immigrPop, transitionMatrix,
                          fixInitStates, varInitStates, initStatesProb,
                          maxAge, simHorizon, fertTr, monthSchoolEnrol,
                          verbose, caller, seed = NULL) {

  # Read factor and Date columns by their values (see standardise_people()).
  initPop   <- standardise_people(initPop)
  immigrPop <- standardise_people(immigrPop)

  # Rate function names are resolved once here, in `caller`, rather than
  # on every event.
  model <- legacy_model(transitionMatrix, fertTr, fixInitStates,
                        varInitStates, initStatesProb, monthSchoolEnrol,
                        caller,
                        probe = forced_probe_grid(maxAge, unname(simHorizon)))

  engine <- new_engine(
    codec      = model$codec,
    registry   = model$registry,
    init_pop   = initPop,
    horizon    = unname(simHorizon),
    max_age    = maxAge,
    subsystems = model$subsystems,
    immigrants = immigrPop,
    verbose    = verbose
  )

  # micSim() has no seed, so the run draws its seed from the caller's
  # stream; micSimParallel() passes its `seeds` to every worker.
  pop <- engine$run(seed = seed)$pop

  # Save the model information needed by convertToLongFormat() as an
  # attribute so it can use the information directly from pop.
  attr(pop, "micsim_spec") <- list(
    absStates  = model$codec$absorbing,
    stateVars  = model$codec$vars,
    simHorizon = unname(simHorizon),
    immigrPop  = immigrPop,
    maxAge     = maxAge
  )
  pop
}


#' Prepare classic model inputs for the engine
#'
#' @description
#' Builds the state codec and rate registry from the transition matrix and
#' prepares the optional fertility and school enrolment subsystems. Only
#' transitions allowed by the matrix can be used by these subsystems.
#'
#' @details
#' The subsystems are the same ones the object API installs with
#' [add_subsystem()], compiled the same way, so the two interfaces cannot
#' simulate different things. Fertility is compiled before school
#' enrolment.
#'
#' @param transitionMatrix,fertTr,fixInitStates,varInitStates,initStatesProb,monthSchoolEnrol
#'   As in [micSim()].
#' @param caller Environment in which to start looking up rate functions.
#' @param probe Age, calendar-year and duration values from
#'   `forced_probe_grid()` used to look for rates that return `Inf`.
#'   With `NULL`, any rate may represent a forced transition.
#' @return A list with `codec`, `registry` and `subsystems`.
#' @keywords internal
#' @noRd
legacy_model <- function(transitionMatrix, fertTr, fixInitStates,
                         varInitStates, initStatesProb, monthSchoolEnrol,
                         caller, probe = NULL) {

  codec <- codec_from_matrix(transitionMatrix)
  registry <- registry_from_matrix(transitionMatrix, codec, caller, probe)
  subsystems <- list()

  if (length(monthSchoolEnrol) > 0L) {
    month <- check_enrolment_month(monthSchoolEnrol)
    enrolment <- legacy_enrolment_transitions(codec)
    declared <- intersect(enrolment,
                          vapply(registry$edges, `[[`, character(1L), "label"))
    # A state space with enrolment pairs that the matrix never uses has
    # nothing to snap; 3.0.0 accepted it silently, and so does this.
    if (length(declared) > 0L) {
      subsystems$school_enrolment <- compile_subsystem(
        school_enrolment(declared, month), codec, registry
      )
    }
  }

  if (length(fertTr) == 0L) {
    return(list(codec = codec, registry = registry,
                subsystems = unname(subsystems)))
  }

  if (length(varInitStates) == 0L || length(initStatesProb) == 0L) {
    micsim_abort(
      "micsim_compat_error",
      paste0(
        "The model has fertility transitions (`fertTr`), so it also needs\n",
        "  to say what state a newborn is in: `varInitStates` lists the\n",
        "  admissible newborn states and `initStatesProb` their weights."
      )
    )
  }

  births <- expand_fertility_transitions(fertTr, codec)
  declared <- vapply(registry$edges, `[[`, character(1L), "label")
  births <- intersect(births, declared)
  if (length(births) == 0L) {
    micsim_abort(
      "micsim_compat_error",
      paste0(
        "None of the `fertTr` patterns matches a transition the transition\n",
        "  matrix allows, so no birth could ever happen. Check that the\n",
        "  patterns name female states (containing \"",
        legacy_roles()$female, "\") and transitions the matrix has."
      )
    )
  }
  # Fertility first, so that the subsystems are in the order 4.0.0 always
  # installed them; see Details.
  subsystems <- c(
    list(fertility = compile_subsystem(
      fertility(births  = births,
                newborn = new_newborn_rule(states  = varInitStates,
                                           prob    = initStatesProb,
                                           inherit = fixInitStates)),
      codec, registry
    )),
    subsystems
  )

  list(codec = codec, registry = registry, subsystems = unname(subsystems))
}


#' The label conventions of the classic interface
#'
#' @description
#' MicSim 3.0.0 attached meaning to five labels. The classic interface
#' keeps them, because a transition matrix has no other way to say what
#' they say; the object API uses none of them.
#'
#' \describe{
#'   \item{`dead`, `emigrated`}{`"dead"` and `"rest"`: the absorbing
#'     states placed first, in that order, among the transition matrix's
#'     absorbing columns. Their order is the state space's order of
#'     absorbing states, which decides which random number goes to
#'     which destination.}
#'   \item{`female`}{`"f"` must appear as a state value on both sides of a
#'     transition for a `fertTr` pattern to match it.}
#'   \item{`no_education`, `low_education`}{`"no"` and `"low"`: with
#'     `monthSchoolEnrol`, a change from the first to the second in one
#'     state variable is a school enrolment.}
#' }
#'
#' @return A named list of labels.
#' @keywords internal
#' @noRd
legacy_roles <- function() {
  list(
    dead          = "dead",
    emigrated     = "rest",
    female        = "f",
    no_education  = "no",
    low_education = "low"
  )
}


#' Order absorbing states as MicSim 3.0.0 did
#'
#' @param labels Absorbing state labels, in the order first named.
#' @return The same labels, with `"dead"` and `"rest"` (where present)
#'   moved to the front, in that order.
#' @keywords internal
#' @noRd
legacy_absorbing_order <- function(labels) {
  roles <- legacy_roles()
  first <- intersect(c(roles$dead, roles$emigrated), labels)
  c(first, setdiff(labels, first))
}


#' Find possible school enrolment transitions
#'
#' @description
#' Finds pairs of transient states that differ only in one variable, which
#' changes from `"no"` to `"low"`. The caller checks which pairs are allowed
#' by the transition matrix.
#'
#' @details
#' Origin and destination labels are compared component by component,
#' matching each state variable's value as a whole.
#'
#' @param codec A `state_codec`.
#' @return Character vector of transition labels, `"from -> to"`.
#' @keywords internal
#' @noRd
legacy_enrolment_transitions <- function(codec) {
  roles <- legacy_roles()
  transient <- codec$transient
  parts <- strsplit(transient, "/", fixed = TRUE)
  labels <- character(0)
  for (i in seq_along(transient)) {
    origin <- parts[[i]]
    for (position in which(origin == roles$no_education)) {
      wanted <- origin
      wanted[position] <- roles$low_education
      j <- match(paste(wanted, collapse = "/"), transient)
      if (!is.na(j)) {
        labels <- c(labels, paste(transient[i], "->", transient[j]))
      }
    }
  }
  if (length(labels) == 0L) {
    abort_enrolment(paste0(
      "School enrolment was requested (`monthSchoolEnrol`), but no transition\n",
      "  in the state space goes from `", roles$no_education, "` to `",
      roles$low_education, "` in one state variable.\n",
      "  In the classic interface, elementary-school enrolment is modelled by\n",
      "  an education variable with those two values. The object API names\n",
      "  the transition instead: school_enrolment(move(var, from, to), month)."
    ))
  }
  labels
}


#' Reconstruct the state space from a transition matrix
#'
#' @description
#' Uses row names as transient states and column names without a matching
#' row as absorbing states. It splits the transient labels at `/` to recover
#' the values of each state variable.
#'
#' @details
#' Absorbing states are taken from the columns. The state variables have no
#' names here, because the matrix does not carry them; they are numbered. 
#' Names only appear in diagnostics, never in results.
#'
#' @param transitionMatrix As in [micSim()].
#' @return A `state_codec`.
#' @keywords internal
#' @noRd
codec_from_matrix <- function(transitionMatrix) {

  if (!is.matrix(transitionMatrix) ||
        is.null(rownames(transitionMatrix)) ||
        is.null(colnames(transitionMatrix))) {
    micsim_abort(
      "micsim_compat_error",
      paste0(
        "`transitionMatrix` must be a matrix with row and column names.\n",
        "  Build it with buildTransitionMatrix()."
      )
    )
  }

  transient <- rownames(transitionMatrix)
  absorbing <- setdiff(colnames(transitionMatrix), transient)
  if (length(absorbing) == 0L) {
    micsim_abort(
      "micsim_compat_error",
      paste0(
        "`transitionMatrix` has no absorbing states: every column is also\n",
        "  a row. A model needs at least one, usually \"dead\", or no life\n",
        "  course can ever end."
      )
    )
  }

  parts <- strsplit(transient, "/", fixed = TRUE)
  widths <- lengths(parts)
  if (length(unique(widths)) != 1L) {
    micsim_abort(
      "micsim_compat_error",
      paste0(
        "The transient states do not all have the same number of state\n",
        "  variables: found ", toString(sort(unique(widths))),
        " components.\n  Every state must name every variable."
      )
    )
  }

  n_vars <- widths[1L]
  state_space <- lapply(seq_len(n_vars), function(position) {
    unique(vapply(parts, `[`, character(1L), position))
  })
  names(state_space) <- paste0("var", seq_len(n_vars))

  codec <- new_state_codec(state_space, absorbing = absorbing)

  # Check that the matrix has a row for every combination of state values
  # recovered from its labels.
  if (!setequal(codec$transient, transient)) {
    micsim_abort(
      "micsim_compat_error",
      paste0(
        "The transient states of `transitionMatrix` are not the full\n",
        "  Cartesian product of their state variables. MicSim's state\n",
        "  space must be a complete grid; build it with expand.grid().\n",
        "  Missing: ", toString(utils_head(setdiff(codec$transient,
                                                  transient), 5L)), "\n",
        "  Unexpected: ", toString(utils_head(setdiff(transient,
                                                     codec$transient), 5L))
      )
    )
  }

  codec
}


#' Build the rate registry from a transition matrix
#'
#' @description
#' Records the origin state, destination state and rate function for each
#' cell whose value differs from `"0"`.
#'
#' @details
#' Within each row, transitions are registered in column order. The engine
#' uses that order to assign random draws to destinations. A function name
#' used in several cells is looked up once per registry.
#'
#' @param transitionMatrix As in [micSim()].
#' @param codec A `state_codec`.
#' @param caller Environment to resolve function names in.
#' @param probe Passed to `new_rate_registry()`.
#' @return A `rate_registry`.
#' @keywords internal
#' @noRd
registry_from_matrix <- function(transitionMatrix, codec, caller,
                                 probe = NULL) {

  transient <- rownames(transitionMatrix)
  destinations <- colnames(transitionMatrix)
  resolved <- new.env(parent = emptyenv())

  transitions <- list()
  for (i in seq_along(transient)) {
    cells <- transitionMatrix[i, ]
    for (j in which(cells != "0")) {
      transitions[[length(transitions) + 1L]] <- list(
        from = transient[i],
        to   = destinations[j],
        rate = resolve_rate(as.character(cells[j]), caller, resolved)
      )
    }
  }

  if (length(transitions) == 0L) {
    micsim_abort(
      "micsim_compat_error",
      paste0(
        "`transitionMatrix` allows no transitions at all: every cell is\n",
        "  \"0\". Nothing can happen to anybody."
      )
    )
  }

  new_rate_registry(transitions, codec, probe = probe)
}


#' Resolve a rate function name to the function itself
#'
#' @description
#' Looks the name up in the environment `micSim()` was called from, and
#' remembers the answer.
#'
#' @details
#' Naming a function in a matrix cell rather than passing the function
#' keeps `transitionMatrix` printable and readable. The cost is that the
#' name has to be resolved somewhere. Doing it once, at model construction,
#' in one named environment, is both fast and testable, and the answer is
#' cached so a name shared by many cells is looked up only once.
#'
#' @param name The function's name, as a string.
#' @param caller Environment to look in.
#' @param cache Environment used to remember names already resolved.
#' @return The function.
#' @keywords internal
#' @noRd
resolve_rate <- function(name, caller, cache) {

  found <- get0(name, envir = cache, inherits = FALSE)
  if (!is.null(found)) return(found)

  fn <- get0(name, envir = caller, mode = "function", inherits = TRUE)
  if (is.null(fn)) {
    micsim_abort(
      "micsim_compat_error",
      paste0(
        "`transitionMatrix` names a rate function `", name, "`, but no\n",
        "  function of that name is visible from where micSim() was\n",
        "  called. Define it before calling, in the same environment."
      )
    )
  }

  assign(name, fn, envir = cache)
  fn
}


#' Expand fertility patterns into explicit transition labels
#'
#' @description
#' Turns the `fertTr` argument, which may name only some of the
#' state variables, into the full list of transitions that create a
#' newborn.
#'
#' @details
#' A pair of states is a birth when all three of these hold:
#'
#' 1. both states are female, i.e. both contain the classic interface's
#'    label `"f"` (see `legacy_roles()`);
#' 2. every value named on the pattern's origin side appears is one of the 
#'    origin state's values (in an variable), and likewise for the destination;
#' 3. what the pattern does not mention is identical on both sides, so
#'    that a birth changes the fertility variable and nothing else.
#'
#' So `"0->1"` expands to a first birth for every nationality, but not to
#' a transition that also changes nationality, and not to one for a man.
#'
#' @section Assumptions:
#' Only women give birth, and in this interface `"f"` is the label that
#' says so.
#'
#' @param fertTr Character vector of `"origin->destination"` patterns.
#' @param codec A `state_codec`.
#' @return Character vector of transition labels in the registry's
#'   `"from -> to"` form.
#' @keywords internal
#' @noRd
expand_fertility_transitions <- function(fertTr, codec) {

  transient <- codec$transient
  parts <- strsplit(transient, "/", fixed = TRUE)
  female <- vapply(parts, function(p) legacy_roles()$female %in% p,
                   logical(1L))

  patterns <- strsplit(as.character(fertTr), "->", fixed = TRUE)
  labels <- character(0)

  for (pattern in patterns) {
    if (length(pattern) != 2L) {
      micsim_abort(
        "micsim_compat_error",
        paste0("`fertTr` entry `", paste(pattern, collapse = "->"),
               "` is not of the form \"origin->destination\".")
      )
    }
    origin_values      <- strsplit(trimws(pattern[1L]), "/", fixed = TRUE)[[1L]]
    destination_values <- strsplit(trimws(pattern[2L]), "/", fixed = TRUE)[[1L]]

    for (i in which(female)) {
      origin <- parts[[i]]
      if (!all(origin_values %in% origin)) next

      for (j in which(female)) {
        destination <- parts[[j]]
        if (!all(destination_values %in% destination)) next

        rest_origin      <- origin[!origin %in% origin_values]
        rest_destination <- destination[!destination %in% destination_values]
        if (!identical(rest_origin, rest_destination)) next

        labels <- c(labels, paste(transient[i], "->", transient[j]))
      }
    }
  }

  unique(labels)
}


#' head() without depending on utils
#'
#' Returns a short prefix of a vector for diagnostic messages.
#'
#' @param x A vector.
#' @param n How many elements to keep.
#' @return The first `n` elements, or all of them if there are fewer.
#' @keywords internal
#' @noRd
utils_head <- function(x, n) {
  x[seq_len(min(n, length(x)))]
}
