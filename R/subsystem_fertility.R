# =====================================================================
# subsystem_fertility.R
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   Adds a child whenever a transition selected as a birth occurs, and
#   records births for the simulation report.
#
# PLACE IN THE PIPELINE
#   fertility(births, newborn)  --add_subsystem()-->  compiled against the
#   model's declared transitions  --each run-->  engine.R fires on_event
#   --> subsystem_fertility.R adds and schedules a newborn --> engine.R
#
# WHAT YOU NEED TO KNOW
#   A subsystem adds behaviour at named points in the engine's event loop.
#   Births are ordinary state transitions; this subsystem supplies their
#   additional consequence of creating a person.
#
#   The newborn rule uses the mother's state after the birth transition.
#   No state label has a built-in meaning: the selected transitions define
#   who can give birth.
#
# GUARANTEES
#   * The birth callback creates one child per selected transition event.
#   * Every birth takes seven uniforms from the mother's birth block
#     (rng_streams.R): six seed the child's own stream, one draws the
#     child's state. A child's random numbers therefore depend on the
#     mother and the birth's position among her births, not on when or in
#     which order births happen in the population.
# =====================================================================


#' Add a child when a birth transition occurs
#'
#' @description
#' Creates a subsystem for [add_subsystem()]: whenever one of the `births`
#' transitions happens, a child is added to the population, with the
#' mother recorded and a state drawn from the `newborn` rule.
#'
#' @details
#' The births themselves are ordinary transitions with ordinary rates,
#' usually a parity variable going up by one:
#'
#' ```
#' transition(move("parity", "0", "1", where = sex == "female"), first_birth)
#' ```
#'
#' `fertility()` only says which of those transitions create a child.
#' Nothing about it depends on what the states are called: there is no
#' reserved label for women. Whoever has a birth transition is a mother.
#' If you want that checked, give `mothers`, and every birth transition
#' must start in a state satisfying it.
#'
#' At each birth the subsystem
#'
#' 1. draws the newborn's state from `newborn`, given the state the mother
#'    has just moved INTO;
#' 2. registers the child, born at the event's date, with `motherID` set
#'    and every duration clock starting at birth;
#' 3. schedules the child's first event.
#'
#' The result's `pop` then carries a `motherID` column.
#'
#' @section Report:
#' After a run, `result$extensions[[name]]` holds
#' \describe{
#'   \item{`n_births`}{the number of children born during the run;}
#'   \item{`births`}{one row per birth: `motherID`, `childID`, `date`
#'     (`yyyymmdd`), `motherAge` (years, two decimals), `transition` (the
#'     mother's birth transition, `"from -> to"`) and `childState` (the
#'     state the child was born into);}
#'   \item{`by_transition`}{births per birth transition, a named integer
#'     vector that lists transitions without births as zero;}
#'   \item{`by_child_state`}{births per state of the newborn.}
#' }
#' After a parallel run the identifiers are those of `pop`.
#'
#' @param births Declared transitions that create a child: a [move()]
#'   or [moves_of()] selector, a list of them, or a character vector of
#'   transition labels such as `"female/0->female/1"`.
#' @param newborn The state a child is born into: a [newborn_rule()], or a
#'   table-form rule from [new_newborn_rule()].
#' @param mothers Optional condition on the state variables, e.g.
#'   `sex == "female"`. When given, every birth transition must start in a
#'   state satisfying it, or the model is refused.
#' @param name The subsystem's name in the model and in
#'   `result$extensions`.
#'
#' @return A subsystem, for [add_subsystem()].
#'
#' @seealso [add_subsystem()], [newborn_rule()], [moves_of()]
#'
#' @examples
#' fertility(
#'   births  = moves_of("parity"),
#'   newborn = newborn_rule(sex = c(male = 0.515, female = 0.485),
#'                          parity = "0"),
#'   mothers = sex == "female"
#' )
#'
#' @export
fertility <- function(births, newborn, mothers = NULL, name = "fertility") {

  if (missing(births)) {
    abort_fertility(paste0(
      "fertility() needs `births`: the transitions that create a child,\n",
      "  e.g. births = moves_of(\"parity\")."
    ))
  }
  if (missing(newborn) ||
        !(inherits(newborn, "micsim_newborn_spec") ||
            inherits(newborn, "newborn_rule"))) {
    abort_fertility(paste0(
      "fertility() needs `newborn`: the state a child is born into, built\n",
      "  with newborn_rule() (or new_newborn_rule() for a full table)."
    ))
  }
  mothers <- new_role(substitute(mothers), parent.frame())

  new_subsystem(
    name   = name,
    type   = "fertility",
    params = list(births = births, newborn = newborn, mothers = mothers),
    spawns = TRUE,
    compile = function(ctx) {
      compile_fertility(births, newborn, mothers, ctx$codec, ctx$registry,
                        name)
    }
  )
}


#' Compile the fertility subsystem against a model
#'
#' @param births,newborn As in [fertility()].
#' @param mothers A captured `micsim_role`, or `NULL`.
#' @param codec A `state_codec`.
#' @param registry A `rate_registry`.
#' @param name The subsystem's name, for messages.
#' @return The compiled form; see `R/subsystem.R`.
#' @keywords internal
#' @noRd
compile_fertility <- function(births, newborn, mothers, codec, registry,
                              name) {

  what <- paste0("fertility(births = ) of subsystem `", name, "`")
  ids <- resolve_claim(births, codec, registry, what)
  edges <- registry$edges[ids]

  into_absorbing <- vapply(edges, function(e) codec$is_absorbing(e$to_idx),
                           logical(1L))
  if (any(into_absorbing)) {
    abort_fertility(paste0(
      what, " includes transition(s) into an absorbing state: ",
      toString(vapply(edges[into_absorbing], `[[`, character(1L), "label")),
      ".\n  A mother must stay in the population to give birth."
    ))
  }

  if (!is.null(mothers)) {
    may_mother <- states_where(codec, mothers,
                               paste0("`mothers` condition of `", name, "`"))
    outside <- !may_mother[vapply(edges, `[[`, integer(1L), "from_idx")]
    if (any(outside)) {
      abort_fertility(paste0(
        "Subsystem `", name, "`: birth transition(s) start in a state that\n",
        "  does not satisfy `mothers = ", deparse1(mothers$expr), "`: ",
        toString(utils_head(vapply(edges[outside], `[[`, character(1L),
                                   "label"), 5L)), "."
      ))
    }
  }

  rule <- resolve_newborn_rule(newborn, codec, name)

  is_birth <- logical(length(registry$edges))
  is_birth[ids] <- TRUE

  list(
    claims   = list(births = ids),
    describe = paste0("fertility (", length(ids), " birth transition(s))"),
    merge    = merge_reports,
    instantiate = function(run) {
      st <- new.env(parent = emptyenv())
      st$n_births <- 0L
      st$birth_ids <- ids
      st$log <- new_event_log(list(motherID = integer(), childID = integer(),
                                   time_days = numeric(),
                                   mother_birth_days = numeric(),
                                   edge_id = integer(), child_idx = integer()))
      list(
        hooks  = fertility_hooks(is_birth, rule, st),
        finish = function() fertility_report(st, registry, codec)
      )
    }
  )
}


#' The fertility subsystem's hooks, for one run
#'
#' @description
#' Returns the `on_event` callback that turns a birth transition into a
#' newborn.
#'
#' @section Assumptions:
#' **Random numbers.** Each birth takes seven uniforms from the mother's
#' birth block: six seed the child's own stream, the seventh draws the
#' child's state. A child's random numbers therefore depend on the mother
#' and on the birth's position among her births, not on the order
#' of events in the population.
#'
#' **The mother's state is the destination, not the origin.** A woman who
#' has just moved from `"f/AT/0"` to `"f/AT/1"` is in `"f/AT/1"` when her
#' child's attributes are inherited. This matters only if an inherited
#' variable is the one the birth changed.
#'
#' @param is_birth Logical vector over the registry's edges: `TRUE` for a
#'   birth.
#' @param rule A `newborn_rule`.
#' @param st The run's state environment; `n_births` is counted in it.
#'
#' @return An `engine_hooks` object with `on_event` set.
#'
#' @keywords internal
#' @noRd
fertility_hooks <- function(is_birth, rule, st) {

  on_event <- function(ctx) {

    if (!is_birth[ctx$edge$id]) {
      return(invisible(NULL))
    }

    mother_slot  <- ctx$slot
    mother_state <- ctx$codec$decode(ctx$to_idx)

    # (1) the birth's random numbers, from the mother's birth block: six
    #     seed the child's stream, the seventh draws the child's state.
    u <- birth_uniforms(ctx$persons, mother_slot, 7L)
    newborn_state <- assign_newborn_state(rule, mother_state, u = u[7L])

    # (2) register the child. Their duration clocks all start at the birth
    #     date, which for someone born during the run is the literal truth
    #     rather than an assumption.
    newborn_idx  <- ctx$codec$encode(newborn_state)
    newborn_id   <- ctx$persons$next_id()
    newborn_slot <- ctx$persons$add(
      id             = newborn_id,
      birth_days     = ctx$now_days,
      state_idx      = newborn_idx,
      init_state_idx = newborn_idx,
      mother_id      = ctx$persons$id_of(mother_slot),
      rng_state      = state_from_uniforms(u[1:6])
    )
    st$n_births <- st$n_births + 1L
    st$log$add(motherID = ctx$persons$id_of(mother_slot),
               childID = as.integer(newborn_id),
               time_days = ctx$now_days,
               mother_birth_days = ctx$persons$birth_days_of(mother_slot),
               edge_id = ctx$edge$id, child_idx = newborn_idx)

    # (3) tell anyone else who is interested.
    fire_hook(ctx$hooks, "on_birth",
              new_hook_context(persons = ctx$persons, queue = ctx$queue,
                               store = ctx$store, codec = ctx$codec,
                               registry = ctx$registry, hooks = ctx$hooks,
                               schedule = ctx$schedule,
                               now_days = ctx$now_days,
                               slot = newborn_slot,
                               mother_slot = mother_slot))

    # (4) bring the child into the running simulation.
    ctx$schedule(newborn_slot, ctx$now_days)

    invisible(NULL)
  }

  new_engine_hooks(on_event = on_event)
}


#' Build the fertility report from a run's birth log
#'
#' @description
#' Builds `result$extensions$fertility` (or whatever the subsystem is
#' named) from the run's event log.
#'
#' @param st Per-run environment containing the birth count, selected edge
#'   IDs and event log.
#' @param registry,codec The model's rate registry and state codec.
#' @return A list with `n_births`, `births`, `by_transition` and
#'   `by_child_state`, as described in [fertility()].
#' @keywords internal
#' @noRd
fertility_report <- function(st, registry, codec) {
  cols <- st$log$columns()
  labels <- vapply(registry$edges, `[[`, character(1L), "label")
  births <- data.frame(
    motherID   = as.integer(cols$motherID),
    childID    = as.integer(cols$childID),
    date       = as.character(getInDateFormat(cols$time_days)),
    motherAge  = round((cols$time_days - cols$mother_birth_days) / DAYS_PER_YEAR, 2),
    transition = labels[cols$edge_id],
    childState = codec$decode(cols$child_idx),
    stringsAsFactors = FALSE
  )
  claimed <- labels[st$birth_ids]
  list(
    n_births       = st$n_births,
    births         = births,
    by_transition  = count_by(births$transition, claimed),
    by_child_state = count_by(births$childState)
  )
}


#' Count occurrences, optionally over a fixed set of categories
#' @param x Character vector.
#' @param levels Categories to list even when absent, or `NULL`.
#' @return Named integer vector, in the order of `levels` (or sorted).
#' @keywords internal
#' @noRd
count_by <- function(x, levels = NULL) {
  levels <- if (is.null(levels)) sort(unique(x)) else unique(c(levels, sort(unique(x))))
  stats::setNames(tabulate(match(x, levels), nbins = length(levels)), levels)
}


#' Raise a fertility error
#' @param problem The message.
#' @return Never returns.
#' @keywords internal
#' @noRd
abort_fertility <- function(problem) {
  micsim_abort("micsim_fertility_error", problem)
}
