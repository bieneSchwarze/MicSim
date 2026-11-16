# =====================================================================
# engine.R
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   Runs the simulation. This is the discrete-event loop: take the soonest
#   scheduled event, make it happen, record it, schedule what comes next,
#   repeat until the horizon closes or nothing is left to do.
#
# PLACE IN THE PIPELINE
#   a model spec
#        |
#        v
#   engine.R  ---- asks ---->  registry.R  (what can happen from here?)
#        |       ---- asks ---->  solver.R   (when would it happen?)
#        |       ---- reads --->  engine_person.R (how old, how long?)
#        |       ---- uses ---->  engine_queue.R  (what happens next?)
#        |       ---- writes -->  engine_store.R  (transition records)
#        |       ---- calls --->  engine_hooks.R  (extension behaviour)
#        v
#   the `pop` data frame
#
# WHAT YOU NEED TO KNOW  (concepts in ?MicSim-concepts)
#   discrete-event simulation, event queue, global clock, competing risks,
#   the three time scales, absorbing state.
#
#   The shape of the algorithm:
#
#     for everyone alive at the start:  schedule their next event
#     loop over the following:
#       take the soonest event off the queue; it sets the global clock
#       stop if that takes us past the end of the simulation
#       record the transition and move the individual into the new state
#       if the new state is not absorbing, schedule their next event
#
#   Everything else in this file is bookkeeping around this.
#
# GUARANTEES
#   * The global clock never moves backwards. No subsystem hook may move
#     an event before the present (on_candidates and on_schedule are
#     both checked), and the event loop stops with an internal error if it
#     ever happens.
#   * Every individual draws from their own random-number stream
#     (rng_streams.R), so a given seed always produces the same run, on any
#     number of cores, and an individual's random numbers do not depend on
#     anybody else or on the order in which events are processed.
#   * The caller's random-number state is restored on exit. Without a
#     seed, the run takes exactly one draw from it, the seed it then uses.
#
# UNITS
#   Event dates, ages and durations are stored in days as produced by
#   getInDays(). Waiting-time calculations and rate functions use years,
#   with DAYS_PER_YEAR for conversion. This conversion only occurs when
#   building a hazard. The suffixes `_days` and `_years` identify these units.
#   see ?MicSim-concepts, "Units in variable names".
#
# EXTENSIONS AND IMMIGRATION
#   Subsystems, such as fertility, add behaviour through hooks and are
#   instantiated at the start of each run; see subsystem.R. Immigration
#   is handled directly here: individuals are registered from an arrival
#   table, and their first transitions are scheduled from arrival.
# =====================================================================


#' Create a simulation engine
#'
#' @description
#' Combines the state mapping, transition rates, population data and
#' simulation settings into an engine whose `run()` method executes a
#' simulation.
#'
#' @details
#' The engine can be run repeatedly with the same inputs. Each call to
#' `run()` creates new population, event and transition stores and new
#' subsystem instances. The supplied codec, registry and additional hooks
#' are retained by the engine; additional hooks may hold their own state.
#'
#' **The initial population** is the `initPop` data frame: columns `ID`,
#' `birthDate` and `initState`, in that order. Birth dates are `yyyymmdd`.
#' Every individual must be born on or before the simulation starts and be
#' younger than `max_age` at that moment; both are checked here, with the
#' offending identifiers named.
#'
#' **Random numbers.** Every individual draws from their own L'Ecuyer-CMRG
#' stream of random numbers, derived from the seed and their ID, or for
#' a newborn from the mother's stream (see `R/rng_streams.R`). Within one
#' scheduling, destinations receive their uniforms in the order of the
#' state space. The order in which individuals are registered or scheduled,
#' and the order in which transitions were declared, do not affect
#' the random numbers.
#'
#' @section Assumptions:
#' Ordinary waiting-time calculations use the shorter of the time to the
#' simulation end and the time to age `max_age - MAX_AGE_MARGIN`.
#' The margin keeps this calculation slightly below `max_age`, which is
#' useful for rate tables indexed by whole age. Subsystem hooks can adjust
#' event times after waiting times have been drawn.
#'
#' @param codec A `state_codec`. Defines the state space.
#' @param registry A `rate_registry` object containing the allowed
#'   transitions and their rate functions.
#' @param init_pop Non-empty starting population data frame with columns
#'   `ID`, `birthDate` and `initState`, in that order. Dates use `yyyymmdd`.
#' @param horizon Length-2 vector `c(start, end)` of `yyyymmdd` dates. The
#'   start must precede the end.
#' @param max_age Maximum age for scheduling transitions, in years.
#' @param subsystems List of compiled subsystems (from
#'   `compile_subsystem()`, as [add_subsystem()] stores them), in the
#'   order they take effect. Each run instantiates them afresh.
#' @param immigrants Data frame of individuals entering the population
#'   after the simulation starts, with columns `ID`, `immigrDate`,
#'   `birthDate` and `immigrInitState`, in that order. Dates are
#'   `yyyymmdd`.
#' @param hooks Additional `engine_hooks` object, combined after subsystem
#'   hooks. Intended for testing and engine development; model extensions
#'   are normally supplied through `subsystems`.
#' @param verbose Logical; whether to report progress through [message()].
#'   The compatibility shim passes `TRUE` and the object API defaults to
#'   `FALSE`.
#'
#' @return An object of class `micsim_engine`, a list with
#'   \describe{
#'     \item{`run(seed)`}{execute the simulation. Returns a list with `pop`,
#'       `newborns` (identifiers and exact birth days of the individuals
#'       born during the run, used to number newborns in multi-core runs),
#'       `diagnostics` (including `seed`, the seed used, and `seed_drawn`)
#'       and `extensions`, the last holding what each subsystem reported,
#'       by name.}
#'     \item{`codec`, `registry`}{the collaborators, for inspection}
#'   }
#'
#' @seealso [MicSim-concepts]; `new_rate_registry()`, `new_person_store()`,
#'   `new_event_queue()`, `new_transition_store()`, `draw_waiting_time()`.
#'
#' @examples
#' codec <- new_state_codec(list(sex = c("m", "f")), absorbing = "dead")
#' mortality <- function(age, calTime) 0.00003 * exp(0.1 * age)
#' registry <- new_rate_registry(
#'   list(list(from = NA, to = "dead", rate = mortality)), codec
#' )
#' init_pop <- data.frame(
#'   ID = 1:3,
#'   birthDate = c("19301231", "19560101", "19910704"),
#'   initState = c("f", "f", "m"),
#'   stringsAsFactors = FALSE
#' )
#' engine <- new_engine(codec, registry, init_pop,
#'                      horizon = c(20000101, 20501231), max_age = 120)
#' result <- engine$run(seed = 234)
#' result$pop
#'
#' @keywords internal
#' @noRd
new_engine <- function(codec,
                       registry,
                       init_pop,
                       horizon,
                       max_age,
                       subsystems = list(),
                       immigrants = NULL,
                       hooks      = new_engine_hooks(),
                       verbose    = FALSE) {

  init_pop   <- standardise_people(init_pop)
  immigrants <- standardise_people(immigrants)
  validate_engine_inputs(codec, registry, init_pop, horizon, max_age, hooks)
  validate_engine_subsystems(subsystems, registry)
  check_initial_population(init_pop, codec, getInDays(horizon[1]), max_age)
  if (!is.null(immigrants)) {
    validate_immigrants(immigrants, horizon, max_age, codec, init_pop$ID)
  }

  # Subsystems extend the engine through hooks that can adjust hazards,
  # candidate times and scheduled events, or respond to transitions.
  # Their hooks are created for each run and combined in subsystem order,
  # followed by any additional hooks supplied directly to the engine.
  
  # Include motherID in the output whenever a subsystem declares that
  # it can add individuals, as the fertility subsystem does.
  include_mother <- any(vapply(subsystems, function(s)
    isTRUE(s$spawns), logical(1L)))

  start_days <- getInDays(horizon[1])
  end_days   <- getInDays(horizon[2])

  # The random numbers of one scheduling are drawn in the order of the
  # destinations in the state space, so that a destination receives the
  # same uniform whatever order the transitions were declared in, and
  # whatever the states are called. budget_position[[i]][j] is the rank,
  # in that order, of the destination of the j-th edge leaving state i
  # (see rng_streams.R).
  budget_position <- lapply(seq_len(codec$n), function(i) {
    edges <- registry$feasible(i)
    if (is.null(edges)) return(NULL)
    to_idx <- vapply(edges, `[[`, integer(1L), "to_idx")
    rank <- integer(length(edges))
    rank[order(to_idx)] <- seq_along(edges)
    rank
  })

  # ------------------------------------------------------------------
  # Schedule the next event for one individual.
  #
  # Draw a waiting time for each allowed transition, apply scheduling
  # hooks and queue the earliest remaining event. See ?MicSim-concepts,
  # "Competing risks".
  #
  # `slot` is the individual's position in the person store. `now_days`
  # is the date from which scheduling begins: the simulation start,
  # arrival, a transition or a date supplied by a subsystem. The stores,
  # counters and run hooks are passed explicitly; model settings are
  # taken from the enclosing new_engine() call.
  # ------------------------------------------------------------------
  schedule_next <- function(persons, queue, counters, hooks, slot,
                            now_days) {

    state_idx <- persons$state_of(slot)
    if (codec$is_absorbing(state_idx)) return(invisible(NULL))

    edges <- registry$feasible(state_idx)
    if (is.null(edges)) return(invisible(NULL))

    birth_days <- persons$birth_days_of(slot)
    age_years  <- (now_days - birth_days) / DAYS_PER_YEAR

    # Limit the ordinary waiting-time calculation to the earlier of the
    # simulation end and the age boundary, including its small margin.
    # Scheduling hooks may subsequently adjust the chosen event time.
    years_to_horizon <- (end_days - now_days) / DAYS_PER_YEAR
    years_to_max_age <- (max_age - MAX_AGE_MARGIN) - age_years
    horizon_years <- min(years_to_horizon, years_to_max_age)

    # One random budget per destination, all drawn now from this person's
    # next substream (rng_streams.R). Every destination gets its budget
    # whether or not it is used (a forced transition, a zero rate or no
    # window left), so nothing else ever shifts.
    budgets <- budgets_from_uniforms(
      next_event_uniforms(persons, slot, length(edges)),
      budget_position[[state_idx]]
    )

    caltime_years <- caltime_in_years(now_days)

    best_years <- Inf
    best_to    <- NA_integer_
    best_id    <- NA_integer_

    # Retain all candidates when an on_candidates hook is active, so it
    # can adjust their times before the winner is selected. Candidate
    # collection itself does not change the waiting-time draws.
    want_candidates <- hook_active(hooks, "on_candidates")
    if (want_candidates) {
      cand_years <- numeric(length(edges))
      cand_to    <- integer(length(edges))
      cand_id    <- integer(length(edges))
      k <- 0L
    }

    for (j in seq_along(edges)) {
      edge <- edges[[j]]
      # Look up duration separately for each transition because its
      # affected state variable determines which duration clock is used. 
      # Inf when the transition has no single affected variable, 
      # which covers every absorbing destination.
      duration_days <- persons$duration_days(slot, edge$affected_var,
                                             now_days)

      hazard <- bind_hazard(edge,
                            age0_years      = age_years,
                            caltime0_years  = caltime_years,
                            duration0_years = duration_days / DAYS_PER_YEAR)
      hazard <- adjust_by_hook(
        hooks, "modify_hazard", hazard,
        new_hook_context(slot = slot, edge = edge, now_days = now_days,
                         persons = persons, codec = codec)
      )

      drawn <- draw_waiting_time(
        lambda        = hazard$lambda,
        horizon_years = horizon_years,
        forced_time   = hazard$forced_time,
        crossings     = hazard$crossings,
        label         = edge$label,
        budget        = budgets[j]
      )

      if (!drawn$reliable) counters$unreliable <- counters$unreliable + 1L
      if (drawn$forced)    counters$forced     <- counters$forced + 1L

      if (want_candidates) {
        k <- k + 1L
        cand_years[k] <- drawn$time
        cand_to[k]    <- edge$to_idx
        cand_id[k]    <- edge$id
      }

      # Replace the current choice only with a strictly earlier draw.
      # Ties therefore retain the first transition in registry order.
      if (drawn$time < best_years) {
        best_years <- drawn$time
        best_to    <- edge$to_idx
        best_id    <- edge$id
      }
    }

    if (want_candidates) {
      # Competing risks are resolved AFTER a subsystem has moved
      # candidates, so a calendar rule that postpones one transition lets
      # every other transition compete over the postponement. Moving a
      # candidate before now is refused: the global clock never goes back.
      cand <- adjust_by_hook(
        hooks, "on_candidates",
        list(time_days = now_days + cand_years * DAYS_PER_YEAR,
             to_idx = cand_to, edge_id = cand_id),
        new_hook_context(slot = slot, now_days = now_days,
                         persons = persons, codec = codec)
      )
      if (any(cand$time_days < now_days)) {
        abort_engine(paste0(
          "A subsystem moved a scheduled event to before the present\n",
          "  (day ", min(cand$time_days), " < day ", now_days, "). Events may\n",
          "  be postponed or brought forward, but never into the past."
        ))
      }
      if (!any(is.finite(cand$time_days))) return(invisible(NULL))
      # Strictly earlier wins, as below: which.min() takes the first.
      win <- which.min(cand$time_days)
      event <- list(time_days = cand$time_days[win],
                    to_idx    = cand$to_idx[win],
                    edge_id   = cand$edge_id[win])
    } else {
      if (!is.finite(best_years)) return(invisible(NULL))

      # `edge_id` tells an adjusting hook WHICH transition was drawn, as an
      # integer it can look up without comparing labels. The queue keeps
      # only the date and the destination.
      event <- list(time_days = now_days + best_years * DAYS_PER_YEAR,
                    to_idx    = best_to,
                    edge_id   = best_id)
    }
    event <- adjust_by_hook(
      hooks, "on_schedule", event,
      new_hook_context(slot = slot, now_days = now_days, persons = persons,
                       queue = queue, codec = codec)
    )
    if (is.null(event)) return(invisible(NULL))
    # The same rule as for on_candidates: an adjusted event may be
    # postponed or brought forward, but never moved into the past.
    if (event$time_days < now_days) {
      abort_engine(paste0(
        "A subsystem moved a scheduled event to before the present\n",
        "  (day ", event$time_days, " < day ", now_days, "). Events may\n",
        "  be postponed or brought forward, but never into the past."
      ))
    }

    queue$push(time_days = event$time_days, slot = slot,
               to_idx = event$to_idx)
    if (queue$size() > counters$peak_queue) {
      counters$peak_queue <- queue$size()
    }

    invisible(NULL)
  }

  # ------------------------------------------------------------------
  # Execute one run.
  # ------------------------------------------------------------------
  run <- function(seed = NULL) {

    started_at <- Sys.time()

    # Every individual draws from their own stream, derived from the seed
    # (rng_streams.R). Without a seed, one is drawn from the caller's
    # stream. Either way the caller's random-number state is put back on
    # exit, so the run itself leaves no trace on it.
    rng <- open_run_rng(seed)
    on.exit(rng$restore(), add = TRUE)

    persons <- new_person_store(n_vars = length(codec$vars))
    queue   <- new_event_queue()
    store   <- new_transition_store()

    # Every subsystem gets fresh hooks and fresh state for this run, so
    # nothing carries over from a previous run of the same model. Core
    # subsystems first, in the order they were added; raw `hooks` last.
    instances <- instantiate_subsystems(
      subsystems, list(codec = codec, registry = registry)
    )
    run_hooks <- do.call(combine_hooks,
                         c(unname(lapply(instances, `[[`, "hooks")),
                           list(hooks)))

    # Share diagnostic counters with scheduling helpers through an
    # environment, so updates affect this run's counters directly.
    counters <- new.env(parent = emptyenv())
    counters$unreliable      <- 0L
    counters$forced          <- 0L
    counters$peak_queue      <- 0L

    # The scheduling primitive, handed to subsystems through `ctx`. It is
    # the engine's own schedule_next() with the run's stores already
    # bound, so a subsystem can bring a new individual into the running
    # simulation without knowing anything about how scheduling works.
    schedule <- function(slot, now_days) {
      schedule_next(persons, queue, counters, run_hooks, slot, now_days)
    }

    if (verbose) message("Registering the starting population ...")
    seed_initial_population(persons, init_pop, codec, start_days, max_age,
                            rng_states = person_stream_states(rng$master,
                                                              init_pop$ID))

    # Each individual draws from their own stream, so the order in which
    # they are scheduled does not affect the random numbers they receive.
    if (verbose) message("Scheduling first events ...")
    for (slot in seq_len(persons$size())) {
      schedule_next(persons, queue, counters, run_hooks, slot, start_days)
    }

    # Register and schedule immigrants after the starting population,
    # preserving immigrant-table row order. Registration happens here
    # even for later arrivals; each first event is drawn from arrival.
    n_immigrants <- 0L
    if (!is.null(immigrants)) {
      if (verbose) message("Registering immigrants ...")
      slots <- seed_immigrant_population(
        persons, immigrants, codec, max_age,
        rng_states = person_stream_states(rng$master, immigrants$ID)
      )
      n_immigrants <- length(slots)
      arrival_days <- getInDays(immigrants$immigrDate)
      for (k in seq_along(slots)) {
        # An immigrant's clock starts at their ARRIVAL, not at the start of
        # the simulation: their age, the calendar time their rates see and
        # the window they have left are all measured from the day they
        # enter the population.
        schedule_next(persons, queue, counters, run_hooks, slots[k],
                      arrival_days[k])
      }
    }

    # ---- the discrete-event loop ----------------------------------
    if (verbose) message("Simulation is running ...")
    clock_days   <- start_days
    reported_year <- getYear(start_days)

    repeat {
      event <- queue$pop_min()
      if (is.null(event)) break                    # no scheduled events remain

      now_days <- event$time_days
      if (now_days > end_days) break               # beyond the simulation end

      # The queue is a min-heap and every drawn event lies in the future;
      # on_candidates may not move an event before now, and on_schedule is
      # checked the same way when the event is scheduled. The clock can
      # therefore only move backwards through a bug in this file, and the
      # loop says so rather than carrying on.
      if (now_days < clock_days) {
        abort_engine(paste0(
          "The global clock moved backwards, from day ", clock_days,
          " to day ", now_days, ".\n",
          "  This is an internal error, not a problem with the model."
        ))
      }
      clock_days <- now_days

      if (verbose) {
        current_year <- getYear(now_days)
        if (current_year > reported_year) {
          message("  year ", current_year)
          reported_year <- current_year
        }
      }

      slot     <- event$slot
      to_idx   <- event$to_idx
      from_idx <- persons$state_of(slot)
      edge     <- find_edge(registry, from_idx, to_idx)

      store$append(slot      = slot,
                   from_idx  = from_idx,
                   to_idx    = to_idx,
                   time_days = now_days,
                   age_days  = persons$age_days(slot, now_days))

      persons$move_to(slot, to_idx, now_days, edge$changed_vars)

      fire_hook(run_hooks, "on_event",
                new_hook_context(persons = persons, queue = queue,
                                 store = store, codec = codec,
                                 registry = registry, hooks = run_hooks,
                                 schedule = schedule, now_days = now_days,
                                 slot = slot, from_idx = from_idx,
                                 to_idx = to_idx, edge = edge))

      if (codec$is_absorbing(to_idx)) {
        fire_hook(run_hooks, "on_absorb",
                  new_hook_context(persons = persons, queue = queue,
                                   store = store, codec = codec,
                                   registry = registry, hooks = run_hooks,
                                   schedule = schedule, now_days = now_days,
                                   slot = slot, from_idx = from_idx,
                                   to_idx = to_idx, edge = edge))
      } else {
        schedule_next(persons, queue, counters, run_hooks, slot, now_days)
      }
    }

    if (verbose) message("Simulation has finished.")

    # Individuals born during the run, with their exact birth times: a
    # multi-core run uses them to number newborns in the order a
    # single-core run would have.
    people <- persons$as_columns()
    born <- !is.na(people$motherID)

    list(
      pop = assemble_pop(store, persons, codec,
                         include_mother = include_mother),
      newborns = list(ID = people$ID[born], birth_days = people$birthDays[born]),
      diagnostics = list(
        n_events          = store$size(),
        n_individuals     = persons$size(),
        n_immigrants      = n_immigrants,
        n_unreliable      = counters$unreliable,
        n_forced          = counters$forced,
        peak_queue        = counters$peak_queue,
        seed              = rng$seed,
        seed_drawn        = rng$drawn,
        rng_kind          = "L'Ecuyer-CMRG",
        seconds           = as.numeric(difftime(Sys.time(), started_at,
                                                units = "secs"))
      ),
      extensions = finish_subsystems(instances)
    )
  }

  structure(
    list(run = run, codec = codec, registry = registry),
    class = "micsim_engine"
  )
}


#' Check the basic inputs required by the engine
#'
#' @description
#' Checks the input classes, starting-population structure, simulation
#' dates and maximum age before the run's data stores are created.
#'
#' @details
#' The codec must contain at least one state variable, and the starting
#' population must be non-empty with the expected columns in the expected
#' order. The simulation dates must be valid and increasing, and `max_age`
#' must be a single positive number.
#'
#' Population values, immigrants and subsystems are checked separately by
#' the other validation helpers called in `new_engine()`. This function
#' does not evaluate transition rates.
#'
#' @param codec,registry,init_pop,horizon,max_age,hooks As in `new_engine()`.
#' @return Invisibly returns `TRUE` if these checks pass. Otherwise,
#'   raises a `micsim_engine_error`.
#' @keywords internal
#' @noRd
validate_engine_inputs <- function(codec, registry, init_pop, horizon,
                                   max_age, hooks) {

  if (!inherits(codec, "state_codec")) {
    abort_engine("`codec` must be a state_codec.")
  }
  if (length(codec$vars) == 0L) {
    abort_engine(paste0(
      "The state space has no named variables.\n",
      "  Build the codec from a NAMED list, e.g.\n",
      "  new_state_codec(list(sex = c('m', 'f')), absorbing = 'dead').\n",
      "  The names are what the duration clocks are indexed by."
    ))
  }
  if (!inherits(registry, "rate_registry")) {
    abort_engine("`registry` must be a rate_registry.")
  }
  if (!inherits(hooks, "engine_hooks")) {
    abort_engine("`hooks` must be an engine_hooks object.")
  }

  if (!is.data.frame(init_pop) ||
        !identical(names(init_pop), c("ID", "birthDate", "initState"))) {
    abort_engine(paste0(
      "`init_pop` must be a data frame with columns ID, birthDate and\n",
      "  initState, in that order."
    ))
  }
  if (nrow(init_pop) == 0L) {
    abort_engine("`init_pop` is empty: there is nobody to simulate.")
  }

  if (length(horizon) != 2L) {
    abort_engine("`horizon` must be c(start, end), both as yyyymmdd.")
  }
  horizon_days <- getInDays(horizon)
  if (anyNA(horizon_days)) {
    abort_engine("`horizon` contains a date that is not a valid yyyymmdd.")
  }
  if (horizon_days[1L] >= horizon_days[2L]) {
    abort_engine("`horizon` must start strictly before it ends.")
  }

  if (!is.numeric(max_age) || length(max_age) != 1L || is.na(max_age) ||
        max_age <= 0) {
    abort_engine("`max_age` must be a single positive number of years.")
  }

  invisible(TRUE)
}


#' Check compiled subsystems before a run
#'
#' @description
#' Checks that `subsystems` is a list of compiled subsystem objects and
#' that every transition identifier they claim is within the registry's
#' range.
#' 
#' @param subsystems List of compiled subsystems.
#' @param registry The `rate_registry` they must have been compiled
#'   against.
#' @return `invisible(TRUE)`, or raises `micsim_engine_error`.
#' @keywords internal
#' @noRd
validate_engine_subsystems <- function(subsystems, registry) {
  if (!is.list(subsystems)) {
    abort_engine("`subsystems` must be a list of compiled subsystems.")
  }
  for (s in subsystems) {
    if (!inherits(s, "micsim_compiled_subsystem")) {
      abort_engine(paste0(
        "Every element of `subsystems` must be a compiled subsystem, as\n",
        "  add_subsystem() stores it. Add subsystems to a model with\n",
        "  add_subsystem(model, fertility(...))."
      ))
    }
    n_edges <- length(registry$edges)
    for (ids in s$claims) {
      if (length(ids) > 0L && (max(ids) > n_edges || min(ids) < 1L)) {
        abort_engine(paste0(
          "Subsystem `", s$name, "` was compiled against a different set of\n",
          "  transitions than the model now has."
        ))
      }
    }
  }
  invisible(TRUE)
}


#' Register the starting population
#'
#' @description
#' Checks the starting population and adds each individual to the person
#' store, preserving input row order.
#'
#' @details
#' Two conditions make an individual impossible, and the engine stops on
#' both, listing the offending identifiers:
#'
#' * born after the simulation starts, which would give them a negative
#'   age;
#' * already older than `max_age` at the start, which means they should
#'   never have been at risk in the first place.
#'
#' The `birthDate` string is kept exactly as supplied. The output contract
#' returns it unchanged, so re-deriving it from the day count would alter
#' the user's own data.
#'
#' @param persons A `person_store` object to which individuals are added.
#' @param init_pop The starting population data frame.
#' @param codec A `state_codec` object used to encode initial states.
#' @param start_days Simulation start date on the `getInDays()` day scale.
#' @param max_age Maximum permitted age at entry, in years.
#' @param rng_states Integer matrix with one row of six per individual,
#'   from `person_stream_states()`; `NULL` leaves the streams unset
#'   (unit tests).
#'
#' @return The updated `persons` object, invisibly.
#'
#' @keywords internal
#' @noRd
seed_initial_population <- function(persons, init_pop, codec, start_days,
                                    max_age, rng_states = NULL) {

  state_idx  <- check_initial_population(init_pop, codec, start_days, max_age)
  birth_days <- getInDays(init_pop$birthDate)

  for (i in seq_len(nrow(init_pop))) {
    persons$add(
      id             = init_pop$ID[i],
      birth_days     = birth_days[i],
      state_idx      = state_idx[i],
      init_state_idx = state_idx[i],
      birth_ymd      = as.character(init_pop$birthDate[i]),
      rng_state      = if (is.null(rng_states)) NULL else rng_states[i, ]
    )
  }

  invisible(persons)
}



#' Read a population or immigrant table by its values.
#'
#' Turns factor columns into character and Date columns into yyyymmdd
#' integers, so that every later step reads the values the user sees. A
#' factor is stored as integer codes, and as.numeric() returns the codes,
#' not the labels: identifiers factor(c("10", "20", "30")) would become 1, 2
#' and 3. A Date is read correctly by getInDays(), but would be copied into
#' the output in its own format ("1980-01-01") instead of yyyymmdd.
#' Applied wherever a table enters: micsim_model(), update(), new_engine(),
#' micSim() and micSimParallel(). Anything that is not a data frame is
#' returned unchanged, for the input checks to report.
#'
#' @param people The data frame to standardise.
#'
#' @return The standardised data frame.
#'
#' @keywords internal
#' @noRd
standardise_people <- function(people) {
  if (!is.data.frame(people)) return(people)
  for (col in names(people)) {
    v <- people[[col]]
    if (is.factor(v)) {
      people[[col]] <- as.character(v)
    } else if (inherits(v, "Date")) {
      people[[col]] <- getInDateFormat(getInDays(v))
    }
  }
  people
}


#' Check the starting population against the model
#'
#' @description
#' Validates the starting population's dates, entry ages, identifiers and
#' initial states, and returns the encoded state indices.
#'
#' @details
#' The conditions, each reported with every offending identifier:
#' an invalid birth date; born after the simulation starts; older than
#' `max_age` at the start; an identifier that is missing, not a whole
#' number, or repeated; an unknown state; an absorbing state. All conditions
#' are checked and reported together. The absorbing-state check needs known
#' labels and is skipped when a label is unknown.
#'
#' @param init_pop The starting population data frame.
#' @param codec A `state_codec`.
#' @param start_days Simulation start date on the `getInDays()` day scale.
#' @param max_age Maximum permitted age at entry, in years.
#' @return The individuals' state indices, invisibly.
#' @keywords internal
#' @noRd
check_initial_population <- function(init_pop, codec, start_days, max_age) {

  # Every check runs and all faults are reported together (see
  # new_problem_collector()).
  found <- new_problem_collector()

  birth_days <- getInDays(init_pop$birthDate)
  if (anyNA(birth_days)) {
    bad <- init_pop$ID[is.na(birth_days)]
    found$add(paste0(
      "The starting population contains birth dates that are not valid\n",
      "  yyyymmdd dates. Affected IDs: ", show_ids(bad), "."
    ))
  }

  # Ages only for the valid birth dates.
  age_days <- start_days - birth_days
  born_later <- !is.na(age_days) & age_days < 0
  if (any(born_later)) {
    bad <- init_pop$ID[born_later]
    found$add(paste0(
      "The starting population contains individuals born after the\n",
      "  simulation begins. Affected IDs: ", show_ids(bad), "."
    ))
  }
  too_old <- !is.na(age_days) & age_days / DAYS_PER_YEAR > max_age
  if (any(too_old)) {
    bad <- init_pop$ID[too_old]
    found$add(paste0(
      "The starting population contains individuals already older than\n",
      "  max_age (", max_age, " years) when the simulation begins.\n",
      "  Affected IDs: ", show_ids(bad), "."
    ))
  }

  found$try(check_identifiers(init_pop$ID, "starting population"))

  # encode() raises on an unknown label, naming it, which is the right
  # error for a typo in a state name. The absorbing check needs the
  # indices, so it runs only when every label is known.
  state_idx <- found$try(codec$encode(as.character(init_pop$initState)))
  if (!is.null(state_idx)) {
    found$try(refuse_absorbing_entry(codec, state_idx, init_pop$ID,
                                     "starting population"))
  }
  found$raise("The starting population")
  invisible(state_idx)
}


#' Collect several faults of one table and report them as one error.
#'
#' @description
#' add(message) records a fault of class `micsim_engine_error`; try(expr)
#' evaluates expr and records the MicSim error it raises, if any, with its
#' class, returning expr's value or NULL. raise(what) does nothing without
#' faults; with one fault it raises it unchanged, with several it raises
#' one error listing all of them, carrying the classes of every fault and
#' the messages as `$problems`.
#'
#' @return A list with three functions: add(message), try(expr) and raise(what).
#'
#' @keywords internal
#' @noRd
new_problem_collector <- function() {
  e <- new.env(parent = emptyenv())
  e$messages <- character(0)
  e$classes <- character(0)
  add <- function(message, class = "micsim_engine_error") {
    e$messages <- c(e$messages, message)
    e$classes <- c(e$classes, class)
    invisible(NULL)
  }
  try_check <- function(expr) {
    tryCatch(expr, error = function(err) {
      cls <- setdiff(class(err), c("error", "condition"))
      if (!any(startsWith(cls, "micsim_"))) stop(err)
      add(conditionMessage(err), cls[1L])
      NULL
    })
  }
  raise <- function(what) {
    n <- length(e$messages)
    if (n == 0L) return(invisible(NULL))
    if (n == 1L) micsim_abort(e$classes, e$messages)
    micsim_abort(
      unique(e$classes),
      paste0(what, " has ", n, " problems:\n",
             paste0("  * ", gsub("\n", "\n  ", e$messages, fixed = TRUE),
                    collapse = "\n")),
      problems = e$messages
    )
  }
  list(add = add, try = try_check, raise = raise)
}


#' Format an identifier list for an error message
#'
#' @description
#' Lists up to ten identifiers and, when the list is longer, appends the
#' number of remaining entries.
#' 
#' @param ids Vector of identifiers to display.
#' @return One string.
#' @keywords internal
#' @noRd
show_ids <- function(ids) {
  more <- length(ids) - 10L
  paste0(toString(utils_head(ids, 10L)),
         if (more > 0L) paste0(" and ", more, " more") else "")
}


#' Check that person identifiers are valid and unique
#'
#' @description
#' Converts identifiers to numeric values for checking. Rejects missing or
#' non-finite values, fractions, values outside the supported integer range
#' and duplicates.
#'
#' @details
#' After numeric conversion, identifiers must be whole numbers with an
#' absolute value no greater than `.Machine$integer.max`. Uniqueness is
#' checked on the converted values. This helper does not change the input
#' identifiers stored elsewhere.
#'
#' @param ids The identifiers.
#' @param what Which table, for the message.
#' @return Invisibly returns `TRUE` if these checks pass. Otherwise,
#'   raises a `micsim_engine_error`.
#' @keywords internal
#' @noRd
check_identifiers <- function(ids, what) {
  numeric_ids <- suppressWarnings(as.numeric(ids))
  bad <- is.na(numeric_ids) | numeric_ids != round(numeric_ids) |
    abs(numeric_ids) > .Machine$integer.max
  bad[is.na(bad)] <- TRUE
  if (any(bad)) {
    abort_engine(paste0(
      "In the ", what, ", identifiers must be whole numbers.\n",
      "  Offending: ", toString(utils_head(ids[bad], 10L)), "."
    ))
  }
  twice <- unique(numeric_ids[duplicated(numeric_ids)])
  if (length(twice) > 0L) {
    abort_engine(paste0(
      "In the ", what, ", identifiers are repeated: ",
      toString(utils_head(twice, 10L)), "."
    ))
  }
  invisible(TRUE)
}


#' Check immigrant data before registration
#'
#' @description
#' Rejects an immigrant population that describes impossible arrivals,
#' naming every offending identifier.
#'
#' @details
#' Four conditions, all checked and reported together:
#'
#' * arriving before being born;
#' * arriving before the simulation starts;
#' * arriving after the simulation ends, which can have no effect;
#' * arriving already older than `max_age`, which means they were never at
#'   risk.
#'
#' Every offending identifier is listed. Identifier, overlap and state
#' checks are reported in the same error.
#'
#' @section Assumptions:
#' Identifiers must not collide with the starting population's. That is
#' not checked here but in the person store, which rejects a duplicate at
#' the moment it is registered.
#'
#' @param immigrants The immigrant data frame.
#' @param horizon The simulation horizon, `c(start, end)` as `yyyymmdd`.
#' @param max_age Maximum permitted age at arrival, in years.
#' @return `invisible(TRUE)`, or raises `micsim_engine_error`.
#' @keywords internal
#' @noRd
validate_immigrants <- function(immigrants, horizon, max_age, codec = NULL,
                                init_ids = NULL) {

  expected <- c("ID", "immigrDate", "birthDate", "immigrInitState")
  if (!is.data.frame(immigrants) ||
        !identical(names(immigrants), expected)) {
    abort_engine(paste0(
      "`immigrants` must be a data frame with columns ", toString(expected),
      ",\n  in that order."
    ))
  }
  if (nrow(immigrants) == 0L) {
    abort_engine(
      "`immigrants` has no rows. Pass NULL when there is no immigration."
    )
  }

  arrival_days <- getInDays(immigrants$immigrDate)
  birth_days   <- getInDays(immigrants$birthDate)
  horizon_days <- getInDays(horizon)

  # Every check runs and all faults are reported together.
  found <- new_problem_collector()
  complain <- function(bad, what) {
    if (!any(bad, na.rm = TRUE)) return(invisible(NULL))
    found$add(paste0(
      "In the immigrant population, ", what, ".\n",
      "  Affected IDs: ", toString(immigrants$ID[which(bad)]), "."
    ))
  }

  complain(is.na(arrival_days) | is.na(birth_days),
           "some dates are not valid yyyymmdd dates")
  complain(arrival_days < birth_days,
           "some individuals arrive before they are born")
  complain(arrival_days < horizon_days[1L],
           paste0("some individuals arrive before the simulation starts;\n",
                  "  anyone present at the start belongs in the starting\n",
                  "  population instead"))
  complain(arrival_days > horizon_days[2L],
           "some individuals arrive after the simulation ends")
  complain((arrival_days - birth_days) / DAYS_PER_YEAR > max_age,
           paste0("some individuals are already older than max_age (",
                  max_age, " years) when they arrive"))

  found$try(check_identifiers(immigrants$ID, "immigrant population"))
  if (!is.null(init_ids)) {
    clash <- intersect(suppressWarnings(as.numeric(immigrants$ID)),
                       suppressWarnings(as.numeric(init_ids)))
    clash <- clash[!is.na(clash)]
    if (length(clash) > 0L) {
      found$add(paste0(
        "Identifiers are used in both the starting and the immigrant\n",
        "  population: ", toString(utils_head(clash, 10L)), "."
      ))
    }
  }
  if (!is.null(codec)) {
    state_idx <- found$try(
      codec$encode(as.character(immigrants$immigrInitState)))
    if (!is.null(state_idx)) {
      found$try(refuse_absorbing_entry(codec, state_idx, immigrants$ID,
                                       "immigrant population"))
    }
  }
  found$raise("The immigrant population")

  invisible(TRUE)
}


#' Register the immigrant population
#'
#' @description
#' Adds immigrants to the person store in input row order and returns their
#' internal storage positions, called slots.
#'
#' @details
#' An immigrant differs from a member of the starting population in
#' exactly three ways, and in nothing else:
#'
#' * they enter the population on their arrival date rather than on the
#'   simulation's start date, so their first event is drawn from that
#'   moment;
#' * they carry `is_migrant = TRUE`, which the person store’s print method
#'   counts in order to left-truncate their first episode at arrival
#'   rather than at the start of the simulation;
#' * their `initState` is `immigrInitState`, the state they hold on
#'   arrival.
#'
#' Duration-clock initialisation is handled by the person store, as it is
#' for the starting population.
#'
#' **No lag arithmetic.** Because the queue stores each event's absolute
#' date, scheduling an immigrant's first event from their arrival date
#' needs no correction term: the arrival date is simply the moment we look
#' forward from.
#'
#' @param persons A `person_store` to fill.
#' @param immigrants The immigrant data frame, already validated.
#' @param codec A `state_codec`.
#' @param max_age Maximum age in years. Accepted for interface consistency;
#'   this helper does not use it. Arrival ages are checked separately by
#'   `validate_immigrants()`.
#' @param rng_states Integer matrix with one row of six per individual,
#'   from `person_stream_states()`; `NULL` leaves the streams unset
#'   (unit tests).
#'
#' @return Integer vector of slots, in row order.
#'
#' @keywords internal
#' @noRd
seed_immigrant_population <- function(persons, immigrants, codec, max_age,
                                      rng_states = NULL) {

  arrival_days <- getInDays(immigrants$immigrDate)
  birth_days   <- getInDays(immigrants$birthDate)
  state_idx    <- codec$encode(as.character(immigrants$immigrInitState))
  refuse_absorbing_entry(codec, state_idx, immigrants$ID, "immigrant population")

  slots <- integer(nrow(immigrants))
  for (i in seq_len(nrow(immigrants))) {
    slots[i] <- persons$add(
      id             = immigrants$ID[i],
      birth_days     = birth_days[i],
      state_idx      = state_idx[i],
      init_state_idx = state_idx[i],
      birth_ymd      = as.character(immigrants$birthDate[i]),
      is_migrant     = TRUE,
      rng_state      = if (is.null(rng_states)) NULL else rng_states[i, ]
    )
  }
  slots
}


#' Reject absorbing states at population entry
#'
#' @description
#' Checks that individuals enter the simulated population in transient
#' states. Absorbing states have no outgoing transitions and are not
#' accepted as initial states for starting individuals or immigrants.
#'
#' @param codec A `state_codec`.
#' @param state_idx Their initial state indices.
#' @param ids Their identifiers.
#' @param what Which table, for the message.
#' @return `invisible(TRUE)`, or raises `micsim_engine_error`.
#' @keywords internal
#' @noRd
refuse_absorbing_entry <- function(codec, state_idx, ids, what) {
  absorbed <- codec$is_absorbing(state_idx)
  if (any(absorbed)) {
    abort_engine(paste0(
      "In the ", what, ", some individuals start in an absorbing state\n",
      "  (", toString(unique(codec$decode(state_idx[absorbed]))), "). Affected IDs: ",
      toString(ids[absorbed]), "."
    ))
  }
  invisible(TRUE)
}


#' Find the transition corresponding to a scheduled event
#'
#' @description
#' Finds the registry entry connecting the individual's current state to
#' the scheduled destination. The engine uses this entry to identify which
#' duration clocks to reset and to supply transition details to hooks.
#'
#' @details
#' The queue stores an individual's slot, event date and destination state.
#' When the event is processed, this helper searches the transitions
#' available from the current state and returns the first matching
#' destination.
#'
#' @param registry A `rate_registry`.
#' @param from_idx,to_idx State indices.
#' @return The edge. Raises a `micsim_engine_error` if no matching
#'   transition exists.
#' @keywords internal
#' @noRd
find_edge <- function(registry, from_idx, to_idx) {
  for (edge in registry$feasible(from_idx)) {
    if (edge$to_idx == to_idx) return(edge)
  }
  abort_engine(paste0(
    "A scheduled event names a transition the registry does not have:\n",
    "  ", registry$codec$decode(from_idx), " -> ",
    registry$codec$decode(to_idx), ".\n",
    "  This is an internal error, not a problem with the model."
  ))
}


#' Raise an engine error
#'
#' @description
#' Stops processing with an error of class `micsim_engine_error`.
#'
#' @param problem Character string containing the error message.
#' @return Does not return; raises a `micsim_engine_error`.
#' @keywords internal
#' @noRd
abort_engine <- function(problem) {
  micsim_abort("micsim_engine_error", problem)
}


#' Print a simulation engine summary
#'
#' @param x A `micsim_engine`.
#' @param ... Additional arguments, not used. Included for compatibility
#'   with the [print()] generic.
#' @return `x`, invisibly, after printing the summary.
#' @exportS3Method print micsim_engine
#' @noRd
print.micsim_engine <- function(x, ...) {
  cat("<micsim_engine>\n")
  cat("  states     : ", x$codec$n, " (", x$codec$K, " transient, ",
      length(x$codec$absorbing), " absorbing)\n", sep = "")
  cat("  transitions: ", length(x$registry$edges), "\n", sep = "")
  cat("  run it with engine$run(), optionally engine$run(seed = 234)\n")
  invisible(x)
}
