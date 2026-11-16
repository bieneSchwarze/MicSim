# =====================================================================
# simulate.R
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   Runs a model, and holds the result together with the model that made
#   it.
#
# PLACE IN THE PIPELINE
#   micsim_model --> simulate() --> micsim_result --> as_long() / as_wide()
#
# WHAT YOU NEED TO KNOW
#   Nothing new about the method. The simulation itself is `engine.R`; this
#   is the object the user holds afterwards.
#
# WHY THE RESULT IS AN OBJECT AND NOT A DATA FRAME
#   A `pop` data frame as micSim() used to return it loses the model
#   that produced it, the seed it used and what happened during the run.
#   Keeping them alongside `pop` means a result can be converted without
#   repeating the state space, horizon and immigrants, can be traced back
#   to the model that made it, and shows whether the numerical integration
#   degraded. `as.data.frame()` still returns the plain data frame.
#
# GUARANTEES
#   * `as.data.frame(result)` is identical to what the old micSim() returns
#     for the same model and seed. This must not change.
#   * Converting a result needs no arguments: everything the converters ask
#     for is in the model the result carries.
# =====================================================================


#' Simulate from a model
#'
#' @description
#' Runs the simulation a [micsim_model()] describes and returns the result
#' together with the model, the seed and what happened during the run.
#'
#' @details
#' This is a method for [stats::simulate()], so it has that generic's
#' argument names. With the default `seed = NULL`, a seed is drawn from R's
#' current random number state (one draw), so `set.seed(234)` before the
#' call makes a run reproducible. The drawn seed is stored in the result.
#' Apart from that draw, the caller's random number state is left exactly
#' as it was, including the generator kind. Every individual in the simulation
#' draws from their own random-number stream, derived from the seed and
#' their ID (a newborn's from its mother's). A life course therefore does not
#' depend on who else is in the population, on the order of the population's rows
#' or of the transitions, or on the number of cores. Two scenarios simulated with
#' the same seed share their random numbers (common random numbers), i.e.
#' an individual whose rates are the same in both lives the same life in both.
#'
#' A model can be simulated repeatedly. Nothing is shared between runs:
#' every subsystem's state is created fresh at the start of each run. If
#' the model was built by a different MicSim version, a `micsim_version_warning`
#' is given.
#'
#' After the run, every transition rate is estimated from the simulated
#' population and compared with the supplied one ([reestimate()]).
#' Using `plot()` on the result shows the comparison.
#' Pass `reestimate = FALSE` to skip it.
#'
#' @param object A `micsim_model`.
#' @param nsim Number of simulations. Only `1` is supported; run the model
#'   again, with a different seed, for more.
#' @param seed Integer seed, or `NULL` to draw one from R's current random number
#' state. The same seed gives the same result on any number of cores.
#' @param cores Number of processor cores. Splitting the population across
#'   cores is exact rather than approximate, because individuals in
#'   MicSim's model are independent of one another. The result does not depend
#'   on the number of cores. No more workers are started than there are
#'   individuals in the starting population.
#' @param verbose Report progress through [message()].
#' @param reestimate Estimate the transition rates from the simulated
#'   population after the run and compare them with the supplied ones
#'   ([reestimate()]).
#' @param ... Not used.
#'
#' @return An object of class `micsim_result`, a list with
#'   \describe{
#'     \item{`pop`}{the simulated population:
#'       one row per transition, one row with `NA` transition
#'       columns for an individual who never moved. [as.data.frame()]
#'       returns it.}
#'     \item{`model`}{the `micsim_model` that was simulated, so the result
#'       is self-describing and [as_long()] needs no arguments.}
#'     \item{`seed`}{the seed the run used; when none was supplied, the one
#'       drawn from R's random number state.}
#'     \item{`diagnostics`}{`n_events`, `n_individuals`, `n_starting`,
#'       `n_immigrants` and `n_newborns`; `n_forced` (transitions forced
#'       by an `Inf` rate), `n_unreliable` (waiting times drawn from a
#'       degraded numerical integration); `peak_queue` (largest number
#'       of pending events; for a multi-core run the largest of any core),
#'       `cores` (worker processes used) and `seconds` (wall time of
#'       the simulation).}
#'     \item{`extensions`}{a named list with what each subsystem reports,
#'       e.g. `result$extensions$fertility$births`; see [fertility()] and
#'       [school_enrolment()].}
#'     \item{`reestimation`}{a `micsim_reestimation` ([reestimate()]), or
#'       `NULL` with `reestimate = FALSE`.}
#'     \item{`meta`}{provenance: `micsim_version`, `r_version`, `platform`,
#'       `rng_kind` (always `"L'Ecuyer-CMRG"`), `seed_drawn`
#'       (`TRUE` if the seed was drawn from R's random number state),
#'       `created` (time the run finished) and `reestimation_seconds`.}
#'   }
#'
#' @seealso [micsim_model()], [as_long()], [as_wide()], [MicSim-concepts]
#'
#' @examples
#' mortality <- function(age, calTime) 0.00003 * exp(0.1 * age)
#' model <- micsim_model(
#'   state_space = state_space(sex = c("m", "f")),
#'   transitions = transitions(transition(NA, "dead", mortality)),
#'   population  = data.frame(ID = 1:10, birthDate = rep("19500101", 10),
#'                            initState = rep("f", 10),
#'                            stringsAsFactors = FALSE),
#'   horizon = c(20000101, 20401231), max_age = 120
#' )
#' result <- simulate(model, seed = 234)
#' result
#'
#' @export
simulate.micsim_model <- function(object, nsim = 1, seed = NULL,
                                  cores = 1L, verbose = FALSE,
                                  reestimate = TRUE, ...) {

  if (!identical(as.integer(nsim), 1L)) {
    micsim_abort(
      "micsim_model_error",
      paste0(
        "`nsim` must be 1. To run a model several times, call simulate()\n",
        "  again with a different seed. Each result carries its own seed."
      )
    )
  }

  check_model_version(object)
  check_seed(seed)
  cores <- check_cores(cores)

  if (cores > 1L) {
    out <- simulate_distributed(object, seed, cores, verbose)
  } else {
    run <- micsim_engine_for(object, verbose = verbose)$run(seed = seed)
    d <- run$diagnostics
    out <- new_micsim_result(
      pop = run$pop, model = object, seed = d$seed, rng_kind = d$rng_kind,
      counts = d, cores = 1L, seconds = d$seconds,
      extensions = run$extensions, seed_drawn = d$seed_drawn
    )
  }
  if (isTRUE(reestimate)) {
    if (verbose) message("Reestimating the transition rates ...")
    t0 <- Sys.time()
    out$reestimation <- reestimate_rates(out$pop, object)
    out$meta$reestimation_seconds <-
      as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }
  out
}


#' Warn when a model was built by a different MicSim version.
#'
#' @description
#' A model stores internal objects (state codec, rate registry, compiled
#' subsystems), not only the user's inputs, so a model saved with
#' saveRDS() under one version may not run correctly under another. The
#' check is a warning, not an error: most models will still run.
#'
#' @param model A `micsim_model`.
#'
#' @return `TRUE` when the versions agree, otherwise
#' `FALSE` after a `micsim_version_warning`, invisibly.
#'
#' @keywords internal
#' @noRd
check_model_version <- function(model) {
  built <- model$micsim_version
  running <- as.character(getNamespaceVersion("MicSim"))
  if (identical(built, running)) return(invisible(TRUE))
  warning(micsim_warning(
    "micsim_version_warning",
    paste0(
      "This model was built with ",
      if (is.null(built)) "a MicSim version that did not record itself"
      else paste0("MicSim ", built),
      ", but MicSim ", running, " is running it.
",
      "  A model stores internal objects whose structure can change between
",
      "  versions. Rebuild it with micsim_model() (and add_subsystem())
",
      "  if the simulation fails or its results look wrong."
    ),
    built = built, running = running
  ))
  invisible(FALSE)
}


#' Check `cores`
#' @description
#' Checks if `cores` is a numeric input and if the number fulfils
#' the requirements. If not, an error is thrown.
#'
#' @param cores Number of cores.
#' @param error_class Condition class of the error. `micSimParallel()`
#'   passes `micsim_parallel_error`.
#'
#' @return Integer number of cores.
#'
#' @keywords internal
#' @noRd
check_cores <- function(cores, error_class = "micsim_model_error") {
  n <- suppressWarnings(as.numeric(cores))
  if (length(n) != 1L || is.na(n) || !is.finite(n) || n != round(n) ||
        n < 1 || n > .Machine$integer.max) {
    micsim_abort(error_class,
                 "`cores` must be a single whole number, at least one.")
  }
  as.integer(n)
}


#' Assemble a simulation result
#'
#' @description
#' The only constructor of `micsim_result`. Used by both single-core and
#' multi-core runs, so that the two always have the same fields.
#'
#' @param pop,model,seed,extensions The parts of the result.
#' @param rng_kind The kind of random number generator used.
#' @param counts A list with `n_unreliable`, `n_forced` and `peak_queue`.
#' @param cores Number of cores.
#' @param seconds Wall time of the simulation.
#' @param seed_drawn TRUE if the seed was drawn from the caller's random stream.
#'
#' @return A `micsim_result`; `reestimation` is filled in by the caller.
#' @keywords internal
#' @noRd
new_micsim_result <- function(pop, model, seed, rng_kind, counts, cores,
                              seconds, extensions, seed_drawn = FALSE) {
  people <- pop[!duplicated(pop$ID), , drop = FALSE]
  n_starting <- sum(people$ID %in% model$population$ID)
  n_immigrants <- if (is.null(model$immigrants)) 0L else
    sum(people$ID %in% model$immigrants$ID)
  structure(
    list(
      pop         = pop,
      model       = model,
      seed        = if (is.null(seed)) NA_integer_ else as.integer(seed),
      diagnostics = list(
        n_events          = sum(!is.na(pop$From)),
        n_individuals     = nrow(people),
        n_starting        = as.integer(n_starting),
        n_immigrants      = as.integer(n_immigrants),
        n_newborns        = as.integer(nrow(people) - n_starting - n_immigrants),
        n_forced          = as.integer(counts$n_forced),
        n_unreliable      = as.integer(counts$n_unreliable),
        peak_queue        = as.integer(counts$peak_queue),
        cores             = as.integer(cores),
        seconds           = as.numeric(seconds)
      ),
      extensions   = extensions,
      reestimation = NULL,
      meta = list(
        micsim_version = as.character(getNamespaceVersion("MicSim")),
        r_version      = R.version.string,
        platform       = R.version$platform,
        rng_kind       = rng_kind,
        # TRUE when no seed was supplied and the run drew `seed` from the
        # caller's random stream.
        seed_drawn     = isTRUE(seed_drawn),
        created        = Sys.time(),
        reestimation_seconds = NA_real_
      )
    ),
    class = "micsim_result"
  )
}


#' Run a model across several cores
#'
#' @description
#' Splits the population, runs the same model on each piece in a separate
#' R process, and recombines the results.
#'
#' @details
#' This is the same mechanism the old [micSimParallel()] used. 
#' Each core runs one worker, a separate R process. Individuals are dealt
#' out to the workers in turn (round-robin), every worker runs with the
#' same seed, and every individual draws from their own stream, so the split
#' changes nothing. Newborns are renumbered in order of birth, which gives
#' them the identifiers a single-core run gives them.
#'
#' Splitting is exact rather than approximate, as individuals in MicSim's
#' model do not depend on one another. This means that simulating them in groups
#' and concatenating the results is a valid draw from the same model. The result
#' is identical to a single-core run with the same seed. No more workers are
#' used than there are individuals in the starting population; immigrants are
#' dealt among the same workers.
#'
#' Each subsystem's per-worker outputs are combined by its own `merge`
#' (counts, for example, are added up), and the run diagnostics
#' are summed across workers.
#'
#' @section Assumptions:
#' Not valid for an eventual linked lives model.
#'
#'
#' @param model A `micsim_model`.
#' @param seed Integer seed, or NULL to draw one from R's
#'   random number state, as in a single-core run.
#' @param cores Number of cores.
#' @param verbose Report progress.
#' @return A `micsim_result`.
#' @keywords internal
#' @noRd
simulate_distributed <- function(model, seed, cores, verbose) {

  unsafe <- names(Filter(function(s) !isTRUE(s$parallel_safe),
                         model$subsystems))
  if (length(unsafe) > 0L) {
    micsim_abort(
      "micsim_model_error",
      paste0(
        "Subsystem(s) ", toString(paste0("`", unsafe, "`")), " cannot be\n",
        "  split across cores. Splitting is only valid while individuals are\n",
        "  independent, and these subsystems link them. Use cores = 1."
      )
    )
  }

  # Settle the seed exactly as a single-core run does: without one, it is
  # drawn from the caller's stream. Every worker then runs with this seed,
  # and since each individual's random numbers depend only on the seed and
  # their identifier (rng_streams.R), the split does not change anything.
  rng <- open_run_rng(seed)
  rng$restore()

  # A worker needs at least one member of the starting population, so
  # there are never more workers than starting individuals. Immigrants are
  # dealt among the same workers.
  cores <- min(cores, nrow(model$population))
  init_chunks <- split_population(model$population, list(), cores, "population")
  immigr_chunks <- if (is.null(model$immigrants)) {
    vector("list", cores)
  } else {
    split_population(model$immigrants, list(), cores, "immigrants")
  }

  started <- Sys.time()
  if (verbose) message("Starting ", cores, " worker(s) ...")
  cluster <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cluster), add = TRUE)

  work <- lapply(seq_len(cores), function(k) {
    piece <- model
    piece$population <- init_chunks[[k]]
    piece$immigrants <- if (NROW(immigr_chunks[[k]]) > 0L) immigr_chunks[[k]]
    piece
  })

  results <- parallel::parLapply(cluster, work, run_distributed_piece,
                                 seed = rng$seed)

  if (verbose) message("Recombining ...")
  inputs <- lapply(work, function(piece) {
    list(initPop = piece$population, immigrPop = piece$immigrants)
  })
  combined <- recombine_populations(
    lapply(results, `[[`, "pop"), inputs,
    n_supplied = 0L, id_maps = TRUE,
    newborns = lapply(results, `[[`, "newborns")
  )
  pop <- combined$pop
  # Subsystem tables name people; renumber them exactly as the population.
  worker_ext <- lapply(seq_along(results), function(k) {
    renumber_extensions(results[[k]]$extensions, combined$id_maps[[k]])
  })

  total <- function(field) {
    sum(vapply(results, function(r) as.numeric(r$diagnostics[[field]]),
               numeric(1L)))
  }

  new_micsim_result(
    pop = pop, model = model, seed = rng$seed, rng_kind = "L'Ecuyer-CMRG",
    counts = list(
      n_unreliable      = total("n_unreliable"),
      n_forced          = total("n_forced"),
      peak_queue        = max(vapply(results, function(r)
        as.numeric(r$diagnostics$peak_queue), numeric(1L)))
    ),
    cores = cores,
    seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    extensions = merge_extensions(model$subsystems, worker_ext),
    seed_drawn = rng$drawn
  )
}


#' Simulate one worker's share of a multi-core run.
#' @description
#' Runs in a worker process. Being a function of the package namespace, it
#' is sent to the worker without the environment of the caller, so only
#' `piece` and `seed` travel; the worker loads MicSim to unserialise it.
#' Returns the population, the newborns' exact birth times (to number them
#' as a single-core run would), the diagnostics and the subsystem reports.
#'
#' @param piece A chunk of a `micsim_model`.
#' @param seed Integer seed of the run, settled on the main process.
#'
#' @return A list with the population, the newborns' exact birth times,
#' the diagnostics and the subsystem reports.
#'
#' @keywords internal
#' @noRd
run_distributed_piece <- function(piece, seed) {
  run <- micsim_engine_for(piece, verbose = FALSE)$run(seed = seed)
  list(pop = run$pop, newborns = run$newborns,
       diagnostics = run$diagnostics, extensions = run$extensions)
}


#' Convert a result to long or wide format
#'
#' @description
#' `as_long()` gives one row per episode, which is the layout event-history and
#' survival analysis expect. `as_wide()` gives one row per individual.
#'
#' @details
#' Neither needs arguments. The result already carries the model, and the model
#' knows the state space, the horizon, the absorbing states and the
#' immigrants, which is everything these two converters need.
#'
#' Both take the result as their first argument, so they can be put in
#' a pipe easily: `simulate(model, seed = 234) |> as_long()`.
#'
#' @param x A `micsim_result`.
#' @param ... Passed to methods.
#'
#' @return A data frame.
#'
#'   `as_long()` gives one row per episode, with the columns `ID`,
#'   `birthDate`, `Tstart`, `Tstop`, `statusEntry`, `statusExit`, `OD`,
#'   `ns`, `Episode`, and one column per state variable (a state such as
#'   `"f/DE/1"` is split into columns for sex, nationality and parity).
#'   Dates use `yyyymmdd` format. `statusEntry` is 1 if the start of the
#'   episode was observed and 0 if it was left truncated; `statusExit` is 1
#'   if its end was observed and 0 if it was right censored. `OD` says how
#'   the episode ended: the name of an absorbing state such as `"dead"`,
#'   `"cens"` if it was still running when observation ended, or the
#'   values that changed, such as `"0->1"`. `ns` is the number of episodes
#'   of the individual, and `Episode` numbers them in order.
#'
#'   `as_wide()` gives one row per individual, with the columns `ID`,
#'   `birthDate`, `initState`, `ns` (the number of transitions, 0 for
#'   someone who never moved) and numbered groups `From.1`, `To.1`,
#'   `transitionTime.1`, `transitionAge.1`, `From.2`, and so on. `motherID`
#'   is not included; join it from `as.data.frame(result)` by `ID` if you
#'   need it.
#'
#' @seealso [simulate.micsim_model()], [MicSim-concepts]
#'
#' @export
as_long <- function(x, ...) {
  UseMethod("as_long")
}


#' @rdname as_long
#' @param migr Whether immigrants' first episodes begin on their arrival date.
#'    Defaults to `TRUE` whenever the model includes immigrants.
#' @export
as_long.micsim_result <- function(x, migr = !is.null(x$model$immigrants),
                                  ...) {
  convert_to_long(
    pop            = x$pop,
    migr           = migr,
    stateSpace     = as_state_space_frame(x$model$state_space),
    absStates      = x$model$codec$absorbing,
    allTransitions = declared_transition_labels(x$model),
    simHorizon     = x$model$horizon,
    immigrPop      = x$model$immigrants,
    maxAge         = x$model$max_age
  )
}


#' @rdname as_long
#' @export
as_wide <- function(x, ...) {
  UseMethod("as_wide")
}


#' @rdname as_long
#' @export
as_wide.micsim_result <- function(x, ...) {
  convert_to_wide(x$pop)
}


#' The simulated population, as a plain data frame
#'
#' @description
#' Returns the `pop` data frame of a result. Includes one row per transition,
#' plus one row with `NA` transition columns for each individual who never
#' experienced a transition.
#'
#' @details
#' Use this wherever code still expects a plain `pop` data frame, such as an
#' existing script, an old workflow or another package. It is the same
#' format MicSim has always returned.
#'
#' @param x A `micsim_result`.
#' @param ... Not used.
#' @return The `pop` data frame.
#' @export
as.data.frame.micsim_result <- function(x, ...) {
  x$pop
}


#' Print a simulation result
#'
#' @param x A `micsim_result`.
#' @param ... Not used; present for compatibility with the [print()] generic.
#' @return `x`, invisibly. Called for the side effect of printing.
#' @export
print.micsim_result <- function(x, ...) {
  d <- x$diagnostics
  big <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)
  cat("<micsim_result>", if (!is.null(x$meta)) paste0(" MicSim ", x$meta$micsim_version), "\n", sep = "")
  cat("  ", big(d$n_events), " event(s) among ", big(d$n_individuals),
      " individual(s)", sep = "")
  parts <- c(if (isTRUE(d$n_immigrants > 0L)) paste(big(d$n_immigrants), "immigrant(s)"),
             if (isTRUE(d$n_newborns > 0L)) paste(big(d$n_newborns), "born during the run"))
  if (length(parts) > 0L) cat(" (", paste(parts, collapse = ", "), ")", sep = "")
  cat("\n")
  cat("  horizon: ", x$model$horizon[1L], " to ", x$model$horizon[2L],
      ", max age ", x$model$max_age, "\n", sep = "")
  cat("  seed   : ",
      if (is.na(x$seed)) "none" else x$seed,
      if (isTRUE(x$meta$seed_drawn)) " (drawn from the caller's random stream)",
      if (isTRUE(d$cores > 1L)) paste0(", ", d$cores, " cores"), "\n", sep = "")
  if (isTRUE(d$n_unreliable > 0L)) {
    cat("  NOTE: ", big(d$n_unreliable),
        " waiting time(s) were drawn from a degraded numerical\n",
        "        integration, typically at a jump of a step-shaped rate.\n",
        "        Compare the supplied and recovered rates with plot().\n", sep = "")
  }
  if (length(x$extensions) > 0L) {
    cat("  extensions  : ", toString(names(x$extensions)),
        " (see $extensions)\n", sep = "")
  }
  if (!is.null(x$reestimation)) {
    r <- x$reestimation
    ratios <- vapply(r, `[[`, numeric(1L), "ratio")
    cat("  reestimation: ", length(r), " transition(s)",
        if (any(!is.na(ratios))) paste0(", observed/expected ",
                                        sprintf("%.2f", min(ratios, na.rm = TRUE)), " to ",
                                        sprintf("%.2f", max(ratios, na.rm = TRUE))),
        " (see $reestimation, plot())\n", sep = "")
  }
  cat("  as.data.frame() for the population, as_long() for episodes\n")
  invisible(x)
}


#' Summarise a simulation result
#'
#' @description
#' Counts what happened: how many of each transition, the ages at which
#' transitions took place separately for each dimension of the
#' state space (each state variable that changes and each absorbing
#' state entered), and how the population was distributed across states
#' at the end.
#'
#' @param object A `micsim_result`.
#' @param ... Not used.
#' @return A list of class `summary.micsim_result`, printed by its own
#'   method.
#' @export
summary.micsim_result <- function(object, ...) {
  pop <- object$pop
  moved <- pop[!is.na(pop$From), , drop = FALSE]

  transitions_made <- if (nrow(moved) == 0L) {
    integer(0)
  } else {
    sort(table(paste(moved$From, "->", moved$To)), decreasing = TRUE)
  }

  # The state each individual was in when observation ended: their last
  # destination, or their initial state if they never moved.
  last_row <- !duplicated(pop$ID, fromLast = TRUE)
  final <- pop[last_row, , drop = FALSE]
  final_state <- ifelse(is.na(final$To), final$initState, final$To)

  structure(
    list(
      n_individuals = object$diagnostics$n_individuals,
      n_events      = object$diagnostics$n_events,
      transitions   = transitions_made,
      final_states  = sort(table(final_state), decreasing = TRUE),
      ages          = ages_by_dimension(moved, object$model$codec),
      diagnostics   = object$diagnostics
    ),
    class = "summary.micsim_result"
  )
}


#' Ages at transition, summarised separately for each dimension of the
#' state space.
#'
#' @description
#' A dimension is a state variable or an absorbing state. A transition
#' belongs to every state variable whose value it changes, regardless of
#' the other variables: all births (0 -> 1, 1 -> 2, ...) are summarised
#' together under the parity variable, whoever gave birth. A transition into
#' an absorbing state belongs to that absorbing state ("exit: dead").
#' Variables that never change and absorbing states never entered have no
#' row.
#'
#' @param moved The rows of `pop` with a transition.
#' @param codec The model's codec.
#' @return A matrix with one row per dimension and the columns n, Min.,
#' 1st Qu., Median, Mean, 3rd Qu., Max.; NULL when nothing moved.
#'
#' @keywords internal
#' @noRd
ages_by_dimension <- function(moved, codec) {
  if (nrow(moved) == 0L) return(NULL)
  ages <- moved$transitionAge
  to <- as.character(moved$To)
  absorbed <- to %in% codec$absorbing
  k <- length(codec$vars)
  from_parts <- split_labels(as.character(moved$From), k)
  to_parts   <- split_labels(ifelse(absorbed, NA_character_, to), k)

  rows <- list()
  for (v in seq_len(k)) {
    changed <- !absorbed & !is.na(from_parts[, v]) & !is.na(to_parts[, v]) &
      from_parts[, v] != to_parts[, v]
    if (any(changed)) rows[[codec$vars[v]]] <- ages[changed]
  }
  for (a in codec$absorbing) {
    hit <- to == a
    if (any(hit)) rows[[paste0("exit: ", a)]] <- ages[hit]
  }
  if (length(rows) == 0L) return(NULL)
  t(vapply(rows, function(x) {
    x <- x[!is.na(x)]
    c(n = length(x), Min. = min(x), `1st Qu.` = stats::quantile(x, 0.25, names = FALSE),
      Median = stats::median(x), Mean = mean(x),
      `3rd Qu.` = stats::quantile(x, 0.75, names = FALSE), Max. = max(x))
  }, numeric(7L)))
}


#' Print a result summary
#'
#' @param x A `summary.micsim_result`.
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.summary.micsim_result <- function(x, ...) {
  cat("Simulated ", x$n_individuals, " individual(s), ", x$n_events,
      " event(s)\n\n", sep = "")

  if (length(x$transitions) > 0L) {
    cat("Transitions:\n")
    print(x$transitions)
    cat("\nAge at transition, by dimension of the state space:\n")
    print(x$ages, digits = 4L)
    cat("\n")
  }

  cat("State at the end of observation:\n")
  print(x$final_states)
  invisible(x)
}


#' Prepare the state space for long-format conversion
#'
#' Converts a [state_space()] object into a data frame containing all
#' combinations of its state-variable values, using `expand.grid()`.
#' "convertToLongFormat()" uses the column names to name the state
#' variables in its output.
#'
#' @param ss A `micsim_state_space`.
#' @return A data frame with one column per state variable.
#' @keywords internal
#' @noRd
as_state_space_frame <- function(ss) {
  do.call(expand.grid, c(unclass(ss), list(stringsAsFactors = FALSE)))
}


#' Prepare declared transition labels for long-format conversion
#'
#' Collects transition labels from the model in the two-column format
#' expected by "convertToLongFormat()". Transitions to absorbing states
#' are excluded because the converter names them separately.
#'
#' Labels from [move()] describe a change in one state variable, such as
#' `"none->primary"`. The converter can use these labels when naming
#' transitions in its output. Full state labels are also included, but
#' labels containing several state variables do not match the
#' converter's comparison of individual variable values.
#'
#' @param model A `micsim_model`.
#' @return A two-column character matrix, or `NULL` when the model has no
#'   transient transitions.
#' @keywords internal
#' @noRd
declared_transition_labels <- function(model) {
  labels <- unlist(lapply(unclass(model$transitions), function(tr) {
    if (!is.null(tr$move)) return(paste0(tr$move$from, "->", tr$move$to))
    if (tr$to %in% model$codec$absorbing) return(NULL)
    paste0(paste(tr$from, collapse = "|"), "->", tr$to)
  }), use.names = FALSE)
  labels <- unique(labels)
  if (length(labels) == 0L) return(NULL)
  cbind(labels, "rate")
}
