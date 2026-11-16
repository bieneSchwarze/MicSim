# =====================================================================
# parallel.R
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   Runs one simulation across several processor cores.
#
# PLACE IN THE PIPELINE
#   micSimParallel() --(splits the population)--> micSim() on each worker
#                    --(recombines, renumbers newborns)--> the `pop` frame
#
# WHAT YOU NEED TO KNOW
#   Splitting is valid when there are no interactions between people assigned
#   to different workers. Every worker runs with the same seed, and every
#   individual draws from their own stream (rng_streams.R), so the split
#   changes nothing but the numbers newborns are given.
#   Newborn IDs are allocated locally and then renumbered above all
#   supplied IDs. Supplied IDs must already be unique across workers.
#
# GUARANTEES
#   * Every ID in the combined result is unique.
#   * The same `seeds` gives the same life courses on any number of cores;
#     only the identifiers of newborns depend on the split.
# =====================================================================


#' Run a microsimulation on several cores
#'
#' @description
#' Splits the population across processor cores, runs [micSim()] on each,
#' and recombines the results into one `pop` data frame.
#'
#' @details
#' Individuals are independent of one another in MicSim's model, so
#' simulating them in separate groups and concatenating the results is
#' exactly equivalent to simulating them together.
#'
#' **Reproducibility.** Every worker runs with the seed `seeds`, and
#' every individual draws from their own L'Ecuyer-CMRG stream derived
#' from it and their ID. The life courses are therefore the same on any
#' number of cores and under any split; only the identifiers of newborns
#' depend on the split (see Identifiers).
#'
#' **Identifiers.** Each worker assigns newborn IDs starting above the
#' largest ID in its own group, so different workers can assign the same
#' ID to different newborns. When the results are combined, newborns
#' receive consecutive IDs starting one above the larger of the highest
#' supplied ID and the total number of people in the starting and
#' immigrant populations. All supplied IDs are preserved.
#'
#' **Rate functions** are exported to the workers along with everything
#' else visible where `micSimParallel()` was called.
#'
#' @param initPop Data frame describing the starting population; see
#'   [micSim()].
#' @param immigrPop Data frame describing immigrants, or `NULL`.
#' @param initPopList,immigrPopList Pre-computed splits, one list element
#'   per core. Supply these to control which individuals go to which core;
#'   otherwise the split is made automatically into near-equal chunks.
#'   Immigrants assigned to a core without a starting population are
#'   simulated on the first core that has one.
#' @param transitionMatrix,absStates,varInitStates,initStatesProb,fixInitStates
#'   As in [micSim()].
#' @param maxAge,simHorizon,fertTr,monthSchoolEnrol As in [micSim()].
#' @param cores Number of processor cores to use.
#' @param seeds Integer seed of the run, used by every worker.
#' @param verbose Report progress through [message()].
#'
#' @return The same data frame [micSim()] returns, for the whole
#'   population.
#'
#' @seealso [micSim()], [MicSim-concepts]
#'
#' @section Deprecated:
#' `micSimParallel()` is deprecated as of MicSim 4.0.0 and will be removed in a
#' future release. It still works, and shows a message once per R session
#' naming its replacement: build the model with [micsim_model()] and run it with
#' `simulate(model, seed = , cores = )` ([simulate.micsim_model()]). Set
#' `options(MicSim.show_deprecation = FALSE)` to turn the message off.
#'
#' @keywords internal
#' @export
micSimParallel <- function(initPop          = NULL,
                           immigrPop        = NULL,
                           initPopList      = list(),
                           immigrPopList    = list(),
                           transitionMatrix,
                           absStates        = NULL,
                           varInitStates    = c(),
                           initStatesProb   = c(),
                           fixInitStates    = c(),
                           maxAge           = 99,
                           simHorizon,
                           fertTr           = c(),
                           monthSchoolEnrol = c(),
                           cores            = 1L,
                           seeds            = 1254L,
                           verbose          = TRUE) {

  deprecate_message("micSimParallel")
  caller <- parent.frame()

  # Force every argument the workers use, HERE. They are otherwise lazy
  # promises: the worker closure below would be serialised with them
  # unevaluated, and a worker would then evaluate e.g. `maxAge = my_max`
  # or `fixInitStates = f(x)` in an environment where `my_max` or `f` does
  # not exist, failing with "object not found" / "could not find function".
  force(transitionMatrix); force(absStates); force(fixInitStates)
  force(varInitStates); force(initStatesProb); force(maxAge)
  force(simHorizon); force(fertTr); force(monthSchoolEnrol)

  cores <- check_cores(cores, error_class = "micsim_parallel_error")

  # Read factor and Date columns by their values (see standardise_people()).
  initPop       <- standardise_people(initPop)
  immigrPop     <- standardise_people(immigrPop)
  initPopList   <- lapply(initPopList, standardise_people)
  immigrPopList <- lapply(immigrPopList, standardise_people)

  init_chunks <- split_population(initPop, initPopList, cores, "initPop")
  n_start <- sum(vapply(init_chunks, nrow, integer(1L)))
  if (n_start == 0L) {
    micsim_abort(
      "micsim_parallel_error",
      "No starting population was given, either as `initPop` or as `initPopList`."
    )
  }

  immigr_chunks <- if (is.null(immigrPop) && length(immigrPopList) == 0L) {
    vector("list", cores)          # all NULL: no immigration
  } else {
    split_population(immigrPop, immigrPopList, cores, "immigrPop")
  }
  immigr_chunks <- rehome_immigrants(init_chunks, immigr_chunks)
  n_immigrants <- sum(vapply(immigr_chunks,
                             function(x) if (is.null(x)) 0L else nrow(x),
                             integer(1L)))

  if (verbose) {
    message("Starting ", cores, " worker(s) ...")
  }

  cluster <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cluster), add = TRUE)

  # Every worker runs with the same seed. Each individual draws from their
  # own stream, derived from the seed and their ID (rng_streams.R), so the
  # split does not change the result.
  check_seed(seeds)
  force(seeds)

  # A worker is a fresh R process. It has never seen the caller's rate
  # functions, and micSim() resolves them by name, so everything visible
  # where micSimParallel() was called is sent across.
  parallel::clusterExport(cluster, varlist = ls(envir = caller),
                          envir = caller)

  # The rate functions themselves are resolved here, on the master, exactly
  # as micSim() resolves them (inheriting up from the caller), and exported
  # by name. ls(caller) alone misses them whenever micSimParallel() is
  # called from inside a function while the rates live in the global
  # environment, which micSim() on one core finds without complaint.
  rate_names <- unique(as.character(transitionMatrix[transitionMatrix != "0"]))
  rate_env <- new.env(parent = emptyenv())
  for (nm in rate_names) resolve_rate(nm, caller, rate_env)
  parallel::clusterExport(cluster, varlist = ls(envir = rate_env),
                          envir = rate_env)

  work <- lapply(seq_len(cores), function(k) {
    list(initPop = init_chunks[[k]], immigrPop = immigr_chunks[[k]])
  })

  results <- parallel::parLapply(cluster, work, function(chunk) {
    if (is.null(chunk$initPop) || nrow(chunk$initPop) == 0L) return(NULL)
    # The internal body of micSim(), so that the deprecation message is
    # not shown on every worker. The worker's own frame is where the rate
    # functions exported above are found, as with micSim() before.
    legacy_micsim(
      initPop          = chunk$initPop,
      immigrPop        = chunk$immigrPop,
      transitionMatrix = transitionMatrix,
      fixInitStates    = fixInitStates,
      varInitStates    = varInitStates,
      initStatesProb   = initStatesProb,
      maxAge           = maxAge,
      simHorizon       = simHorizon,
      fertTr           = fertTr,
      monthSchoolEnrol = monthSchoolEnrol,
      verbose          = FALSE,
      caller           = environment(),
      seed             = seeds
    )
  })

  if (verbose) message("Recombining ...")
  pop <- recombine_populations(results, work, n_start + n_immigrants)
  if (verbose) message("Simulation has finished.")

  # What micSim() records for the converters, for the combined result:
  # rbind() keeps the attributes of none of the pieces.
  spec <- attr(results[[which(!vapply(results, is.null, logical(1L)))[1L]]],
               "micsim_spec")
  if (!is.null(spec)) {
    spec$immigrPop <- if (n_immigrants > 0L) {
      do.call(rbind, Filter(Negate(is.null), immigr_chunks))
    } else {
      NULL
    }
    attr(pop, "micsim_spec") <- spec
  }
  pop
}


#' Rehome immigrants to a starting population
#'
#' @description
#' Give immigrants of a worker without a starting population to the first
#' worker that has one. A run needs a starting population, and since every
#' individual draws from their own stream, which worker simulates an
#' immigrant does not matter.
#'
#' @param init_chunks List of starting populations, one per worker.
#' @param immigr_chunks List of immigrant populations, one per worker.
#'
#' @return The updated list of immigrant populations, with immigrants of
#'   a worker without a starting population given to the first worker that
#'   has one.
#'
#' @keywords internal
#' @noRd
rehome_immigrants <- function(init_chunks, immigr_chunks) {
  has_start <- vapply(init_chunks, function(x) NROW(x) > 0L, logical(1L))
  orphaned <- which(!has_start &
                      vapply(immigr_chunks, function(x) NROW(x) > 0L, logical(1L)))
  if (length(orphaned) == 0L || !any(has_start)) return(immigr_chunks)
  home <- which(has_start)[1L]
  immigr_chunks[[home]] <- do.call(rbind, c(immigr_chunks[home],
                                            immigr_chunks[orphaned]))
  for (k in orphaned) immigr_chunks[k] <- list(NULL)
  immigr_chunks
}


#' Split a population into one chunk per core
#'
#' @description
#' Checks that a supplied split has one element per worker and returns it.
#' Otherwise, distributes the population rows among workers in turn to
#' produce groups of similar size.
#'
#' @details
#' Rows are assigned round-robin. The returned list has exactly `cores`
#' elements, with zero-row data frames for groups receiving no rows.
#' Dealing rows round-robin rather than in contiguous blocks keeps a
#' population sorted by age or state from sending all the old people to
#' one worker and all the young to another.
#'
#' A caller-supplied split must have exactly `cores` elements; those
#' elements are then returned unchanged. Duplicate IDs across chunks
#' are not checked here or by the separate person stores on workers.
#' Check supplied identifiers before splitting by hand.
#'
#' @param pop The population to split, or `NULL`.
#' @param given A pre-computed split, or an empty list.
#' @param cores Number of chunks wanted.
#' @param what Name of the argument, for error messages.
#' @return A list of `cores` data frames, some possibly with no rows.
#' @keywords internal
#' @noRd
split_population <- function(pop, given, cores, what) {

  if (length(given) > 0L) {
    if (length(given) != cores) {
      micsim_abort(
        "micsim_parallel_error",
        paste0("`", what, "List` has ", length(given), " element(s) but ",
               cores, " core(s) were asked for.")
      )
    }
    return(given)
  }

  if (is.null(pop) || nrow(pop) == 0L) {
    return(vector("list", cores))
  }

  # Deal the rows round-robin rather than in blocks, so that a population
  # sorted by age or by state does not give one worker all the old people
  # and another all the young. 
  assignment <- rep_len(seq_len(cores), nrow(pop))
  chunks <- split(seq_len(nrow(pop)), assignment)

  lapply(seq_len(cores), function(k) {
    rows <- chunks[[as.character(k)]]
    if (is.null(rows)) pop[0L, , drop = FALSE] else pop[rows, , drop = FALSE]
  })
}


#' Combine the workers' results and renumber the newborns
#'
#' @description
#' Concatenates each worker's `pop`, giving every individual born during
#' the simulation a fresh identifier so that no two collide.
#'
#' @details
#' Combines workers' life courses. An ID in a worker's output that did not
#' appear in its inputs identifies a newborn. Newborn IDs are reassigned
#' starting one above the larger of `n_supplied` and the largest supplied ID.
#' The same mapping updates `motherID` when the mother was herself born
#' during the simulation.
#'
#' @param results One `pop` data frame per worker, some possibly `NULL`.
#' @param work The inputs each worker was given, to tell newborns from
#'   supplied individuals.
#' @param n_supplied Total number of identifiers the user supplied.
#' @param id_maps If `TRUE`, also return the renumbering, so that the
#'   identifiers in the workers' subsystem reports can be renumbered the
#'   same way.
#' @param newborns List with one element per worker holding `ID` and `birth_days`
#'   of its newborns, or `NULL`. When given, all newborns are numbered together in
#'   order of birth from one above the largest supplied `ID` (or 1), as in a
#'   single-core run, and `n_supplied` is not used.
#'
#' @return A combined data frame sorted by numeric ID. If `id_maps = TRUE`,
#'   returns a list with `pop` and `id_maps`. Each worker's entry in `id_maps`
#'   contains `old` and `new` ID vectors, or is `NULL` if no IDs changed.
#'
#' @keywords internal
#' @noRd
recombine_populations <- function(results, work, n_supplied,
                                  id_maps = FALSE, newborns = NULL) {

  supplied_ids <- unlist(lapply(work, function(w) c(w$initPop$ID, w$immigrPop$ID)),
                         use.names = FALSE)
  next_id <- max(c(n_supplied, supplied_ids), na.rm = TRUE)
  renumbered <- list()
  maps <- vector("list", length(results))

  # With the newborns' exact birth times (simulate() on several cores),
  # all newborns are numbered together in order of birth, from one above
  # the largest supplied ID (or 1): the numbers a single-core run gives
  # them. Simultaneous births are ordered by worker, then by the worker's
  # own numbering.
  if (!is.null(newborns)) {
    born <- do.call(rbind, lapply(seq_along(newborns), function(k) {
      b <- newborns[[k]]
      if (is.null(b) || length(b$ID) == 0L) return(NULL)
      data.frame(worker = k, old = b$ID, birth_days = b$birth_days)
    }))
    if (!is.null(born)) {
      born <- born[order(born$birth_days, born$worker, born$old), , drop = FALSE]
      first <- max(c(0, supplied_ids), na.rm = TRUE)
      born$new <- as.integer(first + seq_len(nrow(born)))
      maps <- lapply(seq_along(results), function(k) {
        b <- born[born$worker == k, , drop = FALSE]
        if (nrow(b) == 0L) NULL else list(old = b$old, new = b$new)
      })
    }
  }

  for (k in seq_along(results)) {
    piece <- results[[k]]
    if (is.null(piece) || nrow(piece) == 0L) next

    supplied <- c(work[[k]]$initPop$ID, work[[k]]$immigrPop$ID)
    present <- unique(piece$ID)

    if (!is.null(newborns)) {
      map <- maps[[k]]
    } else {
      born_here <- setdiff(present, supplied)
      map <- NULL
      if (length(born_here) > 0L) {
        map <- list(old = born_here, new = next_id + seq_along(born_here))
        next_id <- next_id + length(born_here)
        maps[[k]] <- map
      }
    }

    if (!is.null(map)) {
      newborns_k <- map$old
      replacement <- map$new
      hit <- match(piece$ID, newborns_k)
      piece$ID[!is.na(hit)] <- replacement[hit[!is.na(hit)]]
      # A newborn's mother keeps her own identifier, but a newborn who
      # went on to have a child of her own is somebody's mother too.
      if ("motherID" %in% names(piece)) {
        hit <- match(piece$motherID, newborns_k)
        piece$motherID[!is.na(hit)] <- replacement[hit[!is.na(hit)]]
      }
    }

    renumbered[[length(renumbered) + 1L]] <- piece
  }

  if (length(renumbered) == 0L) {
    micsim_abort("micsim_parallel_error",
                 "Every worker returned an empty result.")
  }

  pop <- do.call(rbind, renumbered)
  pop <- pop[order(as.numeric(pop$ID)), , drop = FALSE]
  row.names(pop) <- NULL
  if (id_maps) list(pop = pop, id_maps = maps) else pop
}
