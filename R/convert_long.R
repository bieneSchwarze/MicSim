# =====================================================================
# convert_long.R
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   Converts simulation output into one row per episode: a continuous
#   period spent in one state. This layout records both transitions and
#   time spent at risk for event-history and survival analysis.
#
# PLACE IN THE PIPELINE
#   engine_store.R --(the `pop` data frame)--> convert_long.R --> the user
#
# WHAT YOU NEED TO KNOW  (concepts in ?MicSim-concepts)
#   episode, left truncation, right censoring, occurrence-exposure rate.
#
#   `pop` answers "what happened, and when". Long format answers a
#   different question: "who was in which state, from when to when, and how
#   did that spell end". One row per SPELL rather than per EVENT, with a
#   start date, an end date, and two flags saying whether we actually saw
#   the spell begin and end.
#
#   These flags distinguish time observed in a state from transitions actually
#   observed. For someone in the starting population who quit smoking before
#   the simulation began, the abstinence spell is left-truncated. The person
#   contributes time at risk of relapse from the start of observation, with
#   `statusEntry = 0` because the transition into abstinence was not observed.
#   If the spell is still ongoing when observation ends, it is right-censored
#   (`statusExit = 0`).
#
#     statusEntry   1 = we saw the spell begin, 0 = left truncated
#     statusExit    1 = we saw it end,          0 = right censored
#
# WORKED EXAMPLE
#   A woman born 1985, childless at the start of a 2014-2024 simulation,
#   who has a child in 2016 and is still alive at the end:
#
#     ID  Tstart    Tstop     statusEntry statusExit OD       Episode ns
#     47  20140101  20160414  0           1          0->1     1       2
#     47  20160414  20241231  1           0          cens     2       2
#
#   Her first spell is LEFT TRUNCATED: she was childless long before 2014,
#   so `Tstart` is the start of observation and `statusEntry` is 0. Her
#   second spell is RIGHT CENSORED: it was still running when we stopped 
#   looking.
#
# GUARANTEES
#   * Every individual in `pop` appears with at least one episode.
#   * Consecutive episodes for an individual share a boundary date:
#     the next episode's `Tstart` is the preceding episode's `Tstop`.
#   * Nothing is read from the global environment; everything comes from
#     arguments or from what micSim() recorded on `pop`.
# =====================================================================


#' Convert simulation output to long format
#'
#' @description
#' Reshapes the transitions returned by [micSim()] into one row per
#' episode.
#'
#' @details
#' Each row describes one uninterrupted spell in one state. Each row records
#' its start and end dates, the transition that ended it, and flags
#' indicating whether entry and exit were observed. Each state variable
#' has its own column. For example, the state `"f/DE/1"` is split into
#' separate columns for sex, nationality and parity.
#'
#' The `OD` column describes how the episode ended. A transition to an
#' absorbing state is labelled with that state's name, such as `"dead"`.
#' An episode still in progress when observation ends is labelled
#' `"cens"`. Other transitions are labelled with the values that changed.
#' For example, `"f/DE/0->f/DE/1"` becomes `"0->1"`.
#'
#' Supplying `allTransitions` lets a SINGLE-VARIABLE pattern be reported
#' under the user's own name. It cannot do more than that: the match is
#' made component by component, so a pattern naming two or more variables
#' can never equal one of the candidates and the changed-variable form is
#' used instead.
#'
#' The output can be used to select episodes in a relevant origin states and
#' count observed exits (`statusExit == 1`) to the destination of interest. 
#' Episode lengths in can also be calculated `Tstop - Tstart`.
#' 
#' The population table is the first argument, so the function supports
#' the native R pipe:
#' `pop |> convertToLongFormat(stateSpace = stateSpace)`.
#'
#' @section Episode boundaries and observation flags:
#' **Left truncation.** For individuals born before the simulation
#' starts, the first episode begins at the simulation start and has
#' `statusEntry = 0`: entry into their initial state was not observed
#' during the simulation. For individuals born on or after the start
#' date, the first episode begins at birth and has `statusEntry = 1`.
#'
#' **Immigrants.** With `migr = TRUE`, an immigrant's first episode begins
#' on their arrival date rather than at the start of observation, and
#' counts as observed. Their earlier life happened elsewhere and is
#' outside the study.
#'
#' **Right censoring.** Anyone who has neither entered an absorbing state
#' nor already got a censored row gets one more episode, running from
#' their last transition to the end of observation: the end of the
#' simulation, or the day they reach `maxAge` if that comes first. The
#' simulation stops putting a person at risk at `maxAge - 0.01 years`, 
#'
#' @param pop The population data frame returned by [micSim()].
#' @param migr Logical. If `TRUE`, start immigrants' first episodes at
#'   arrival and mark their entries as observed. Requires `immigrPop`,
#'   either supplied explicitly or recorded on `pop`.
#' @param stateSpace The state space, as passed to
#'   [buildTransitionMatrix()]. Only its column NAMES are used, to name the
#'   state-variable columns of the result. When omitted, the names
#'   recorded by [micSim()] are used, which are `var1`, `var2`, ... because
#'   a transition matrix does not carry them.
#' @param absStates The absorbing states. Taken from `pop` when omitted.
#' @param allTransitions Optional transient-transition patterns in the
#'   format supplied to [buildTransitionMatrix()]. Only the first column,
#'   containing the patterns, is used. See Details for how patterns are
#'   matched to output labels.
#' @param simHorizon A two-element vector `c(startDate, endDate)` with
#'   dates in `yyyymmdd` format. Taken from `pop` when omitted.
#' @param immigrPop The immigrant population. Taken from `pop` when
#'   omitted. Needed only when `migr` is `TRUE`.
#' @param maxAge The maximum age of the simulation, in years. Taken from
#'   `pop` when omitted. If neither is available, final episodes end at the
#'   simulation end without an age limit.
#'
#' @return A data frame with the columns `ID`, `birthDate`, `Tstart`,
#'   `Tstop`, `statusEntry`, `statusExit`, `OD`, `ns`, `Episode`, and one
#'   column per state variable. Dates use `yyyymmdd` format. `ns` gives the total
#'   number of episodes for the individual; `Episode` numbers those
#'   episodes in start-date order. 
#'
#' @seealso [micSim()], `convertToWideFormat()`, [MicSim-concepts]
#'
#' @examples
#' stateSpace <- expand.grid(sex = c("m", "f"), fert = c("0", "1"))
#' fertRates <- function(age, calTime) ifelse(age > 15 & age < 45, 0.1, 0)
#' mortRates <- function(age, calTime) 0.00003 * exp(0.1 * age)
#'
#' tm <- buildTransitionMatrix(
#'   allTransitions = cbind("f/0->f/1", "fertRates"),
#'   absTransitions = cbind("dead", "mortRates"),
#'   stateSpace     = stateSpace
#' )
#'
#' set.seed(1)
#' pop <- micSim(
#'   initPop = data.frame(ID = 1:20, birthDate = rep("19850101", 20),
#'                        initState = rep("f/0", 20),
#'                        stringsAsFactors = FALSE),
#'   transitionMatrix = tm,
#'   varInitStates  = rbind(c("m", "0"), c("f", "0")),
#'   initStatesProb = c(0.515, 0.485),
#'   maxAge = 100, simHorizon = c(20140101, 20241231),
#'   fertTr = "f/0->f/1", verbose = FALSE
#' )
#' head(convertToLongFormat(pop, stateSpace = stateSpace))
#'
#' # pipe-friendly:
#' pop |> convertToLongFormat(stateSpace = stateSpace) |> head()
#'
#' @section Deprecated:
#' `convertToLongFormat()` is deprecated as of MicSim 4.0.0 and will be removed in a
#' future release. It still works, and shows a message once per R session
#' naming its replacement: use [as_long()] on the result of
#' [simulate.micsim_model()]. Set
#' `options(MicSim.show_deprecation = FALSE)` to turn the message off.
#'
#' @keywords internal
#' @export
convertToLongFormat <- function(pop,
                                migr           = FALSE,
                                stateSpace     = NULL,
                                absStates      = NULL,
                                allTransitions = NULL,
                                simHorizon     = NULL,
                                immigrPop      = NULL,
                                maxAge         = NULL) {
  deprecate_message("convertToLongFormat")
  convert_to_long(pop = pop, migr = migr, stateSpace = stateSpace,
                  absStates = absStates, allTransitions = allTransitions,
                  simHorizon = simHorizon, immigrPop = immigrPop,
                  maxAge = maxAge)
}


#' Convert to long format without the deprecation message
#'
#' @description
#' The body of [convertToLongFormat()]. `convertToLongFormat()` shows the
#' deprecation message and calls this; [as_long()] calls it directly.
#'
#' @param pop,migr,stateSpace,absStates,allTransitions,simHorizon,immigrPop,maxAge
#'   As in [convertToLongFormat()].
#' @return The long-format data frame, as [convertToLongFormat()] documents
#'   it.
#' @keywords internal
#' @noRd
convert_to_long <- function(pop, migr = FALSE, stateSpace = NULL,
                            absStates = NULL, allTransitions = NULL,
                            simHorizon = NULL, immigrPop = NULL,
                            maxAge = NULL) {

  spec <- micsim_spec(pop)

  absStates  <- from_spec(absStates,  spec$absStates,  "absStates")
  simHorizon <- from_spec(simHorizon, spec$simHorizon, "simHorizon")
  if (migr) {
    immigrPop <- from_spec(immigrPop, spec$immigrPop, "immigrPop")
  }
  if (is.null(maxAge)) maxAge <- spec$maxAge
  state_names <- state_variable_names(stateSpace, spec, pop)
  n_vars <- length(state_names)

  if (!is.null(allTransitions) && is.vector(allTransitions)) {
    allTransitions <- matrix(allTransitions, ncol = 2L, nrow = 1L)
  }

  # `motherID` is about an individual, not about an episode, so it is dropped
  if ("motherID" %in% colnames(pop)) {
    pop <- pop[, colnames(pop) != "motherID", drop = FALSE]
  }

  # ------------------------------------------------------------------
  # One row per transition becomes one row per episode.
  #
  # The trick that makes this work: a transition row already describes the
  # END of a spell. Its `From` is the state the individual was in, and its
  # `transitionTime` is when they left it. So the transition rows ARE the
  # episodes, once each one is given a start date and the individuals who
  # never left their last state are given a final censored row.
  # ------------------------------------------------------------------
  long <- data.frame(
    ID        = pop$ID,
    birthDate = pop$birthDate,
    birthyear = ifelse(is.na(pop$birthDate), NA,
                       1970 + getInDays(pop$birthDate) / DAYS_PER_YEAR),
    Tstop     = pop$transitionTime,
    .to       = as.character(pop$To),
    stringsAsFactors = FALSE
  )

  # Use the state occupied during each episode. Individuals without a
  # transition have no `From` value, so use their initial state instead.
  # Replace these missing labels before splitting them to give every
  # row one value per state variable.
  episode_state <- ifelse(is.na(pop$From),
                          as.character(pop$initState),
                          as.character(pop$From))
  from_state <- do.call(rbind, strsplit(episode_state, "/", fixed = TRUE))
  never_moved <- which(is.na(pop$From))
  for (i in seq_len(n_vars)) {
    long[[state_names[i]]] <- from_state[, i]
  }

  # Initially mark entries and exits as observed. A row without a transition
  # has no observed exit, so that is adjusted here. Initial-entry flags are
  # adjusted later in the code once the episodes are grouped by individual.
  long$statusEntry <- rep(1L, nrow(long))
  long$statusExit  <- rep(1L, nrow(long))
  long$statusExit[never_moved] <- 0L

  # Fill missing end dates with the simulation end or age limit,
  # whichever comes first.
  open_spell <- is.na(long$Tstop)
  long$Tstop[open_spell] <- end_of_observation(long$birthDate[open_spell],
                                               simHorizon, maxAge)

  long$OD <- vapply(
    seq_len(nrow(pop)),
    function(i) name_transition(pop$From[i], pop$To[i], absStates,
                                allTransitions),
    character(1L)
  )

  long <- long[order(as.numeric(long$ID)), , drop = FALSE]

  # ------------------------------------------------------------------
  # The final censored episode.
  #
  # Someone who ended the simulation alive, and whose last row records a
  # transition rather than a censoring, is still living the spell that
  # transition started. That spell is not in `pop`, because nothing
  # happened to end it, so it is added here.
  # ------------------------------------------------------------------
  long <- append_final_censoring(long, absStates, simHorizon, state_names,
                                 maxAge)

  # ------------------------------------------------------------------
  # Episode start dates.
  #
  # Each episode begins when the preceding episode ended. For the first
  # episode, use the simulation start or birth date, then adjust for
  # immigrant arrival when requested.
  # ------------------------------------------------------------------
  long <- long[order(as.numeric(long$ID)), , drop = FALSE]
  long$Tstart <- c(NA_character_, long$Tstop[-nrow(long)])
  first_of_individual <- !duplicated(long$ID)
  long$Tstart[first_of_individual] <- NA_character_

  start_year <- 1970 + getInDays(simHorizon[1L]) / DAYS_PER_YEAR
  born_during <- which(is.na(long$Tstart) & long$birthyear >= start_year)
  born_before <- which(is.na(long$Tstart) & long$birthyear <  start_year)

  long$Tstart[born_during] <- long$birthDate[born_during]
  long$Tstart[born_before] <- as.character(simHorizon[1L])
  long$statusEntry[born_before] <- 0L

  if (migr) {
    long <- adjust_immigrant_entry(long, immigrPop)
  }

  # ------------------------------------------------------------------
  # Number the episodes.
  # ------------------------------------------------------------------
  long <- long[order(long$Tstart), , drop = FALSE]
  long <- long[order(as.numeric(long$ID)), , drop = FALSE]

  counts <- table(long$ID)
  long$ns <- as.integer(counts[as.character(long$ID)])
  long$Episode <- sequence(as.integer(counts[unique(as.character(long$ID))]))

  long <- long[, c("ID", "birthDate", "Tstart", "Tstop", "statusEntry",
                   "statusExit", "OD", "ns", "Episode", state_names),
               drop = FALSE]
  row.names(long) <- NULL
  long
}


#' Name the transition that ended an episode
#'
#' @description
#' Works out what to put in the `OD` column: the user's own transition
#' pattern where one matches, the absorbing state's name where the spell
#' ended in one, and `"cens"` where the spell did not end at all.
#'
#' @details
#' For a transition between transient states, split the two labels into
#' state-variable values and compare each position. `"f/DE/0"` and
#' `"f/DE/1"` give `"f->f"`, `"DE->DE"` and `"0->1"`. If exactly one
#' of those patterns was supplied in `allTransitions`, use it. Otherwise
#' build a label from the values that changed.
#'
#' Note what follows from the component-by-component comparison: only a
#' pattern naming ONE state variable can ever match. A model that declares
#' its first birth as `"f/DE/0->f/DE/1"` gets `"0->1"` in `OD`, because
#' the full label is not among the candidates. This behaviour is fixed and
#' relied upon by downstream analysis.
#'
#' @param from,to State labels, either of which may be `NA`.
#' @param absStates The absorbing states.
#' @param allTransitions The user's transient transitions, or `NULL`.
#' @return A single string.
#' @keywords internal
#' @noRd
name_transition <- function(from, to, absStates, allTransitions) {

  if (is.na(from) || is.na(to)) return("cens")

  from <- as.character(from)
  to   <- as.character(to)
  if (to %in% absStates) {
    return(to)
  }

  from_parts <- strsplit(from, "/", fixed = TRUE)[[1L]]
  to_parts   <- strsplit(to,   "/", fixed = TRUE)[[1L]]

  if (!is.null(allTransitions)) {
    candidates <- paste(from_parts, to_parts, sep = "->")
    declared <- which(allTransitions[, 1L] %in% candidates)
    if (length(declared) == 1L) {
      return(allTransitions[declared, 1L])
    }
  }

  differing <- which(from_parts != to_parts)
  if (length(differing) == 0L) return("cens")
  paste(paste(from_parts[differing], to_parts[differing], sep = "->"),
        collapse = "/")
}


#' Add the final, censored episode
#'
#' @description
#' Adds a final episode for each individual whose last recorded
#' transition led to a transient state. The new episode describes their
#' time in that state until observation ended.
#'
#' @details
#' The last transition row describes the state the individual left.
#' Its internal `.to` column holds the destination state, which is split
#' into the new episode's state-variable columns. 
#' 
#' Individuals who reached an absorbing state are excluded: their last
#' spell ended, in whichever absorbing state it was. So are individuals
#' who already have a censored row, which is everyone who never
#' transitioned at all.
#'
#' @param long The episode table so far, ordered by ID.
#' @param absStates The absorbing states.
#' @param simHorizon The simulation horizon.
#' @param state_names Names of the state-variable columns.
#' @param maxAge The maximum age in years, or `NULL`.
#' @return The episode table with the extra rows appended.
#' @keywords internal
#' @noRd
append_final_censoring <- function(long, absStates, simHorizon, state_names,
                                   maxAge = NULL) {

  ended <- long$ID[long$OD %in% absStates | long$statusExit == 0L]
  still_running <- setdiff(unique(long$ID), ended)
  if (length(still_running) == 0L) return(long)

  theirs <- long[long$ID %in% still_running, , drop = FALSE]
  last_row <- c(diff(as.numeric(theirs$ID)), 1) != 0     # last row per ID
  extra <- theirs[last_row, , drop = FALSE]

  for (i in seq_len(nrow(extra))) {
    values <- strsplit(extra$.to[i], STATE_SEP, fixed = TRUE)[[1L]]
    if (length(values) == length(state_names)) {
      extra[i, state_names] <- as.list(values)
    }
  }

  extra$statusExit <- 0L
  extra$OD <- "cens"
  extra$Tstop <- end_of_observation(extra$birthDate, simHorizon, maxAge)

  rbind(long, extra)
}


#' Find the end date for episodes without an observed exit
#'
#' @description
#' The end of the simulation, or the day the individual reaches `maxAge`
#' minus the engine's margin, whichever is earlier. Past that day the
#' simulation puts nobody at risk, so no exposure may be counted there.
#'
#' @param birthDate A vector of birth dates in `yyyymmdd` format.
#' @param simHorizon The simulation horizon.
#' @param maxAge The maximum age in years, or `NULL` for "the end of the
#'   simulation".
#' @return A character vector of end dates in `yyyymmdd` format, with one
#'   date per supplied birth date.
#' @keywords internal
#' @noRd
end_of_observation <- function(birthDate, simHorizon, maxAge) {
  end_days <- getInDays(simHorizon[2L])
  if (length(birthDate) == 0L) return(character(0))
  if (is.null(maxAge)) return(rep(as.character(simHorizon[2L]), length(birthDate)))
  capped <- pmin(end_days,
                 getInDays(birthDate) + (maxAge - MAX_AGE_MARGIN) * DAYS_PER_YEAR)
  out <- as.character(getInDateFormat(capped))
  # The horizon's own spelling where it is the binding limit, so that a
  # result without maxAge censoring is unchanged character for character.
  out[capped >= end_days] <- as.character(simHorizon[2L])
  out
}


#' Start an immigrant's first episode at their arrival
#'
#' @description
#' Sets each immigrant's first episode start date to their arrival date
#' and marks that entry as observed (`statusEntry = 1`). Later episodes
#' are unchanged.
#'
#' @details
#' This convention treats arrival as entry into the study population.
#' It does not imply that the individual's current state began on that
#' date; time spent in that state before arrival is outside observation.
#'
#' @param long The episode table, with `Tstart` filled in.
#' @param immigrPop The immigrant population.
#' @return The episode table with immigrants' first episodes adjusted.
#' @keywords internal
#' @noRd
adjust_immigrant_entry <- function(long, immigrPop) {
  is_first <- !duplicated(long$ID) & long$ID %in% immigrPop$ID
  arrival <- immigrPop$immigrDate[match(long$ID[is_first], immigrPop$ID)]
  long$Tstart[is_first] <- as.character(arrival)
  long$statusEntry[is_first] <- 1L
  long
}


#' Read the model context MicSim attached to its output
#'
#' @description
#' `micSim()` records the few things the converters need on its return
#' value, so that the common case needs no extra arguments.
#'
#' @details
#' Some transformations and export/import steps remove custom attributes.
#' If `micsim_spec` is lost, the required conversion settings must be
#' supplied explicitly. This helper reads only the attribute; it does not
#' look for settings in the user's global environment.
#'
#' @param pop The population data frame returned by [micSim()].
#' @return A named list, possibly empty.
#' @keywords internal
#' @noRd
micsim_spec <- function(pop) {
  spec <- attr(pop, "micsim_spec")
  if (is.null(spec)) list() else spec
}


#' Take an argument, or fall back to what MicSim recorded
#'
#' @param given The value the user passed, or `NULL`.
#' @param recorded The value from the spec attribute, or `NULL`.
#' @param name The argument's name, for the error message.
#' @return `given` if it is not `NULL`, otherwise `recorded`.
#' @keywords internal
#' @noRd
from_spec <- function(given, recorded, name) {
  if (!is.null(given)) return(given)
  if (!is.null(recorded)) return(recorded)
  micsim_abort(
    "micsim_convert_error",
    paste0(
      "`", name, "` was not supplied and `pop` does not record it.\n",
      "  Pass it explicitly."
    )
  )
}


#' Choose names for the state-variable columns
#'
#' @description
#' Uses names from `stateSpace` when supplied, otherwise names recorded
#' in `spec$stateVars`. If neither source is available, counts the values
#' in the first initial-state label and generates `var1`, `var2`, and
#' so on.
#'
#' @details
#' A transition matrix carries the state SPACE in its dimension names but
#' not the NAMES of the state variables, so a result converted without
#' `stateSpace` gets `var1`, `var2` and so on. Passing the same
#' `stateSpace` that built the transition matrix restores `sex`, `nat`,
#' `fert`.
#'
#' @param stateSpace The user's state space, or `NULL`.
#' @param spec A list of model settings recorded on the simulation output.
#' @param pop The population data frame returned by [micSim()].
#' @return A character vector of state-variable column names.
#' @keywords internal
#' @noRd
state_variable_names <- function(stateSpace, spec, pop) {

  if (!is.null(stateSpace)) {
    if (is.null(dim(stateSpace))) {
      # A one-variable state space may be given as a bare vector carrying
      # its name as an attribute.
      nm <- attr(stateSpace, "name")
      return(if (is.null(nm)) "var1" else nm)
    }
    return(colnames(stateSpace))
  }

  if (!is.null(spec$stateVars)) return(spec$stateVars)

  width <- length(strsplit(as.character(pop$initState[1L]), "/",
                           fixed = TRUE)[[1L]])
  paste0("var", seq_len(width))
}
