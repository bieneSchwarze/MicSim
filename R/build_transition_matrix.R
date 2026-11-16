# =====================================================================
# build_transition_matrix.R
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   Turns a short, readable description of which transitions a model allows
#   into the full matrix of every state pair, with the name of the rate
#   function in each cell.
#
# PLACE IN THE PIPELINE
#   user-defined stateSpace, transition patterns and rate-function names
#       --> buildTransitionMatrix() [this file] --> transitionMatrix
#       --> micSim() [compat.R] --> rate lookup [registry.R] --> engine.R
#
#   The user defines stateSpace, usually with expand.grid(), and writes
#   the rate functions. This helper builds the matrix; micSim() later
#   looks up the named functions. The model-object interface instead uses
#   transition() declarations and does not need this matrix helper.
#
# WHAT YOU NEED TO KNOW  (concepts in ?MicSim-concepts)
#   composite state, absorbing state.
#
#   A state combines one value from each of the user's state variables.
#   For example, a model might include:
#
#       sex:        "m", "f"
#       parity:     "0", "1", "2"
#       education:  "no", "low", "medium", "high"
#                   (none, elementary, secondary, higher)
#
#   These combinations give 24 transient states. With one absorbing
#   state ("dead"), the matrix has 24 rows, 25 columns and 600 cells.
#
#   Rather than fill every cell by hand, the user supplies patterns such as:
#
#       "no->low"      education changes from none to elementary;
#                      sex and parity stay the same
#       "f/0->f/1"     parity changes from zero to one for women;
#                      education stays the same
#       "dead"         applies a mortality rate to every transient state
#       "f/low/dead"   restricts that mortality specification to
#                      states containing both "f" and "low"
#
#   Each pattern is supplied with the name of its rate function, and
#   this file fills the corresponding matrix cells.
#
#   The result has one row per transient state and one column per transient
#   or absorbing state. A cell contains "0" when a transition is not allowed,
#   and otherwise the name of the function giving its rate.
#
# GUARANTEES
#   * Row names are the transient states; column names are the transient
#     states followed by the absorbing ones. That column ORDER is a fixed
#     contract: it is the order of the state space, which decides which of
#     an individual's random numbers each destination receives, so changing
#     it may change the result of a run for a given seed.
# =====================================================================


#' Build a transition matrix from transition patterns
#'
#' @description
#' Expands a short description of which transitions a model allows into
#' the full matrix that [micSim()] expects.
#'
#' @details
#' **Transient transitions** are given as `"origin->destination"`.
#' Each side names the state-variable values needed to describe the change.
#' The helper finds origins containing the values on the left and destinations
#' containing the values on the right. Values not mentioned in the pattern
#' must be the same before and after the transition. So in a model of sex, 
#' education and parity, `"0->1"` expands to a first birth for every 
#' combination of sex and education, while `"f/0->f/1"` restricts it to women.
#'
#' **Absorbing transitions** are given as a state name rather than a pair,
#' because the origin is "anywhere". `"dead"` makes every transient state
#' able to reach death; a compound label such as `"f/low/dead"` restricts
#' it to the states whose other variables match.
#'
#' The second column of each matrix is the NAME of the function supplying
#' the rate, as a string. Storing names rather than the functions
#' themselves is the convention [micSim()] expects; it resolves each name
#' once, in the environment it was called from.
#'
#' **Any label can name an absorbing state.** `"dead"`, `"rest"`,
#' `"abroad"` or `"institution"` are all fine, and a model may have as
#' many as it needs. In a compound label the absorbing part is the
#' component that is not a value of any state variable, so `"f/IT/abroad"`
#' is the state `"abroad"`, reachable from Italian women.
#'
#' @section Assumptions:
#' The absorbing states are the columns after the transient ones, and
#' their ORDER is the order of the absorbing states in the state space,
#' which decides which of an individual's random numbers each destination
#' receives. It is `"dead"` first and `"rest"` second, where present,
#' then every other absorbing state in the order `absTransitions`
#' first names it. `"rest"` conventionally means emigration.
#'
#' A label that could be read two ways is refused: an absorbing state
#' that is also the value of a state variable, or a compound label with
#' no part, or more than one part, that is not a state value.
#'
#' A pattern that matches nothing is not an error. It produces a matrix
#' with no cell set, which is usually a typo in a state value; check the
#' result with `table(transitionMatrix != "0")` if a model behaves
#' unexpectedly.
#'
#' @param allTransitions Two-column matrix of transitions between
#'   transient states: the pattern `"origin->destination"` and the name of
#'   the rate function. A single transition may be given as a length-2
#'   vector. `NULL` when the model has none.
#' @param absTransitions Two-column matrix of transitions into absorbing
#'   states: the absorbing state's label and the name of the rate
#'   function. A single one may be given as a length-2 vector.
#' @param stateSpace User-defined combinations of state-variable values, as
#'   produced by [base::expand.grid()]: one row per transient state and one
#'   column per variable. A vector is accepted for a single variable.
#'
#' @return A character matrix with one row per transient state and one
#'   column per state, holding `"0"` where no transition is possible and
#'   otherwise the name of the rate function.
#'
#' @seealso [micSim()], [MicSim-concepts]
#'
#' @examples
#' stateSpace <- expand.grid(sex = c("m", "f"), fert = c("0", "1"))
#' buildTransitionMatrix(
#'   allTransitions = cbind("f/0->f/1", "fertRates"),
#'   absTransitions = cbind("dead", "mortRates"),
#'   stateSpace     = stateSpace
#' )
#'
#' @section Deprecated:
#' `buildTransitionMatrix()` is deprecated as of MicSim 4.0.0 and will be removed in a
#' future release. It still works, and shows a message once per R session
#' naming its replacement: declare the transitions with [transition()], [move()]
#' and [transitions()] and pass them to [micsim_model()]. Set
#' `options(MicSim.show_deprecation = FALSE)` to turn the message off.
#'
#' @keywords internal
#' @export
buildTransitionMatrix <- function(allTransitions, absTransitions,
                                  stateSpace) {
  deprecate_message("buildTransitionMatrix")

  # ------------------------------------------------------------------
  # Accept the shorthand forms the classic interface accepts: a single
  # transition as a bare length-2 vector, a one-variable state space as a
  # bare character vector.
  # ------------------------------------------------------------------
  if (is.vector(allTransitions) && !is.null(allTransitions)) {
    allTransitions <- matrix(allTransitions, ncol = 2L, nrow = 1L)
  }
  if (is.vector(absTransitions)) {
    absTransitions <- matrix(absTransitions, ncol = 2L, nrow = 1L)
  }
  if (is.null(dim(stateSpace))) {
    stateSpace <- matrix(stateSpace, ncol = 1L)
  }

  if (is.null(absTransitions) || ncol(absTransitions) != 2L) {
    micsim_abort(
      "micsim_transition_matrix_error",
      paste0(
        "`absTransitions` must have two columns: the absorbing state and\n",
        "  the name of the function giving its rate. Every model needs at\n",
        "  least one absorbing state, usually \"dead\"."
      )
    )
  }
  if (!is.null(allTransitions) && ncol(allTransitions) != 2L) {
    micsim_abort(
      "micsim_transition_matrix_error",
      paste0(
        "`allTransitions` must have two columns: the transition pattern\n",
        "  \"origin->destination\" and the name of its rate function."
      )
    )
  }

  absStates <- absTransitions[, 1L]
  state_labels <- apply(stateSpace, 1L, paste, collapse = "/")
  state_values <- unique(as.character(unlist(stateSpace)))

  # Patterns match on values, whichever variable holds them, so a value
  # used by two variables makes them ambiguous: with sex coded 1/2 and
  # parity 0/1/2, "1->2" is both a sex change and a second birth. This
  # is accepted with a warning; distinct values avoid this ambiguity.
  per_variable <- lapply(seq_len(ncol(stateSpace)), function(j)
    unique(as.character(stateSpace[, j])))
  shared <- unique(unlist(per_variable)[duplicated(unlist(per_variable))])
  if (length(shared) > 0L) {
    warning(micsim_warning(
      "micsim_ambiguous_values_warning",
      paste0(
        "State value(s) ", toString(shared), " are used by more than one\n",
        "  state variable, so transition patterns naming them are ambiguous\n",
        "  and may match transitions you did not mean. Give each variable\n",
        "  its own values (e.g. \"m\"/\"f\" rather than 1/2 for sex)."
      )
    ))
  }

  # Which absorbing states does this model use? Whatever absTransitions
  # names: in each label, the one component that is not a state value.
  # legacy_absorbing_order() puts "dead" and "rest" first when present.
  absorbing_part <- vapply(absStates, function(label) {
    parts <- unlist(strsplit(label, "/", fixed = TRUE))
    candidates <- parts[!parts %in% state_values]
    if (length(candidates) != 1L) {
      micsim_abort(
        "micsim_transition_matrix_error",
        paste0(
          "Cannot tell which part of `", label, "` is the absorbing state.\n",
          "  An absorbing label is one name that is not a value of any state\n",
          "  variable, optionally preceded by state values that restrict who\n",
          "  can reach it, e.g. \"dead\" or \"f/IT/dead\"."
        )
      )
    }
    candidates
  }, character(1L), USE.NAMES = FALSE)
  absorbing_names <- legacy_absorbing_order(unique(absorbing_part))

  # The empty matrix. Column order is transient states first, absorbing
  # states last; this order affects random-number allocation (see header).
  transitionMatrix <- matrix(
    "0",
    nrow = nrow(stateSpace),
    ncol = nrow(stateSpace) + length(absorbing_names),
    dimnames = list(state_labels, c(state_labels, absorbing_names))
  )

  # ------------------------------------------------------------------
  # Absorbing transitions.
  #
  # A bare label such as "dead" applies to every transient state. A
  # compound one such as "f/low/dead" applies only to the states whose
  # other variables match, so the label is split into the absorbing part
  # and the restriction.
  # ------------------------------------------------------------------
  for (i in seq_len(nrow(absTransitions))) {
    parts <- unlist(strsplit(absStates[i], "/", fixed = TRUE))
    column <- which(colnames(transitionMatrix) == absorbing_part[i])

    if (length(parts) == 1L) {
      transitionMatrix[, column] <- absTransitions[i, 2L]
    } else {
      restriction <- parts[parts != absorbing_part[i]]
      rows <- which(apply(stateSpace, 1L, state_contains, values = restriction))
      transitionMatrix[rows, column] <- absTransitions[i, 2L]
    }
  }

  # ------------------------------------------------------------------
  # Transitions between transient states.
  #
  # For each pattern, find every origin state containing the pattern's
  # origin values and every destination state containing its destination
  # values, then keep the pairs that agree on everything the pattern did
  # not mention. That last test is what makes "0->1" mean "a first birth,
  # with sex and nationality unchanged".
  # ------------------------------------------------------------------
  if (!is.null(allTransitions)) {
    ends <- do.call(rbind, strsplit(allTransitions[, 1L], "->", fixed = TRUE))

    for (i in seq_len(nrow(ends))) {
      origin_values      <- unlist(strsplit(ends[i, 1L], "/", fixed = TRUE))
      destination_values <- unlist(strsplit(ends[i, 2L], "/", fixed = TRUE))

      origins <- stateSpace[
        apply(stateSpace, 1L, state_contains, values = origin_values), ,
        drop = FALSE
      ]
      destinations <- stateSpace[
        apply(stateSpace, 1L, state_contains, values = destination_values), ,
        drop = FALSE
      ]

      for (j in seq_len(nrow(origins))) {
        origin <- as.character(unlist(origins[j, ]))
        for (k in seq_len(nrow(destinations))) {
          destination <- as.character(unlist(destinations[k, ]))

          # What the pattern did not mention, on each side.
          rest_origin      <- origin[!origin %in% origin_values]
          rest_destination <- destination[!destination %in% destination_values]

          if (setequal(rest_origin, rest_destination)) {
            row    <- which(rownames(transitionMatrix) ==
                              paste(origin, collapse = "/"))
            column <- which(colnames(transitionMatrix) ==
                              paste(destination, collapse = "/"))
            transitionMatrix[row, column] <- allTransitions[i, 2L]
          }
        }
      }
    }
  }

  transitionMatrix
}


#' Does a state contain all of these values?
#'
#' @description
#' The matching test behind every pattern in
#' [buildTransitionMatrix()]: a state matches a pattern when every value
#' the pattern names appears somewhere in the state.
#'
#' @details
#' Note what it does not do: it does not care which state variable a value
#' belongs to. `"f"` matches whichever variable happens to hold it. That
#' is what lets a user write `"0->1"` without naming the parity variable,
#' and it is also why distinct values across variables are recommended.
#' A value used by two variables can match either; the matrix builder warns
#' about this ambiguity but does not reject it.
#'
#'
#' @param state One row of the state space.
#' @param values Character vector of values the pattern names.
#' @return `TRUE` or `FALSE`.
#' @keywords internal
#' @noRd
state_contains <- function(state, values) {
  all(values %in% as.character(unlist(state)))
}
