# =====================================================================
# convert_wide.R
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   Reshapes simulation output into one row per individual, with numbered
#   columns describing each recorded transition.
#
# PLACE IN THE PIPELINE
#   engine_store.R --(the `pop` data frame)--> convert_wide.R --> the user
#
# WHAT YOU NEED TO KNOW
#   The simulation output records each transition on a separate row.
#   Wide format places these transitions side by side: From.1, To.1,
#   transitionTime.1, transitionAge.1, From.2, and so on.
#
#   This layout is useful for questions about individuals, such as their
#   age at first birth. For analyses of time spent in different states,
#   use convertToLongFormat().
#
#   Transition numbering follows the existing row order within each
#   individual. The recorded state labels, dates and ages are retained.
#
# GUARANTEES
#   * Exactly one row per individual in `pop`.
#   * `ns` gives the number of recorded transitions, including 0 for
#     individuals who never transitioned.
# =====================================================================


#' Convert simulation output to wide format
#'
#' @description
#' Reshapes the population table returned by [micSim()] into one row per
#' individual, with numbered columns for their recorded transitions.
#'
#' @details
#' Each transition is represented by four columns: its origin state,
#' destination state, date and age. The first transition appears in
#' `From.1`, `To.1`, `transitionTime.1` and `transitionAge.1`; the second
#' appears in the corresponding `.2` columns, and so on.
#'
#' Transitions are numbered in their existing order within each
#' individual, so the input rows should already be in chronological
#' order. Individuals with fewer transitions have `NA` in the remaining
#' numbered columns.
#'
#' `ns` counts an individual's transitions. Someone who never transitioned
#' occupies one row of `pop`, with `NA` in all four transition columns,
#' and is reported here with `ns = 0`.
#'
#' The output omits `motherID`. If parent links are needed, retain a
#' separate table of `ID` and `motherID` from `pop`, then join it to the
#' wide result by `ID`.
#'
#' The function takes a single data frame, so it is pipe-friendly:
#' `pop |> convertToWideFormat()`.
#'
#' @param pop The population data frame returned by [micSim()].
#'
#' @return A data frame with one row per individual and columns `ID`,
#'   `birthDate`, `initState`, `ns`, and numbered groups of `From`, `To`,
#'   `transitionTime` and `transitionAge`.
#'
#' @seealso [micSim()], [convertToLongFormat()], [MicSim-concepts]
#'
#' @examples
#' stateSpace <- c("m", "f")
#' mortRates <- function(age, calTime) 0.00003 * exp(0.1 * age)
#' tm <- buildTransitionMatrix(NULL, cbind("dead", "mortRates"), stateSpace)
#'
#' set.seed(1)
#' pop <- micSim(
#'   initPop = data.frame(ID = 1:10, birthDate = rep("19500101", 10),
#'                        initState = rep("f", 10),
#'                        stringsAsFactors = FALSE),
#'   transitionMatrix = tm, maxAge = 120,
#'   simHorizon = c(20000101, 20401231), verbose = FALSE
#' )
#' convertToWideFormat(pop)
#'
#' # convertToWideFormat() is pipe-friendly:
#' pop |> convertToWideFormat()
#'
#' @section Deprecated:
#' `convertToWideFormat()` is deprecated as of MicSim 4.0.0 and will be removed in a
#' future release. It still works, and shows a message once per R session
#' naming its replacement: use [as_wide()] on the result of
#' [simulate.micsim_model()]. Set
#' `options(MicSim.show_deprecation = FALSE)` to turn the message off.
#'
#' @keywords internal
#' @export
convertToWideFormat <- function(pop) {
  deprecate_message("convertToWideFormat")
  convert_to_wide(pop)
}


#' Convert to wide format without the deprecation message
#'
#' @description
#' The body of [convertToWideFormat()]. `convertToWideFormat()` shows the
#' deprecation message and calls this; [as_wide()] calls it directly.
#'
#' @param pop As in [convertToWideFormat()].
#' @return The wide-format data frame, as [convertToWideFormat()] documents
#'   it.
#' @keywords internal
#' @noRd
convert_to_wide <- function(pop) {

  if ("motherID" %in% colnames(pop)) {
    pop <- pop[, colnames(pop) != "motherID", drop = FALSE]
  }

  pop <- pop[order(as.numeric(pop$ID)), , drop = FALSE]

  # Transitions per individual, and their numbering within the individual.
  # `table()` counts rows, and a row is a transition except for the single
  # NA row of someone who never transitioned; that case is corrected at
  # the end, once the reshape has used the row.
  counts <- table(pop$ID)
  pop$ns <- as.integer(counts[as.character(pop$ID)])
  pop$Episode <- sequence(as.integer(counts[unique(as.character(pop$ID))]))

  pop <- pop[, c("ID", "birthDate", "initState", "ns", "Episode",
                 "From", "To", "transitionTime", "transitionAge"),
             drop = FALSE]

  wide <- stats::reshape(
    pop,
    timevar   = "Episode",
    idvar     = "ID",
    direction = "wide",
    v.names   = c("From", "To", "transitionTime", "transitionAge")
  )

  # Someone who never transitioned still occupied a row, so `table()`
  # counted one. Their life course contains no transitions.
  wide$ns[is.na(wide$transitionTime.1)] <- 0L

  row.names(wide) <- NULL
  wide
}
