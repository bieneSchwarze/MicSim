# =====================================================================
# engine_store.R
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   Records every transition as it happens, cheaply, and turns the record
#   into a data frame for users.
#
# PLACE IN THE PIPELINE
#   engine.R --(one row per event)--> engine_store.R --(the `pop` frame)--> result
#                                            ^
#                                            |
#                       engine_person.R --(who these slots are)
#
# WHAT YOU NEED TO KNOW  (concepts in ?MicSim-concepts)
#   episode, transition, absorbing state.
#
#   THE OUTPUT FORMAT. `pop` has the columns ID, birthDate, initState, From, 
#   To, transitionTime, transitionAge, and motherID when fertility is modelled.
#   It has one row per transition, plus one row for every individual who never 
#   transitioned at all.  Everything downstream -- as_long and as_wide, 
#   every estimator, every user script -- reads that shape, so
#   it must not change.
#
# GUARANTEES
#   * Appending a row is amortised O(1). A run of E events costs O(E).
#   * Nothing is decoded, formatted or merged until the run is over.
#     During the run the store holds integers and day counts only.
#   * assemble_pop() reproduces the output format as documented above, including the
#     NA row for individuals who never moved.
#
# UNITS
#   `time_days` and `age_days` are day counts on the getInDays() scale,
#   possibly fractional. Only assemble_pop() turns them into a yyyymmdd string 
#   and an age in years rounded to two decimals.
# =====================================================================


#' Create an empty transition store
#'
#' @description
#' Collects one row per transition during a run. Each row is stored as
#' five numbers per event; turning them into labels, dates and a data frame is 
#' left to `assemble_pop()` at the end of a run.
#'
#' @details
#' **Why so little per row.** The engine writes to this store on every
#' event, so anything done here is done tens of millions of times in a
#' large run. Turning a state index into `"f/DE/1"` and a day count into
#' `"20170704"` is work that only has to happen once per row *at the end*. 
#' During the run the store holds an integer slot, two integer state indices 
#' and two day counts.
#'
#' **Why slots rather than IDs.** A slot is a row number in the
#' person store, so the user-facing `ID` is one vector lookup away and
#' never has to be carried around. The event queue uses slots for the same
#' reason.
#'
#' **Growth by doubling.** When the arrays fill they double, so appending
#' E rows costs O(E) copying in total rather than the O(E^2) of growing a
#' matrix one `rbind()` at a time.
#'
#' @param capacity Initial number of rows. The store doubles when full, so
#'   this only affects performance.
#'
#' @return An object of class `transition_store`: a list of functions.
#'   \describe{
#'     \item{`append(slot, from_idx, to_idx, time_days, age_days)`}{
#'       record one transition}
#'     \item{`as_columns()`}{the trimmed columns, as a named list}
#'     \item{`size()`}{number of recorded transitions}
#'   }
#'
#' @seealso `assemble_pop()` turns a store into the contract data frame;
#'   [MicSim-concepts].
#'
#' @examples
#' store <- new_transition_store()
#' store$append(slot = 1L, from_idx = 3L, to_idx = 4L,
#'              time_days = 17351, age_days = 13699)
#' store$size()
#'
#' @keywords internal
#' @noRd
new_transition_store <- function(capacity = 1024L) {

  e <- new.env(parent = emptyenv())
  e$n   <- 0L
  e$cap <- as.integer(capacity)

  e$slot      <- integer(e$cap)   # row number in the person store
  e$from      <- integer(e$cap)   # state index left
  e$to        <- integer(e$cap)   # state index entered
  e$time_days <- numeric(e$cap)   # when, on the getInDays() scale
  e$age_days  <- numeric(e$cap)   # the individual's age then, in days

  grow <- function() {
    extra_rows <- e$cap
    e$slot      <- c(e$slot,      integer(extra_rows))
    e$from      <- c(e$from,      integer(extra_rows))
    e$to        <- c(e$to,        integer(extra_rows))
    e$time_days <- c(e$time_days, numeric(extra_rows))
    e$age_days  <- c(e$age_days,  numeric(extra_rows))
    e$cap <- e$cap + extra_rows
  }

  append_row <- function(slot, from_idx, to_idx, time_days, age_days) {
    if (e$n == e$cap) grow()
    i <- e$n + 1L

    # Coerce first, so nothing below can fail half-way through.
    sl_v <- as.integer(slot)
    fr_v <- as.integer(from_idx)
    to_v <- as.integer(to_idx)
    tm_v <- as.numeric(time_days)
    ag_v <- as.numeric(age_days)

    # WRITE IN PLACE via the "steal-back" idiom, exactly as in
    # engine_person.R. `e$col[i] <- v` looks like an O(1) write, but the
    # column is reachable from the environment `e`, so R duplicates the
    # whole column on every call, and recording E events costs O(E^2)
    # (measured: 1.4 s for 20,000 events, 21 s for 80,000). Dropping the
    # store's reference first leaves the local as the sole owner, so the
    # write is in place; then the column is handed back.
    v <- e$slot;      e$slot      <- NULL; v[i] <- sl_v; e$slot      <- v
    v <- e$from;      e$from      <- NULL; v[i] <- fr_v; e$from      <- v
    v <- e$to;        e$to        <- NULL; v[i] <- to_v; e$to        <- v
    v <- e$time_days; e$time_days <- NULL; v[i] <- tm_v; e$time_days <- v
    v <- e$age_days;  e$age_days  <- NULL; v[i] <- ag_v; e$age_days  <- v

    e$n <- i
    invisible(NULL)
  }

  as_columns <- function() {
    keep <- seq_len(e$n)
    list(
      slot      = e$slot[keep],
      from      = e$from[keep],
      to        = e$to[keep],
      time_days = e$time_days[keep],
      age_days  = e$age_days[keep]
    )
  }

  structure(
    list(
      append     = append_row,
      as_columns = as_columns,
      size       = function() e$n,
      state      = e             # exposed for tests and diagnostics only
    ),
    class = "transition_store"
  )
}


#' Assemble the `pop` data frame
#'
#' @description
#' Turns the recorded transitions and the population register into the
#' data frame that `micSim()` returns, in exactly the format the file header
#' describes.
#' 
#' @details
#' The complete format, because it is easy to get subtly wrong and
#' everything downstream depends on it:
#'
#' \describe{
#'   \item{one row per transition}{`From` and `To` name the states left
#'     and entered, `transitionTime` the date, `transitionAge` the
#'     individual's age in years at that date, rounded to two decimals.}
#'   \item{one row per individual who never transitioned}{with `NA` in all
#'     four of those columns. This is how the long-format converter
#'     recognises an episode that was never completed, so these rows carry
#'     information and must not be dropped.}
#'   \item{no entry rows}{Joining the population, whether by being a member of
#'   the starting population, an immigrant or a newborn, does not create a row
#'   of its own.}
#'   \item{columns and their order}{`ID`, `birthDate`, `initState`,
#'     `From`, `To`, `transitionTime`, `transitionAge`, then `motherID`
#'     when fertility is modelled.}
#'   \item{row order}{by `ID`, then chronologically. `transitionAge` is
#'     rounded to two decimals, so two events of one individual less than
#'     ~3.65 days apart can share it; ordering by the rounded age alone
#'     would then leave them in arbitrary order. Ties are
#'     therefore broken by the order in which the engine processed the
#'     events, which is chronological.}
#' }
#'
#' `motherID` is taken straight from the person store, which already
#' records each individual's mother (`NA` for anyone not born during the
#' run), rather than merged in from a separate table of births.
#'
#' @section Assumptions:
#' `birthDate` and `transitionTime` are returned as **character**
#' `yyyymmdd`. Normalising `birthDate` to character is a deliberate
#' tidy-up, so its type does not depend on whether any birth occurred.
#'
#' @param store A `transition_store`.
#' @param persons A `person_store`. Supplies every individual who took
#'   part, including immigrants and those born during the run.
#' @param codec A `state_codec`, used to turn state indices back into
#'   labels.
#' @param include_mother Add the `motherID` column. The engine passes
#'   `TRUE` when the model has fertility transitions.
#'
#' @return A data frame with the columns described above.
#'
#' @seealso `new_transition_store()`, `new_person_store()`,
#'   [MicSim-concepts]
#'
#' @keywords internal
#' @noRd
assemble_pop <- function(store, persons, codec, include_mother = FALSE) {

  tr <- store$as_columns()
  people <- persons$as_columns()

  # Who took part: the starting population plus every immigrant and
  # newborn added to it during the run.
  init_pop <- data.frame(
    ID        = people$ID,
    birthDate = people$birthDate,
    initState = codec$decode(people$initState),
    stringsAsFactors = FALSE
  )

  # What happened. Decoding and formatting happen here, once per row, and
  # nowhere else. With no transitions recorded this is a zero-row frame,
  # and the merge below still produces one NA row per individual, which is
  # the correct result for a run in which nobody moved.
  transitions_out <- data.frame(
    ID             = people$ID[tr$slot],
    From           = codec$decode(tr$from),
    To             = codec$decode(tr$to),
    transitionTime = as.character(getInDateFormat(tr$time_days)),
    transitionAge  = round(tr$age_days / DAYS_PER_YEAR, 2),
    stringsAsFactors = FALSE
  )
  # Processing order of each event (the store appends in the order the
  # engine processed events, i.e. chronologically). Used only to break ties
  # in the final ordering, then dropped. merge() does NOT preserve the
  # order of rows within a key unless its input is sorted by that key,
  # which the store is not, so the order cannot be left to merge().
  transitions_out$.seq <- seq_len(nrow(transitions_out))

  pop <- merge(init_pop, transitions_out, all = TRUE, by = "ID")

  if (isTRUE(include_mother)) {
    pop$motherID <- people$motherID[match(pop$ID, people$ID)]
  }

  pop <- pop[order(as.numeric(pop$ID), as.numeric(pop$transitionAge),
                   pop$.seq), , drop = FALSE]
  pop$.seq <- NULL

  # Row names after a merge-then-reorder are a meaningless permutation of
  # integers, and users see them on every head(pop). Resetting them is
  # cosmetic.
  row.names(pop) <- NULL
  pop
}


#' Print a transition store
#'
#' Reports how many transitions have been recorded and the window they
#' span. Deliberately does not print any rows: the store is an internal
#' accumulator, and what a user should look at is the assembled `pop`.
#'
#' @param x A `transition_store`.
#' @param ... Not used; present for compatibility with the [print()] generic.
#' @return `x`, invisibly. Called for the side effect of printing.
#' @exportS3Method print transition_store
#' @noRd
print.transition_store <- function(x, ...) {
  n <- x$size()
  cat("<transition_store>\n")
  cat("  transitions recorded: ", n, "\n", sep = "")
  if (n > 0L) {
    cols <- x$as_columns()
    cat("  first event         : ", getInDateFormat(min(cols$time_days)),
        "\n", sep = "")
    cat("  last event          : ", getInDateFormat(max(cols$time_days)),
        "\n", sep = "")
    cat("  individuals involved: ", length(unique(cols$slot)), "\n", sep = "")
  }
  invisible(x)
}
