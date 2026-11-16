# =====================================================================
# engine_queue.R
# ---------------------------------------------------------------------
# WHAT THIS FILE DOES
#   Holds the events that have been scheduled but have not happened yet,
#   and hands back the earliest one on demand.
#
# PLACE IN THE PIPELINE
#   engine.R --(schedule an event)--> engine_queue.R --(earliest event)--> engine.R
#
#   This file contains NO demographic knowledge. It does not know what a
#   state is, what a transition means, or what a person is. All it knows is
#   "something is due to happen to slot 17 on day 20315". That separation
#   is deliberate, as it makes the queue testable on its own against an
#   obviously correct reference (sort()).
#
# WHAT YOU NEED TO KNOW  (concepts in ?MicSim-concepts)
#   discrete-event simulation, event queue, global clock.
#
#   A BINARY MIN-HEAP is an array that behaves like an almost-sorted list.
#   The rule is: the entry at position i is never later than the entries at
#   positions 2i and 2i+1. Draw the array as a binary tree by putting entry
#   1 at the top, entries 2 and 3 below it, entries 4..7 below those, and so
#   on; the rule says every parent comes before both of its children. That
#   single rule guarantees the earliest entry sits at position 1, so reading
#   it costs nothing, while inserting or removing costs only log2(Q) swaps,
#   because a value only ever has to travel up or down one path of the tree.
#
#   We never need the queue fully sorted. We only ever need "what is next",
#   which is exactly what a heap is for.
#
# GUARANTEES
#   * HEAP PROPERTY: for every occupied position i > 1, the entry at i does
#     not precede the entry at its parent i %/% 2.
#   * AT MOST ONE LIVE EVENT PER PERSON. Pushing a second event for the same
#     person silently replaces the first. An individual never has one pending 
#     event at any time.
#   * DETERMINISM: two events falling on the same day are returned in the
#     order they were pushed. Without this a fixed seed would not give a
#     fixed result.
#   * No `<<-`. All mutation writes into an environment held by the
#     closures, and every write is an ordinary assignment.
#
# UNITS
#   `time_days` is a day count on the same scale as getInDays(): days since
#   1970-01-01, fractional, because an event can fall part-way through a
#   day. `slot` and `to_idx` are integer indices, not labels.
# =====================================================================


#' Create an empty event queue
#'
#' @description
#' Builds a priority queue of scheduled events, ordered by the date on
#' which they are due. The earliest event can always be read or removed in
#' logarithmic time, which is what lets the engine's main loop jump from
#' one event straight to the next.
#'
#' @details
#' **What is stored.** One entry per scheduled event, holding three
#' numbers: the date it is due (`time_days`), which individual it belongs
#' to (`slot`, a row index into the person store), and which state that
#' individual will move into (`to_idx`, a state index from the codec).
#' Anything else the engine needs about the individual,
#' such as their birth date or current state, lives in the person store
#' and is looked up by `slot`.
#'
#' **Cancelling and rescheduling.**
#' Replacing or withdrawing a pending event means removing an arbitrary entry 
#' from the middle of a heap. That is possible but fiddly, because it needs 
#' a second index from slot to heap position that has to stay correct through 
#' every swap.
#'
#' This implementation uses the standard cheaper alternative, **lazy
#' deletion**. Every entry is stamped with a version number when it is
#' pushed, and the queue separately records, per individual, the version
#' of that individual's one live event. Cancelling just forgets that
#' version number, which is O(1) and touches nothing in the heap. The
#' stale entry stays in the array but is recognised on the way out,
#' discarded, and the next entry is examined instead. Rescheduling is the
#' same mechanism: pushing a new event for an individual who already has
#' one supersedes it.
#'
#' This costs memory, as cancelled entries occupy space until they surface.
#' The bound is the total number of pushes, which is the number of events
#' in the run, so it grows linearly with the run and is acceptable. 
#' `n_nodes()` reports the array occupancy and `size()` the number of live 
#' events, so a run that accumulates an unexpected amount of dead weight is visible.
#'
#' **Why the version stamp also breaks ties.** Version numbers come from a 
#' single counter that only increases, so they *are* the
#' insertion order. Comparing them when two events fall on the same day
#' breaks ties by insertion order, which is what makes the queue fully
#' deterministic: among equal dates the event pushed first comes out first.
#'
#' @section Assumptions:
#' Events may be pushed with a date earlier than the date of the event
#' currently at the front. Nothing here forbids it, and the heap handles it
#' correctly, but the engine's global clock must never move backwards, so
#' the engine is responsible for never scheduling into the past. The queue is a 
#' data structure, the clock is a modelling concept.
#'
#' @param capacity Initial size of the internal arrays. They grow by
#'   doubling when full, so this is a performance hint only. The default
#'   holds a thousand events without any reallocation.
#' @param n_slots Initial number of individuals the queue is prepared to
#'   track. Grows automatically as newborns and immigrants appear.
#'
#' @return An object of class `event_queue`: a list of functions.
#'   \describe{
#'     \item{`push(time_days, slot, to_idx)`}{schedule an event, replacing
#'       any event already pending for `slot`. Returns invisibly.}
#'     \item{`pop_min()`}{remove and return the earliest live event as a
#'       list with `time_days`, `slot` and `to_idx`, or `NULL` when the
#'       queue holds no live events.}
#'     \item{`peek_min()`}{the same, without removing it.}
#'     \item{`cancel(slot)`}{withdraw the event pending for `slot`, if any.
#'       Returns `TRUE` if something was withdrawn.}
#'     \item{`has_pending(slot)`}{whether `slot` currently has a live
#'       event.}
#'     \item{`size()`}{number of live events.}
#'     \item{`n_nodes()`}{array occupancy, live plus stale. Diagnostic.}
#'   }
#'
#' @seealso [MicSim-concepts]
#'
#' @examples
#' q <- new_event_queue()
#' q$push(time_days = 20315, slot = 1L, to_idx = 4L)
#' q$push(time_days = 20100, slot = 2L, to_idx = 7L)
#' q$pop_min()$slot   # 2: its event is the earlier one
#' q$size()           # 1
#'
#' @keywords internal
#' @noRd
new_event_queue <- function(capacity = 1024L, n_slots = 1024L) {

  # ------------------------------------------------------------------
  # Private state.
  #
  # `parent = emptyenv()` means this environment has no parent to fall
  # back on, so a typo cannot silently resolve to a variable somewhere
  # else in the package. It makes the state genuinely private.
  #
  # The heap is held as four parallel vectors rather than as a list of
  # records ("struct of arrays" rather than "array of structs"). Position
  # i of each vector describes the same entry. In R this matters: four
  # atomic vectors are four contiguous blocks of memory, whereas a list of
  # 100,000 three-element lists is 100,000 separate allocations.
  # ------------------------------------------------------------------
  e <- new.env(parent = emptyenv())

  e$time    <- numeric(capacity)   # when the event is due, days since epoch
  e$slot    <- integer(capacity)   # which individual it belongs to
  e$to      <- integer(capacity)   # which state they move into
  e$version <- integer(capacity)   # stamp identifying this entry

  e$n <- 0L                        # entries currently in the heap array

  # pending[s] is the version of slot s's one live event, or 0 for "no
  # live event". An entry in the heap is live if and only if its version
  # equals pending[its slot].
  e$pending <- integer(n_slots)

  e$next_version <- 0L             # monotone counter: version AND tie-break
  e$n_live <- 0L                   # live events, i.e. the queue's true size

  # ------------------------------------------------------------------
  # Growth. Both arrays double, so appending n entries costs O(n) copying
  # in total rather than O(n^2). This is the same amortised-doubling idea
  # the transition store uses.
  # ------------------------------------------------------------------
  grow_heap <- function() {
    old <- length(e$time)
    e$time    <- c(e$time,    numeric(old))
    e$slot    <- c(e$slot,    integer(old))
    e$to      <- c(e$to,      integer(old))
    e$version <- c(e$version, integer(old))
  }

  grow_pending <- function(needed) {
    old <- length(e$pending)
    new_size <- max(needed, old * 2L)
    e$pending <- c(e$pending, integer(new_size - old))
  }

  # ------------------------------------------------------------------
  # Ordering. "a precedes b" means a is due earlier, or they are due at
  # the same moment and a was pushed first. Version numbers are unique, so
  # this is a total order and the queue is fully deterministic. The test
  # is written out inline in the loops below rather than called as a
  # function, because it runs O(log Q) times per push and per pop.
  #
  # WRITES ARE IN PLACE. Every operation below first takes the four heap
  # vectors OUT of the environment into locals, works on the locals, and
  # puts them back. Writing `e$time[i] <- x` directly looks O(1) but is
  # not: the vector is reachable from the environment, so R duplicates it
  # on every such write, which made each push and pop cost O(Q) instead of
  # O(log Q) (measured: 20,000 pushes and pops took 21 s). With the vector
  # held by one local only, R mutates it in place. Same "steal-back" idiom
  # as engine_person.R and engine_store.R.
  # ------------------------------------------------------------------

  # ------------------------------------------------------------------
  # Public interface
  # ------------------------------------------------------------------

  push <- function(time_days, slot, to_idx) {
    slot   <- as.integer(slot)
    t_new  <- as.numeric(time_days)
    d_new  <- as.integer(to_idx)
    if (slot > length(e$pending)) grow_pending(slot)
    if (e$n == length(e$time)) grow_heap()

    tm <- e$time;    e$time    <- NULL
    sl <- e$slot;    e$slot    <- NULL
    to <- e$to;      e$to      <- NULL
    vs <- e$version; e$version <- NULL
    pend <- e$pending; e$pending <- NULL

    # One live event per individual. If this slot already has one, the old
    # entry becomes stale the moment we hand out a new version number;
    # the count of live events is unchanged because one replaced one.
    if (pend[slot] == 0L) e$n_live <- e$n_live + 1L
    e$next_version <- e$next_version + 1L
    v_new <- e$next_version
    pend[slot] <- v_new

    # Sift up: carry the new entry and slide parents down into the gap
    # until its place is found. Same result as swapping, half the writes.
    i <- e$n + 1L
    while (i > 1L) {
      parent <- i %/% 2L
      tp <- tm[parent]
      if (!(t_new < tp || (t_new == tp && v_new < vs[parent]))) break
      tm[i] <- tp
      sl[i] <- sl[parent]
      to[i] <- to[parent]
      vs[i] <- vs[parent]
      i <- parent
    }
    tm[i] <- t_new
    sl[i] <- slot
    to[i] <- d_new
    vs[i] <- v_new

    e$n <- e$n + 1L
    e$time <- tm; e$slot <- sl; e$to <- to; e$version <- vs
    e$pending <- pend
    invisible(NULL)
  }

  # Remove entries from the root until a live one is there, and -- when
  # `consume` -- remove that one too and return it. Stale entries (whose
  # version no longer matches their individual's live version, because
  # the event was cancelled or superseded) are discarded on the way out:
  # that is what "lazy deletion" means.
  front <- function(consume) {
    tm <- e$time;    e$time    <- NULL
    sl <- e$slot;    e$slot    <- NULL
    to <- e$to;      e$to      <- NULL
    vs <- e$version; e$version <- NULL
    pend <- e$pending; e$pending <- NULL
    n <- e$n
    event <- NULL

    repeat {
      if (n == 0L) break
      live <- vs[1L] == pend[sl[1L]]
      if (live) {
        event <- list(time_days = tm[1L], slot = sl[1L], to_idx = to[1L])
        if (!consume) break
        pend[sl[1L]] <- 0L                 # this individual is now unscheduled
        e$n_live <- e$n_live - 1L
      }

      # Remove the root: move the last entry into the gap and let it sink,
      # comparing with the earlier of its two children at each step.
      t_new <- tm[n]; s_new <- sl[n]; d_new <- to[n]; v_new <- vs[n]
      n <- n - 1L
      i <- 1L
      repeat {
        left <- 2L * i
        if (left > n) break
        child <- left
        right <- left + 1L
        if (right <= n &&
              (tm[right] < tm[left] ||
                 (tm[right] == tm[left] && vs[right] < vs[left]))) {
          child <- right
        }
        tc <- tm[child]
        if (!(tc < t_new || (tc == t_new && vs[child] < v_new))) break
        tm[i] <- tc
        sl[i] <- sl[child]
        to[i] <- to[child]
        vs[i] <- vs[child]
        i <- child
      }
      if (n > 0L) {
        tm[i] <- t_new; sl[i] <- s_new; to[i] <- d_new; vs[i] <- v_new
      }

      if (live) break
    }

    e$n <- n
    e$time <- tm; e$slot <- sl; e$to <- to; e$version <- vs
    e$pending <- pend
    event
  }

  peek_min <- function() front(consume = FALSE)

  pop_min <- function() front(consume = TRUE)

  cancel <- function(slot) {
    slot <- as.integer(slot)
    if (slot > length(e$pending) || e$pending[slot] == 0L) {
      return(invisible(FALSE))
    }
    # Forgetting the version is the whole cancellation. The entry itself
    # is left in the heap and will be discarded when it reaches the front.
    pend <- e$pending; e$pending <- NULL
    pend[slot] <- 0L
    e$pending <- pend
    e$n_live <- e$n_live - 1L
    invisible(TRUE)
  }

  has_pending <- function(slot) {
    slot <- as.integer(slot)
    slot <= length(e$pending) && e$pending[slot] != 0L
  }

  structure(
    list(
      push        = push,
      pop_min     = pop_min,
      peek_min    = peek_min,
      cancel      = cancel,
      has_pending = has_pending,
      size        = function() e$n_live,
      n_nodes     = function() e$n,
      state       = e            # exposed for tests and diagnostics only
    ),
    class = "event_queue"
  )
}


#' Print an event queue
#'
#' Reports how many events are live and how much dead weight lazy deletion
#' is currently carrying. A stale count that grows keeps growing over a run
#' means events are being cancelled or rescheduled far more often than they
#' are executed. That costs memory and may be a bug to check.
#'
#' @param x An `event_queue`.
#' @param ... Not used; present for compatibility with the [print()] generic.
#' @return `x`, invisibly. Called for the side effect of printing.
#' @exportS3Method print event_queue
#' @noRd
print.event_queue <- function(x, ...) {
  live <- x$size()
  nodes <- x$n_nodes()
  cat("<event_queue>\n")
  cat("  live events : ", live, "\n", sep = "")
  cat("  heap entries: ", nodes, " (", nodes - live,
      " stale, awaiting lazy removal)\n", sep = "")
  nxt <- x$peek_min()
  if (is.null(nxt)) {
    cat("  next event  : none\n")
  } else {
    cat("  next event  : slot ", nxt$slot, " -> state ", nxt$to_idx,
        " on ", getInDateFormat(nxt$time_days), "\n", sep = "")
  }
  invisible(x)
}
