# =====================================================================
# test-queue.R
# ---------------------------------------------------------------------
# The event queue is a pure data structure, so it can be checked against
# an obviously correct reference rather than against itself: whatever
# sequence of pushes we make, popping everything out must reproduce
# sort(). That is the strongest kind of test available in this package and
# it costs nothing to run.
#
# Tests are grouped by what they establish:
#   1. ordering                (the heap property holds)
#   2. determinism             (ties break by insertion order)
#   3. growth                  (arrays double correctly)
#   4. one live event per slot (rescheduling supersedes)
#   5. cancellation            (lazy deletion, needed by linked lives)
#   6. bookkeeping             (size vs n_nodes)
#   7. empty-queue behaviour
#   8. a model-based stress test
# =====================================================================


# ---------------------------------------------------------------------
# 1. Ordering
# ---------------------------------------------------------------------

test_that("popping everything reproduces sort()", {
  # The defining property. Slots are distinct because the queue holds at
  # most one live event per individual; using one slot per event keeps
  # this test purely about ordering.
  set.seed(11L)
  n <- 200L
  times <- round(runif(n, 18000, 22000), 3)

  q <- new_event_queue(capacity = 8L)
  for (i in seq_len(n)) {
    q$push(time_days = times[i], slot = i, to_idx = i %% 7L + 1L)
  }

  popped <- numeric(n)
  for (i in seq_len(n)) {
    popped[i] <- q$pop_min()$time_days
  }

  expect_identical(popped, sort(times))
  expect_null(q$pop_min())
})


test_that("the payload travels with the key", {
  # An entry's slot and destination must still be attached to its own date
  # after every swap the heap performs.
  q <- new_event_queue()
  q$push(time_days = 30, slot = 3L, to_idx = 33L)
  q$push(time_days = 10, slot = 1L, to_idx = 11L)
  q$push(time_days = 20, slot = 2L, to_idx = 22L)

  first  <- q$pop_min()
  second <- q$pop_min()
  third  <- q$pop_min()

  expect_identical(c(first$slot, second$slot, third$slot), c(1L, 2L, 3L))
  expect_identical(c(first$to_idx, second$to_idx, third$to_idx),
                   c(11L, 22L, 33L))
  expect_identical(c(first$time_days, second$time_days, third$time_days),
                   c(10, 20, 30))
})


test_that("events may be pushed out of chronological order", {
  # The heap does not care; only the engine's global clock does. Pushing a
  # date earlier than everything already queued must still come out first.
  q <- new_event_queue()
  q$push(time_days = 500, slot = 1L, to_idx = 1L)
  q$push(time_days = 400, slot = 2L, to_idx = 1L)
  q$push(time_days = 100, slot = 3L, to_idx = 1L)
  expect_identical(q$pop_min()$slot, 3L)
})


# ---------------------------------------------------------------------
# 2. Determinism
# ---------------------------------------------------------------------

test_that("simultaneous events come out in the order they were pushed", {
  # Without this a fixed seed would not give a fixed result. MicSim 3.0.0
  # got the same behaviour from order(), which is stable, so this also
  # keeps the two engines comparable.
  q <- new_event_queue()
  for (s in 1:5) {
    q$push(time_days = 12345, slot = s, to_idx = 1L)
  }
  order_out <- vapply(1:5, function(i) q$pop_min()$slot, integer(1L))
  expect_identical(order_out, 1:5)
})


test_that("ties are broken by insertion order even after reheaping", {
  # Interleave equal and unequal dates so that the tied entries are
  # actually moved around inside the array before they surface.
  q <- new_event_queue(capacity = 4L)
  q$push(time_days = 100, slot = 1L, to_idx = 1L)
  q$push(time_days = 50,  slot = 2L, to_idx = 1L)
  q$push(time_days = 100, slot = 3L, to_idx = 1L)
  q$push(time_days = 10,  slot = 4L, to_idx = 1L)
  q$push(time_days = 100, slot = 5L, to_idx = 1L)

  order_out <- vapply(1:5, function(i) q$pop_min()$slot, integer(1L))
  expect_identical(order_out, c(4L, 2L, 1L, 3L, 5L))
})


# ---------------------------------------------------------------------
# 3. Growth
# ---------------------------------------------------------------------

test_that("the heap grows past its initial capacity", {
  # Doubling growth is what keeps total insertion cost linear; the test
  # only checks correctness, since asserting on timings is unreliable.
  q <- new_event_queue(capacity = 2L)
  n <- 100L
  for (i in seq_len(n)) {
    q$push(time_days = n - i + 1, slot = i, to_idx = 1L)
  }
  expect_identical(q$size(), n)
  popped <- vapply(seq_len(n), function(i) q$pop_min()$time_days, numeric(1L))
  expect_identical(popped, as.numeric(seq_len(n)))
})


test_that("the per-individual table grows when a new slot appears", {
  # Newborns and immigrants arrive during the run with slots beyond
  # anything allocated at the start.
  q <- new_event_queue(capacity = 4L, n_slots = 4L)
  q$push(time_days = 10, slot = 1L, to_idx = 1L)
  q$push(time_days = 5,  slot = 5000L, to_idx = 2L)

  expect_true(q$has_pending(5000L))
  expect_identical(q$pop_min()$slot, 5000L)
})


# ---------------------------------------------------------------------
# 4. One live event per individual
# ---------------------------------------------------------------------

test_that("pushing again for the same individual supersedes the old event", {
  # This is rescheduling. The engine reschedules whenever an individual's
  # circumstances change before their pending event fires.
  q <- new_event_queue()
  q$push(time_days = 900, slot = 1L, to_idx = 4L)
  q$push(time_days = 100, slot = 1L, to_idx = 9L)

  expect_identical(q$size(), 1L)

  event <- q$pop_min()
  expect_identical(event$time_days, 100)
  expect_identical(event$to_idx, 9L)
  expect_null(q$pop_min())          # the superseded entry never surfaces
})


test_that("a later reschedule also supersedes an earlier one", {
  # The superseded entry is the one at the FRONT here, so the stale check
  # has to fire on the first pop rather than on a later one.
  q <- new_event_queue()
  q$push(time_days = 100, slot = 1L, to_idx = 4L)
  q$push(time_days = 900, slot = 1L, to_idx = 9L)

  expect_identical(q$size(), 1L)
  event <- q$pop_min()
  expect_identical(event$time_days, 900)
  expect_null(q$pop_min())
})


# ---------------------------------------------------------------------
# 5. Cancellation
# ---------------------------------------------------------------------

test_that("a cancelled event never surfaces", {
  # Required by linked lives: on separation, the partner's pending event
  # is withdrawn and replaced by the induced separation transition.
  q <- new_event_queue()
  q$push(time_days = 10, slot = 1L, to_idx = 1L)
  q$push(time_days = 20, slot = 2L, to_idx = 1L)
  q$push(time_days = 30, slot = 3L, to_idx = 1L)

  expect_true(q$cancel(2L))
  expect_false(q$has_pending(2L))
  expect_identical(q$size(), 2L)

  expect_identical(q$pop_min()$slot, 1L)
  expect_identical(q$pop_min()$slot, 3L)
  expect_null(q$pop_min())
})


test_that("cancelling the front event works", {
  q <- new_event_queue()
  q$push(time_days = 10, slot = 1L, to_idx = 1L)
  q$push(time_days = 20, slot = 2L, to_idx = 1L)

  q$cancel(1L)
  expect_identical(q$pop_min()$slot, 2L)
  expect_null(q$pop_min())
})


test_that("cancelling a slot with no pending event is a no-op", {
  q <- new_event_queue()
  q$push(time_days = 10, slot = 1L, to_idx = 1L)

  expect_false(q$cancel(2L))        # never scheduled
  expect_false(q$cancel(99999L))    # beyond the allocated table
  expect_identical(q$size(), 1L)
})


test_that("a cancelled individual can be scheduled again", {
  q <- new_event_queue()
  q$push(time_days = 10, slot = 1L, to_idx = 1L)
  q$cancel(1L)
  q$push(time_days = 50, slot = 1L, to_idx = 7L)

  expect_identical(q$size(), 1L)
  event <- q$pop_min()
  expect_identical(event$time_days, 50)
  expect_identical(event$to_idx, 7L)
})


test_that("popping clears the pending flag", {
  # After an event fires the individual is unscheduled until the engine
  # schedules their next one. Without this, a reschedule would be treated
  # as a supersession and the live count would drift.
  q <- new_event_queue()
  q$push(time_days = 10, slot = 1L, to_idx = 1L)
  q$pop_min()

  expect_false(q$has_pending(1L))
  expect_identical(q$size(), 0L)

  q$push(time_days = 20, slot = 1L, to_idx = 1L)
  expect_identical(q$size(), 1L)
})


# ---------------------------------------------------------------------
# 6. Bookkeeping
# ---------------------------------------------------------------------

test_that("size counts live events and n_nodes counts array occupancy", {
  # The gap between the two is the dead weight lazy deletion carries. It
  # is bounded by the number of pushes, but it should be visible.
  q <- new_event_queue()
  q$push(time_days = 10, slot = 1L, to_idx = 1L)
  q$push(time_days = 20, slot = 2L, to_idx = 1L)
  q$push(time_days = 30, slot = 1L, to_idx = 1L)   # supersedes the first

  expect_identical(q$size(), 2L)
  expect_identical(q$n_nodes(), 3L)                # one entry is stale

  q$pop_min()                                      # slot 2 at day 20
  q$pop_min()                                      # slot 1 at day 30
  expect_identical(q$size(), 0L)
  expect_identical(q$n_nodes(), 0L)                # the stale one was dropped
})


test_that("peek_min does not remove the event", {
  q <- new_event_queue()
  q$push(time_days = 10, slot = 1L, to_idx = 5L)
  q$push(time_days = 20, slot = 2L, to_idx = 6L)

  seen <- q$peek_min()
  expect_identical(seen$slot, 1L)
  expect_identical(q$size(), 2L)
  expect_identical(q$pop_min()$slot, 1L)
})


test_that("peek_min skips stale entries", {
  q <- new_event_queue()
  q$push(time_days = 10, slot = 1L, to_idx = 1L)
  q$push(time_days = 20, slot = 2L, to_idx = 1L)
  q$cancel(1L)

  expect_identical(q$peek_min()$slot, 2L)
})


# ---------------------------------------------------------------------
# 7. Empty-queue behaviour
# ---------------------------------------------------------------------

test_that("an empty queue returns NULL rather than erroring", {
  # The engine's main loop stops on NULL, so this is the normal end of a
  # run and must not be an exceptional condition.
  q <- new_event_queue()
  expect_null(q$pop_min())
  expect_null(q$peek_min())
  expect_identical(q$size(), 0L)
  expect_identical(q$n_nodes(), 0L)
})


test_that("a queue emptied by cancellation returns NULL", {
  q <- new_event_queue()
  q$push(time_days = 10, slot = 1L, to_idx = 1L)
  q$cancel(1L)
  expect_null(q$pop_min())
  expect_identical(q$size(), 0L)
})


# ---------------------------------------------------------------------
# 8. Model-based stress test
# ---------------------------------------------------------------------

test_that("interleaved pushes, pops and cancellations match a reference model", {
  # A plain R model of what the queue should contain, run alongside the
  # real thing over a long random sequence of operations. This catches the
  # class of bug that hand-written cases miss: a sift that is correct in
  # isolation but leaves the array inconsistent after a particular
  # interleaving.
  skip_on_cran()
  set.seed(2026L)

  q <- new_event_queue(capacity = 4L, n_slots = 4L)

  # The model: for each slot, the date and insertion order of its one live
  # event. A slot missing from `model_time` has no live event.
  model_time    <- numeric(0)
  model_version <- numeric(0)
  version <- 0L

  model_min_slot <- function() {
    if (length(model_time) == 0L) return(NA_character_)
    earliest <- min(model_time)
    tied <- names(model_time)[model_time == earliest]
    tied[which.min(model_version[tied])]
  }

  for (step in seq_len(600L)) {
    action <- sample(c("push", "pop", "cancel"), 1L,
                     prob = c(0.55, 0.35, 0.10))
    slot <- sample.int(60L, 1L)
    key <- as.character(slot)

    if (action == "push") {
      when <- round(runif(1L, 0, 1000), 2L)
      version <- version + 1L
      q$push(time_days = when, slot = slot, to_idx = 1L)
      model_time[key]    <- when
      model_version[key] <- version

    } else if (action == "pop") {
      expected <- model_min_slot()
      got <- q$pop_min()
      if (is.na(expected)) {
        expect_null(got)
      } else {
        expect_identical(as.character(got$slot), expected)
        expect_identical(got$time_days, unname(model_time[expected]))
        model_time    <- model_time[names(model_time) != expected]
        model_version <- model_version[names(model_version) != expected]
      }

    } else {
      q$cancel(slot)
      model_time    <- model_time[names(model_time) != key]
      model_version <- model_version[names(model_version) != key]
    }

    expect_identical(q$size(), length(model_time))
  }
})
