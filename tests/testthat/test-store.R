# =====================================================================
# test-store.R
# ---------------------------------------------------------------------
# Two things are tested here, and the second matters far more than the
# first.
#
# The STORE is a plain accumulator: append rows, get them back, grow when
# full. Little can go wrong and the tests are correspondingly short.
#
# ASSEMBLE_POP produces the data frame that every downstream consumer of
# this package reads -- both converters, every estimator, every user
# script written since 2011. Its shape is the one part of the redesign
# that may not change, and the mock-up changed it (deviation M1) by
# manufacturing an "entry event" row for every individual. Most of this
# file pins that contract down.
#
# Tests are grouped by what they establish:
#   1. the store accumulates
#   2. growth
#   3. extra columns (linked-lives scaffolding)
#   4. the pop contract: columns
#   5. the pop contract: rows
#   6. the pop contract: values and types
#   7. the pop contract: ordering
#   8. motherID
# =====================================================================


# A two-variable model: states are "m/0", "m/1", "f/0", "f/1", plus
# "dead". The codec enumerates them in that order, so "f/0" is index 3,
# "f/1" is 4 and "dead" is 5.
store_codec <- function() {
  new_state_codec(
    state_space = list(sex = c("m", "f"), fert = c("0", "1")),
    absorbing   = "dead"
  )
}

# Three people:
#   ID 1  a woman in the starting population, who has a child and dies
#   ID 2  a man in the starting population, who does nothing at all
#   ID 3  her child, born during the run, who also does nothing
#
# The point of ID 2 and ID 3 is that they must still appear in `pop`,
# exactly once each, with NA in every transition column.
store_people <- function() {
  p <- new_person_store(n_vars = 2L)
  p$add(1L, getInDays(19800101), 3L, 3L,
        birth_ymd = "19800101")
  p$add(2L, getInDays(19700101), 1L, 1L,
        birth_ymd = "19700101")
  p$add(3L, getInDays(20160701), 1L, 1L,
        birth_ymd = "20160701", mother_id = 1L)
  p
}

store_with_two_events <- function(p) {
  s <- new_transition_store()
  # first birth, then death
  s$append(slot = 1L, from_idx = 3L, to_idx = 4L,
           time_days = getInDays(20160701),
           age_days  = getInDays(20160701) - getInDays(19800101))
  s$append(slot = 1L, from_idx = 4L, to_idx = 5L,
           time_days = getInDays(20200315),
           age_days  = getInDays(20200315) - getInDays(19800101))
  s
}


# ---------------------------------------------------------------------
# 1. The store accumulates
# ---------------------------------------------------------------------

test_that("appended rows come back unchanged", {
  s <- new_transition_store()
  expect_identical(s$size(), 0L)

  s$append(slot = 7L, from_idx = 3L, to_idx = 4L,
           time_days = 17351, age_days = 13699)
  s$append(slot = 2L, from_idx = 1L, to_idx = 5L,
           time_days = 18000.5, age_days = 20000.5)

  cols <- s$as_columns()
  expect_identical(s$size(), 2L)
  expect_identical(cols$slot, c(7L, 2L))
  expect_identical(cols$from, c(3L, 1L))
  expect_identical(cols$to,   c(4L, 5L))
  expect_identical(cols$time_days, c(17351, 18000.5))
  expect_identical(cols$age_days,  c(13699, 20000.5))
})


test_that("an empty store yields empty columns", {
  s <- new_transition_store(capacity = 32L)
  cols <- s$as_columns()
  expect_identical(unname(lengths(cols)), rep(0L, 5L))
})


test_that("the store records rows in the order they happened", {
  # The engine appends in event order, and the assembled output is sorted
  # afterwards, so the store itself must not reorder anything.
  s <- new_transition_store()
  for (i in 1:5) {
    s$append(slot = 1L, from_idx = i, to_idx = i + 1L,
             time_days = 100 * (6 - i), age_days = 1)
  }
  expect_identical(s$as_columns()$from, 1:5)
})


# ---------------------------------------------------------------------
# 2. Growth
# ---------------------------------------------------------------------

test_that("the store grows past its initial capacity", {
  # Doubling growth is what keeps recording E events O(E) rather than the
  # O(E^2) of the legacy per-event rbind() (defect D9).
  s <- new_transition_store(capacity = 4L)
  n <- 500L
  for (i in seq_len(n)) {
    s$append(slot = i, from_idx = 1L, to_idx = 2L,
             time_days = i, age_days = i * 2)
  }
  cols <- s$as_columns()
  expect_identical(s$size(), n)
  expect_identical(cols$slot, seq_len(n))
  expect_identical(cols$time_days, as.numeric(seq_len(n)))
})


# ---------------------------------------------------------------------
# 3. Recorded columns
# ---------------------------------------------------------------------

test_that("the store records exactly five columns per event", {
  s <- new_transition_store()
  expect_identical(names(s$as_columns()),
                   c("slot", "from", "to", "time_days", "age_days"))
})


# ---------------------------------------------------------------------
# 4. The pop contract: columns
# ---------------------------------------------------------------------

test_that("pop has the contract columns in the contract order", {
  # Every downstream consumer reads these by name and several read them by
  # position, so both must hold.
  codec <- store_codec()
  p <- store_people()
  pop <- assemble_pop(store_with_two_events(p), p, codec)

  expect_identical(
    names(pop),
    c("ID", "birthDate", "initState", "From", "To",
      "transitionTime", "transitionAge")
  )
})


# ---------------------------------------------------------------------
# 5. The pop contract: rows
# ---------------------------------------------------------------------

test_that("someone who never transitioned gets exactly one NA row", {
  # This is the heart of the contract. convertToLongFormat() recognises an
  # episode that was never completed by is.na(From), so these rows carry
  # information and must be present, must be single, and must be NA.
  codec <- store_codec()
  p <- store_people()
  pop <- assemble_pop(store_with_two_events(p), p, codec)

  for (id in c(2L, 3L)) {
    rows <- pop[pop$ID == id, , drop = FALSE]
    expect_identical(nrow(rows), 1L)
    expect_true(is.na(rows$From))
    expect_true(is.na(rows$To))
    expect_true(is.na(rows$transitionTime))
    expect_true(is.na(rows$transitionAge))
    # but their identity is still there
    expect_identical(rows$initState, "m/0")
  }
})


test_that("there are no synthetic entry rows", {
  # The mock-up recorded From == To == initState for every individual and
  # every newborn, which changes the row count and the meaning of both
  # columns (deviation M1). The real package never does.
  codec <- store_codec()
  p <- store_people()
  pop <- assemble_pop(store_with_two_events(p), p, codec)

  # 2 transitions + 1 row each for the two who never moved
  expect_identical(nrow(pop), 4L)
  moved <- !is.na(pop$From)
  expect_false(any(pop$From[moved] == pop$To[moved]))
})


test_that("a run in which nobody moves still lists everybody", {
  # MicSim 3.0.0 handled this as a special case and printed a note. Here
  # it falls out of the merge: an empty transition table joined to the
  # population gives one NA row per person.
  codec <- store_codec()
  p <- store_people()
  pop <- assemble_pop(new_transition_store(), p, codec)

  expect_identical(nrow(pop), 3L)
  expect_identical(pop$ID, c(1L, 2L, 3L))
  expect_true(all(is.na(pop$From)))
  expect_true(all(is.na(pop$transitionAge)))
  expect_identical(pop$initState, c("f/0", "m/0", "m/0"))
})


test_that("immigrants and newborns appear like anyone else", {
  # They are appended to the population register during the run, and the
  # register is what pop is built from, so nothing special is needed.
  codec <- store_codec()
  p <- store_people()
  pop <- assemble_pop(store_with_two_events(p), p, codec)

  expect_true(3L %in% pop$ID)
  expect_identical(pop$birthDate[pop$ID == 3L], "20160701")
})


# ---------------------------------------------------------------------
# 6. The pop contract: values and types
# ---------------------------------------------------------------------

test_that("states are decoded back to their labels", {
  codec <- store_codec()
  p <- store_people()
  pop <- assemble_pop(store_with_two_events(p), p, codec)
  moved <- pop[!is.na(pop$From), , drop = FALSE]

  expect_identical(moved$From, c("f/0", "f/1"))
  expect_identical(moved$To,   c("f/1", "dead"))
})


test_that("transitionTime is a character yyyymmdd", {
  # MicSim 3.0.0 produced this column with sprintf(), so it has always
  # been character, and convertToLongFormat() feeds it straight back into
  # getInDays(). See dev/design-notes.md, "Output assembly".
  codec <- store_codec()
  p <- store_people()
  pop <- assemble_pop(store_with_two_events(p), p, codec)
  moved <- pop[!is.na(pop$From), , drop = FALSE]

  expect_type(moved$transitionTime, "character")
  expect_identical(moved$transitionTime, c("20160701", "20200315"))
  expect_type(pop$birthDate, "character")
})


test_that("transitionAge is the age in years, rounded to two decimals", {
  codec <- store_codec()
  p <- store_people()
  pop <- assemble_pop(store_with_two_events(p), p, codec)
  moved <- pop[!is.na(pop$From), , drop = FALSE]

  expected <- round(
    (c(getInDays(20160701), getInDays(20200315)) - getInDays(19800101)) /
      365.25,
    2
  )
  expect_identical(moved$transitionAge, expected)
})


test_that("a fractional event time is reported on the day it falls in", {
  # Event times are continuous. An event at 20160701 plus a third of a day
  # is still an event on 1 July 2016.
  codec <- store_codec()
  p <- store_people()
  s <- new_transition_store()
  s$append(slot = 1L, from_idx = 3L, to_idx = 4L,
           time_days = getInDays(20160701) + 1 / 3,
           age_days  = getInDays(20160701) + 1 / 3 - getInDays(19800101))

  pop <- assemble_pop(s, p, codec)
  expect_identical(pop$transitionTime[!is.na(pop$From)], "20160701")
})


# ---------------------------------------------------------------------
# 7. The pop contract: ordering
# ---------------------------------------------------------------------

test_that("rows are ordered by ID and then by age at transition", {
  codec <- store_codec()
  p <- store_people()
  pop <- assemble_pop(store_with_two_events(p), p, codec)

  expect_identical(pop$ID, c(1L, 1L, 2L, 3L))
  ages <- pop$transitionAge[pop$ID == 1L]
  expect_true(all(diff(ages) > 0))
})


test_that("events sharing a rounded transitionAge stay in chronological order", {
  # transitionAge is rounded to 0.01 years (~3.65 days). Two events of one
  # individual inside the same bin tie on it, and merge() permutes rows
  # within an ID when the transition table is not sorted by ID -- which the
  # store never is, because it is in processing order across everybody.
  # Before the fix, `pop` could list "S -> W" before "W -> S" (defect E7).
  codec <- new_state_codec(list(health = c("W", "S")), absorbing = "dead")
  n <- 200L
  p <- new_person_store(n_vars = 1L)
  for (id in seq_len(n)) {
    p$add(id, getInDays(19700101), 1L, 1L,
          birth_ymd = "19700101")
  }
  birth <- getInDays(19700101)
  s <- new_transition_store()
  # Every individual falls ill and recovers on the same day (same rounded
  # age), and the individuals are processed in an interleaved order, as
  # the engine's global clock produces them.
  set.seed(1L)
  day <- getInDays(20150112) + sample.int(3000L, n, replace = TRUE)
  ord <- order(day)
  for (k in ord) {
    s$append(slot = k, from_idx = 1L, to_idx = 2L,
             time_days = day[k], age_days = day[k] - birth)
  }
  for (k in ord) {
    s$append(slot = k, from_idx = 2L, to_idx = 1L,
             time_days = day[k] + 0.4, age_days = day[k] + 0.4 - birth)
  }

  pop <- assemble_pop(s, p, codec)
  # every individual's rows form a chain: W -> S, then S -> W
  breaks <- 0L
  for (id in seq_len(n)) {
    r <- pop[pop$ID == id, ]
    breaks <- breaks + sum(r$From != c("W", "S") | r$To != c("S", "W"))
  }
  expect_identical(breaks, 0L)
  expect_false(".seq" %in% names(pop))
})


test_that("an individual with no transitions stays with their own ID", {
  # order() compares key by key and sorts NA last WITHIN a tie, so the NA
  # row of a low-numbered individual stays in its own ID group rather than
  # being pushed to the end of the frame. The expression is a literal
  # transcription of MicSim 3.0.0's, so whatever order() does, both
  # engines do the same.
  codec <- store_codec()
  p <- store_people()
  s <- new_transition_store()
  s$append(slot = 2L, from_idx = 1L, to_idx = 5L,
           time_days = getInDays(20180101),
           age_days  = getInDays(20180101) - getInDays(19700101))

  pop <- assemble_pop(s, p, codec)
  expect_identical(pop$ID, c(1L, 2L, 3L))
  expect_true(is.na(pop$From[pop$ID == 1L]))
  expect_identical(pop$To[pop$ID == 2L], "dead")
})


test_that("row names are reset", {
  # After a merge and a reorder they are a meaningless permutation, and
  # the user sees them on every head(pop).
  codec <- store_codec()
  p <- store_people()
  pop <- assemble_pop(store_with_two_events(p), p, codec)
  expect_identical(attr(pop, "row.names"), 1:4)
})


# ---------------------------------------------------------------------
# 8. motherID
# ---------------------------------------------------------------------

test_that("motherID is added only when fertility is modelled", {
  codec <- store_codec()
  p <- store_people()
  s <- store_with_two_events(p)

  without <- assemble_pop(s, p, codec, include_mother = FALSE)
  expect_false("motherID" %in% names(without))

  with_mother <- assemble_pop(s, p, codec, include_mother = TRUE)
  expect_identical(names(with_mother)[8], "motherID")
})


test_that("motherID is NA for everyone not born during the run", {
  codec <- store_codec()
  p <- store_people()
  pop <- assemble_pop(store_with_two_events(p), p, codec,
                      include_mother = TRUE)

  expect_identical(pop$motherID[pop$ID == 3L], 1L)
  expect_true(all(is.na(pop$motherID[pop$ID %in% c(1L, 2L)])))
})


test_that("motherID is repeated across all of an individual's rows", {
  # A newborn who goes on to experience transitions gets one row per
  # transition, and each must carry the mother's identifier. MicSim 3.0.0
  # achieved this with merge(..., all.x = TRUE); here it is a lookup.
  codec <- store_codec()
  p <- store_people()
  s <- new_transition_store()
  s$append(slot = 3L, from_idx = 1L, to_idx = 2L,
           time_days = getInDays(20180101),
           age_days  = getInDays(20180101) - getInDays(20160701))
  s$append(slot = 3L, from_idx = 2L, to_idx = 5L,
           time_days = getInDays(20190101),
           age_days  = getInDays(20190101) - getInDays(20160701))

  pop <- assemble_pop(s, p, codec, include_mother = TRUE)
  expect_identical(pop$motherID[pop$ID == 3L], c(1L, 1L))
})
