# =====================================================================
# test-person.R
# ---------------------------------------------------------------------
# The person store answers the three questions the engine asks on every
# event. Two of them are trivial subtractions; the third, DURATION, is the
# most intricate piece of logic in MicSim 3.0.0 and the one the mock-up
# got wrong (deviation M2). Most of this file is about that.
#
# Tests are grouped by what they establish:
#   1. registration and slots
#   2. identifier lookup
#   3. age
#   4. duration clocks                 <- the substance
#   5. equivalence with MicSim 3.0.0's fallback
#   6. growth
#   7. identifiers for newborns
#   8. output columns
# =====================================================================


# A convenient fixture: a two-variable model, e.g. sex and parity, so that
# states look like "f/0" and there are two duration clocks per person.
jan_1980 <- getInDays(19800101)
jan_2014 <- getInDays(20140101)
jul_2016 <- getInDays(20160701)
jul_2017 <- getInDays(20170701)


# ---------------------------------------------------------------------
# 1. Registration and slots
# ---------------------------------------------------------------------

test_that("slots are handed out consecutively from one", {
  # A slot is a row number, and the event queue stores slots rather than
  # identifiers precisely so that a pop needs no lookup.
  p <- new_person_store(n_vars = 2L)
  expect_identical(p$size(), 0L)

  expect_identical(p$add(1L, jan_1980, 3L, 3L), 1L)
  expect_identical(p$add(2L, jan_1980, 4L, 4L), 2L)
  expect_identical(p$add(7L, jan_1980, 1L, 1L), 3L)
  expect_identical(p$size(), 3L)
})


test_that("attributes come back exactly as they went in", {
  p <- new_person_store(n_vars = 2L)
  s <- p$add(id = 42L, birth_days = jan_1980, state_idx = 3L,
             init_state_idx = 3L,
             birth_ymd = "19800101", is_migrant = TRUE, mother_id = 7L)

  expect_identical(p$id_of(s), 42L)
  expect_identical(p$birth_days_of(s), as.numeric(jan_1980))
  expect_identical(p$birth_ymd_of(s), "19800101")
  expect_identical(p$state_of(s), 3L)
  expect_identical(p$init_state_of(s), 3L)
  expect_true(p$is_migrant(s))
  expect_identical(p$mother_id_of(s), 7L)
})


test_that("the birth date is derived when it is not supplied", {
  # Members of the starting population bring their birthDate as the user
  # wrote it, and the output contract requires it back unchanged.
  # Newborns have no such string, so it is generated.
  p <- new_person_store(n_vars = 1L)
  s <- p$add(1L, jan_1980, 1L, 1L)
  expect_identical(p$birth_ymd_of(s), "19800101")
})


test_that("a model must have at least one state variable", {
  expect_error(new_person_store(n_vars = 0L),
               class = "micsim_person_store_error")
})


# ---------------------------------------------------------------------
# 2. Identifier lookup
# ---------------------------------------------------------------------

test_that("slot_of inverts the identifier", {
  p <- new_person_store(n_vars = 1L)
  p$add(11L, jan_1980, 1L, 1L)
  p$add(22L, jan_1980, 1L, 1L)

  expect_identical(p$slot_of(11L), 1L)
  expect_identical(p$slot_of(22L), 2L)
})


test_that("an unknown identifier gives NA rather than an error", {
  p <- new_person_store(n_vars = 1L)
  p$add(1L, jan_1980, 1L, 1L)
  expect_true(is.na(p$slot_of(999L)))
})


test_that("duplicate identifiers are rejected", {
  # Identifiers must be unique across the starting population, the
  # immigrant population, and everyone born during the run. MicSim 3.0.0
  # only warned about collisions, and only in the parallel path.
  p <- new_person_store(n_vars = 1L)
  p$add(5L, jan_1980, 1L, 1L)
  expect_error(p$add(5L, jan_1980, 1L, 1L),
               class = "micsim_person_store_error")
})


test_that("large and sparse identifiers cost nothing", {
  # The lookup is a hashed environment, not a named vector, so an
  # identifier of one million in a population of three is not a problem.
  # This is the mock-up's M4 defect: it indexed named atomic vectors by
  # as.character(id), and R scans the names of an atomic vector linearly.
  p <- new_person_store(n_vars = 1L)
  p$add(1L,       jan_1980, 1L, 1L)
  p$add(1000000L, jan_1980, 1L, 1L)
  p$add(50L,      jan_1980, 1L, 1L)

  expect_identical(p$slot_of(1000000L), 2L)
  expect_identical(p$size(), 3L)
})


# ---------------------------------------------------------------------
# 3. Age
# ---------------------------------------------------------------------

test_that("age is the exact day count since birth", {
  p <- new_person_store(n_vars = 1L)
  s <- p$add(1L, jan_1980, 1L, 1L)

  expect_identical(p$age_days(s, jan_2014),
                   as.numeric(jan_2014 - jan_1980))
  # sanity in the unit a demographer thinks in
  expect_equal(p$age_days(s, jan_2014) / 365.25, 34, tolerance = 0.01)
})


test_that("age is computed at a fractional event time", {
  # Event dates are continuous; birth dates are whole days.
  p <- new_person_store(n_vars = 1L)
  s <- p$add(1L, jan_1980, 1L, 1L)
  expect_identical(p$age_days(s, jan_2014 + 0.5),
                   as.numeric(jan_2014 - jan_1980) + 0.5)
})


# ---------------------------------------------------------------------
# 4. Duration clocks
# ---------------------------------------------------------------------

test_that("duration is measured from the date the variable last changed", {
  # A woman childless at simulation start has her first birth on
  # 1 July 2016. One year later, the duration that a second-birth rate
  # sees must be one year, not the time since some other transition.
  p <- new_person_store(n_vars = 2L)          # sex, parity
  s <- p$add(1L, jan_1980, 3L, 3L)

  p$move_to(s, to_idx = 4L, now_days = jul_2016, changed_vars = 2L)

  expect_identical(p$duration_days(s, var_idx = 2L, now_days = jul_2017),
                   as.numeric(jul_2017 - jul_2016))
  expect_equal(p$duration_days(s, 2L, jul_2017) / 365.25, 1,
               tolerance = 0.01)
})


test_that("an untouched variable keeps measuring from birth", {
  # Sex never changes, so its clock never moves. A rate that depends on
  # duration in the sex variable would see the person's age, which is the
  # correct reading of "time since this variable last took a new value".
  p <- new_person_store(n_vars = 2L)
  s <- p$add(1L, jan_1980, 3L, 3L)
  p$move_to(s, to_idx = 4L, now_days = jul_2016, changed_vars = 2L)

  expect_identical(p$duration_days(s, var_idx = 1L, now_days = jul_2017),
                   p$age_days(s, jul_2017))
})


test_that("the clocks of different variables are independent", {
  # Two transitions, each touching a different variable. Neither may
  # disturb the other's clock. This is precisely what the mock-up got
  # wrong: it kept ONE entry time per person, so any transition reset the
  # duration of every rate (deviation M2).
  p <- new_person_store(n_vars = 3L)
  s <- p$add(1L, jan_1980, 1L, 1L)

  p$move_to(s, to_idx = 2L, now_days = jan_2014 + 100, changed_vars = 2L)
  p$move_to(s, to_idx = 3L, now_days = jan_2014 + 700, changed_vars = 3L)

  now <- jan_2014 + 1000
  expect_identical(p$duration_days(s, 2L, now), 900)
  expect_identical(p$duration_days(s, 3L, now), 300)
  expect_identical(p$duration_days(s, 1L, now), p$age_days(s, now))
})


test_that("a transition changing several variables resets all their clocks", {
  # Permitted as long as no duration-dependent rate is attached to it; the
  # registry rejects that combination. Where it is permitted, every clock
  # the transition disturbs must move.
  p <- new_person_store(n_vars = 3L)
  s <- p$add(1L, jan_1980, 1L, 1L)

  p$move_to(s, to_idx = 5L, now_days = jan_2014 + 200,
            changed_vars = c(1L, 3L))

  now <- jan_2014 + 500
  expect_identical(p$duration_days(s, 1L, now), 300)
  expect_identical(p$duration_days(s, 3L, now), 300)
  expect_identical(p$duration_days(s, 2L, now), p$age_days(s, now))
})


test_that("a move into an absorbing state disturbs no clock", {
  # An absorbing destination is an atomic label such as "dead", so there
  # are no components to compare and `changed_vars` is empty. The
  # individual never schedules again, so nothing needs resetting.
  p <- new_person_store(n_vars = 2L)
  s <- p$add(1L, jan_1980, 3L, 3L)
  before <- p$cov_change_of(s)

  p$move_to(s, to_idx = 9L, now_days = jul_2016, changed_vars = integer(0))

  expect_identical(p$state_of(s), 9L)
  expect_identical(p$cov_change_of(s), before)
})


test_that("a transition with no single affected variable gives Inf", {
  # MicSim 3.0.0 initialised durSinceLastCovCh to Inf and only replaced it
  # when exactly one covariate changed. That is why a mortality rate may
  # declare a `duration` argument and simply ignore it.
  p <- new_person_store(n_vars = 2L)
  s <- p$add(1L, jan_1980, 3L, 3L)
  expect_identical(p$duration_days(s, var_idx = NA_integer_,
                                   now_days = jul_2017), Inf)
})


test_that("move_to updates the state as well as the clocks", {
  p <- new_person_store(n_vars = 2L)
  s <- p$add(1L, jan_1980, 3L, 3L)
  p$move_to(s, to_idx = 4L, now_days = jul_2016, changed_vars = 2L)
  expect_identical(p$state_of(s), 4L)
})


# ---------------------------------------------------------------------
# 5. Equivalence with MicSim 3.0.0's fallback
# ---------------------------------------------------------------------

test_that("before any observed change, duration equals age", {
  # This is the legacy fallback, reproduced by construction rather than by
  # a special case. MicSim 3.0.0 reconstructed the individual's history,
  # found the entry date unknown for anyone not born during the run, and
  # set `durSinceLastCovCh <- currAge`. Here the clock simply still sits
  # at the date of birth, which gives the same number.
  #
  # The three cases the legacy engine distinguished, all agreeing:
  p <- new_person_store(n_vars = 2L)

  from_start <- p$add(1L, jan_1980, 3L, 3L)
  immigrant  <- p$add(2L, jan_1980, 3L, 3L, is_migrant = TRUE)
  newborn    <- p$add(3L, jul_2016, 3L, 3L, mother_id = 1L)

  for (s in c(from_start, immigrant, newborn)) {
    expect_identical(p$duration_days(s, 1L, jul_2017),
                     p$age_days(s, jul_2017))
  }

  # And specifically for the immigrant: measured from BIRTH, not from the
  # immigration date. That is the legacy behaviour, and it is a stated
  # modelling assumption rather than an oversight.
  expect_false(identical(p$duration_days(immigrant, 1L, jul_2017),
                         as.numeric(jul_2017 - jul_2016)))
})


test_that("a newborn starts every clock at zero", {
  p <- new_person_store(n_vars = 2L)
  s <- p$add(3L, jul_2016, 3L, 3L, mother_id = 1L)
  expect_identical(p$duration_days(s, 1L, jul_2016), 0)
  expect_identical(p$duration_days(s, 2L, jul_2016), 0)
})


# ---------------------------------------------------------------------
# 6. Growth
# ---------------------------------------------------------------------

test_that("the store grows past its initial capacity", {
  # Doubling growth is what keeps registering N people O(N) rather than
  # the O(N^2) of the legacy per-newborn rbind() (defect D9).
  p <- new_person_store(n_vars = 2L, capacity = 4L)
  n <- 300L
  for (i in seq_len(n)) {
    p$add(i, jan_1980 + i, 1L, 1L)
  }

  expect_identical(p$size(), n)
  expect_identical(p$slot_of(n), n)
  expect_identical(p$birth_days_of(n), as.numeric(jan_1980 + n))
  # the clocks survived the reallocation
  expect_identical(p$duration_days(n, 1L, jan_2014),
                   p$age_days(n, jan_2014))
})


# ---------------------------------------------------------------------
# 7. Identifiers for newborns
# ---------------------------------------------------------------------

test_that("next_id continues past the largest identifier seen", {
  # MicSim 3.0.0 kept a `maxId` counter in a superassigned variable and
  # incremented it on every birth.
  p <- new_person_store(n_vars = 1L)
  p$add(1L,  jan_1980, 1L, 1L)
  p$add(99L, jan_1980, 1L, 1L)

  expect_identical(p$max_id(), 99L)
  expect_identical(p$next_id(), 100L)
  expect_identical(p$next_id(), 101L)
})


test_that("an identifier handed out by next_id can be registered", {
  p <- new_person_store(n_vars = 1L)
  p$add(1L, jan_1980, 1L, 1L)
  new_id <- p$next_id()
  s <- p$add(new_id, jul_2016, 1L, 1L, mother_id = 1L)
  expect_identical(p$slot_of(new_id), s)
})


# ---------------------------------------------------------------------
# 8. Output columns
# ---------------------------------------------------------------------

test_that("as_columns returns only the rows in use", {
  # The arrays are over-allocated; the output must not carry the padding.
  p <- new_person_store(n_vars = 1L, capacity = 64L)
  p$add(1L, jan_1980, 1L, 1L)
  p$add(2L, jan_1980, 2L, 2L, is_migrant = TRUE)

  cols <- p$as_columns()
  expect_identical(lengths(cols), c(ID = 2L, birthDate = 2L, birthDays = 2L,
                                    initState = 2L, state = 2L,
                                    isMigrant = 2L,
                                    motherID = 2L))
  expect_identical(cols$ID, c(1L, 2L))
  expect_identical(cols$isMigrant, c(FALSE, TRUE))
})


test_that("an empty store produces empty columns", {
  p <- new_person_store(n_vars = 1L)
  cols <- p$as_columns()
  expect_identical(unname(lengths(cols)), rep(0L, length(cols)))
})
