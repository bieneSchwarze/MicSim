# =====================================================================
# test-enrolment.R
# ---------------------------------------------------------------------
# School enrolment is kind 4 in the transition taxonomy: the transition
# itself is ordinary, but its DATE is replaced by a calendar convention
# after being drawn and before competing risks are resolved. It is the
# only thing in core MicSim that uses the `on_candidates` hook, which
# makes it the test of that hook the way fertility is the test of
# `on_event`.
#
# The calendar rule is small enough to test on its own, and most of this
# file does exactly that before going near an engine.
#
# Tests are grouped by what they establish:
#   1. the calendar rule
#   2. recognising an enrolment transition
#   3. construction and refusal
#   4. end to end
# =====================================================================


edu_codec <- function() {
  new_state_codec(
    list(sex = c("m", "f"), edu = c("no", "low", "med")),
    absorbing = "dead"
  )
}

#' The user half of school enrolment: a rate that is Inf at age seven.
#'
#' The `Inf` makes the transition forced, so the solver returns that exact
#' date instead of a random one; the subsystem then moves it to the
#' enrolment month. Note the exact equality: MicSim 3.0.0's forced-rate
#' probe evaluates at whole-number ages, which is why this works at all.
#' See dev/design-notes.md, "Truncated arguments".
enrol_rate <- function(age) ifelse(age == 7, Inf, 0)

edu_registry <- function(codec = edu_codec()) {
  new_rate_registry(
    list(
      list(from = "m/no", to = "m/low", rate = enrol_rate),
      list(from = "f/no", to = "f/low", rate = enrol_rate)
    ),
    codec
  )
}

#' Children born in 1996, aged three when the simulation starts in 2000.
#'
#' (In MicSim 3.0.0 this age was the only one that worked: its forced
#' probe evaluated `trunc(age) + k / 365.25`, a whole number only every
#' 1461 days, so a child had to be aged 3 or 7 when scheduled. The test
#' "children of every age enrol at seven" below pins the fix.)
edu_init_pop <- function() {
  data.frame(
    ID        = 1:6,
    birthDate = c("19960115", "19960320", "19960601", "19960915",
                  "19961101", "19961231"),
    initState = rep(c("m/no", "f/no"), 3L),
    stringsAsFactors = FALSE
  )
}

edu_engine <- function(month = 9L, horizon = c(20000101, 20101231), ...) {
  codec <- edu_codec()
  registry <- edu_registry(codec)
  new_engine(codec, registry, edu_init_pop(),
             horizon = horizon, max_age = 120,
             subsystems = list(compiled(
               school_enrolment(move("edu", "no", "low"), month = month),
               codec, registry
             )), ...)
}


# ---------------------------------------------------------------------
# 1. The calendar rule
# ---------------------------------------------------------------------

test_that("a draw at or before the enrolment month snaps to that year", {
  # March 2018, enrolment in September: school starts 1 September 2018.
  expect_identical(enrolment_date(getInDays(20180320), 9L),
                   as.numeric(getInDays(20180901)))
  # 1 September itself stays put.
  expect_identical(enrolment_date(getInDays(20180901), 9L),
                   as.numeric(getInDays(20180901)))
})


test_that("a draw after the enrolment month snaps to the next year", {
  # October 2018: the school year has already begun, so 1 September 2019.
  expect_identical(enrolment_date(getInDays(20181015), 9L),
                   as.numeric(getInDays(20190901)))
  expect_identical(enrolment_date(getInDays(20181231), 9L),
                   as.numeric(getInDays(20190901)))
})


test_that("the rule can move an event backwards within the enrolment month", {
  # 20 September snaps back to 1 September of the same year, nineteen days
  # in the past. MicSim 3.0.0 does this too, and its global clock steps
  # backwards silently when it happens. The subsystem clamps the date at
  # the present instead (see the test of the snap below).
  expect_lt(enrolment_date(getInDays(20180920), 9L), getInDays(20180920))
})


test_that("any enrolment month works", {
  expect_identical(enrolment_date(getInDays(20180101), 1L),
                   as.numeric(getInDays(20180101)))
  expect_identical(enrolment_date(getInDays(20180615), 1L),
                   as.numeric(getInDays(20190101)))
  expect_identical(enrolment_date(getInDays(20181101), 12L),
                   as.numeric(getInDays(20181201)))
})


# ---------------------------------------------------------------------
# 2. Recognising an enrolment transition
# ---------------------------------------------------------------------

test_that("the model names the enrolment transition; no label is reserved", {
  # The same school system with labels of the modeller's own choosing, and
  # an education variable that is not called `edu`. Nothing in the
  # package knows "no" or "low" any more.
  codec <- new_state_codec(
    list(sex = c("male", "female"), school = c("none", "primary", "secondary")),
    absorbing = "deceased"
  )
  registry <- new_rate_registry(
    list(list(move = move("school", "none", "primary"), rate = enrol_rate)),
    codec
  )
  sub <- compiled(school_enrolment(move("school", "none", "primary"), 9L),
                  codec, registry)
  expect_identical(length(sub$claims$timing), 2L)        # one per sex
  claimed <- vapply(registry$edges[sub$claims$timing], `[[`, character(1L),
                    "label")
  expect_setequal(claimed, c("male/none -> male/primary",
                             "female/none -> female/primary"))
})


test_that("renaming every label leaves the simulation unchanged", {
  # The point of removing reserved labels: they carried no information.
  renamed <- c(m = "male", f = "female", no = "none", low = "primary",
               med = "secondary")
  codec2 <- new_state_codec(
    list(sex = c("male", "female"), edu = c("none", "primary", "secondary")),
    absorbing = "deceased"
  )
  registry2 <- new_rate_registry(
    list(list(move = move("edu", "none", "primary"), rate = enrol_rate)),
    codec2
  )
  pop2 <- edu_init_pop()
  pop2$initState <- vapply(strsplit(pop2$initState, "/"), function(p)
    paste(renamed[p], collapse = "/"), character(1L))
  engine2 <- new_engine(codec2, registry2, pop2,
                        horizon = c(20000101, 20101231), max_age = 120,
                        subsystems = list(compiled(
                          school_enrolment(move("edu", "none", "primary"), 9L),
                          codec2, registry2)))

  a <- edu_engine(month = 9L)$run(seed = 8L)$pop
  b <- engine2$run(seed = 8L)$pop
  expect_identical(a$transitionTime, b$transitionTime)
  expect_identical(a$transitionAge, b$transitionAge)
  expect_identical(b$To[!is.na(b$To)],
                   unname(vapply(strsplit(a$To[!is.na(a$To)], "/"), function(p)
                     paste(renamed[p], collapse = "/"), character(1L))))
})


test_that("a model may have two education variables, each with its own rule", {
  codec <- new_state_codec(
    list(school = c("none", "primary"), music = c("untrained", "pupil"))
  )
  registry <- new_rate_registry(
    list(list(move = move("school", "none", "primary"), rate = enrol_rate),
         list(move = move("music", "untrained", "pupil"), rate = enrol_rate)),
    codec
  )
  a <- compiled(school_enrolment(move("school", "none", "primary"), 9L,
                                 name = "school"), codec, registry)
  b <- compiled(school_enrolment(move("music", "untrained", "pupil"), 1L,
                                 name = "music"), codec, registry)
  expect_length(intersect(a$claims$timing, b$claims$timing), 0L)
  pop <- data.frame(ID = 1:4, birthDate = c("19960115", "19960320",
                                           "19960601", "19961101"),
                    initState = "none/untrained", stringsAsFactors = FALSE)
  run <- new_engine(codec, registry, pop, horizon = c(20000101, 20101231),
                    max_age = 120, subsystems = list(a, b))$run(seed = 2L)
  moved <- run$pop[!is.na(run$pop$To), ]
  school_rows <- moved$From != moved$To &
    sub("/.*$", "", moved$From) != sub("/.*$", "", moved$To)
  expect_true(all(getMonth(getInDays(moved$transitionTime[school_rows])) == 9L))
  expect_true(all(getMonth(getInDays(moved$transitionTime[!school_rows])) == 1L))
  expect_identical(run$extensions$school$n_enrolments, 4L)
  expect_identical(run$extensions$music$n_enrolments, 4L)
})


test_that("the classic interface keeps 3.0.0's rule: 'no' to 'low' in one variable", {
  codec <- edu_codec()
  found <- legacy_enrolment_transitions(codec)
  expect_setequal(found, c("m/no -> m/low", "f/no -> f/low"))
  # "m/no" -> "f/low" changes sex as well, so it is not an enrolment.
  for (label in found) {
    ends <- trimws(strsplit(label, "->", fixed = TRUE)[[1L]])
    expect_identical(changed_state_variables(ends[1L], ends[2L]), 2L)
  }
})


test_that("the classic rule is component-wise, not a substring match", {
  # MicSim 3.0.0 deleted the literal substring "no" from one label and
  # "low" from the other and compared the remainders, which misfires on a
  # state value that merely contains those letters.
  plain <- new_state_codec(list(sex = c("m", "f")), absorbing = "dead")
  expect_error(legacy_enrolment_transitions(plain),
               class = "micsim_enrolment_error")
  contains_no <- new_state_codec(
    list(nat = c("Norway", "Italy"), fert = c("0", "1"))
  )
  expect_error(legacy_enrolment_transitions(contains_no),
               class = "micsim_enrolment_error")
})


# ---------------------------------------------------------------------
# 3. Construction and refusal
# ---------------------------------------------------------------------

test_that("the compiled subsystem is an ordinary hook set with on_candidates", {
  codec <- edu_codec()
  sub <- compiled(school_enrolment(move("edu", "no", "low"), 9L), codec,
                  edu_registry(codec))
  expect_s3_class(sub, "micsim_compiled_subsystem")
  inst <- sub$instantiate(list(codec = codec))
  expect_s3_class(inst$hooks, "engine_hooks")
  expect_setequal(inst$hooks$active, c("on_event", "on_candidates"))
  expect_null(sub$moves_events_back)         # the flag no longer exists
  expect_true(sub$parallel_safe)
})


test_that("an impossible enrolment month is refused", {
  e <- move("edu", "no", "low")
  expect_error(school_enrolment(e, 0L), class = "micsim_enrolment_error")
  expect_error(school_enrolment(e, 13L), class = "micsim_enrolment_error")
  expect_error(school_enrolment(e, c(1L, 2L)), class = "micsim_enrolment_error")
  expect_error(school_enrolment(month = 9L), class = "micsim_enrolment_error")
})


test_that("a rule for a transition the model does not declare is refused", {
  codec <- edu_codec()
  registry <- edu_registry(codec)
  expect_error(compiled(school_enrolment(move("edu", "low", "med"), 9L),
                        codec, registry),
               class = "micsim_selector_error")
  expect_error(compiled(school_enrolment(move("edu", "no", "high"), 9L),
                        codec, registry),
               class = "micsim_selector_error")
  expect_error(compiled(school_enrolment(move("school", "no", "low"), 9L),
                        codec, registry),
               class = "micsim_selector_error")
})


# ---------------------------------------------------------------------
# 4. End to end
# ---------------------------------------------------------------------

test_that("every enrolment lands on the first of the enrolment month", {
  # The property the whole feature exists for.
  set.seed(1L)
  run <- edu_engine(month = 9L)$run()
  pop <- run$pop
  enrolments <- pop[!is.na(pop$To) & grepl("low", pop$To, fixed = TRUE), ,
                    drop = FALSE]

  expect_gt(nrow(enrolments), 0L)
  expect_identical(run$extensions$school_enrolment$n_enrolments,
                   nrow(enrolments))
  dates <- getInDays(enrolments$transitionTime)
  expect_true(all(getMonth(dates) == 9L))
  expect_true(all(getDay(dates) == 1L))
})


test_that("enrolment follows the seventh birthday and the school calendar", {
  # Forced at exactly age 7, then moved to 1 September: a birthday on or
  # before the end of September enrols that year, a later one the year
  # after. So ages at enrolment lie between 7 minus a month and 8.
  set.seed(2L)
  pop <- edu_engine(month = 9L)$run()$pop
  enrolments <- pop[!is.na(pop$To) & grepl("low", pop$To, fixed = TRUE), ,
                    drop = FALSE]
  expect_identical(nrow(enrolments), 6L)
  expect_true(all(enrolments$transitionAge > 6.9))
  expect_true(all(enrolments$transitionAge < 8))
  born <- getMonth(getInDays(enrolments$birthDate))
  year <- getYear(getInDays(enrolments$transitionTime))
  expect_identical(year, ifelse(born <= 9L, 2003L, 2004L))
})


test_that("children of every age enrol at seven, newborns included", {
  # MicSim 3.0.0's forced probe fired only for children aged 3 or 7
  # (truncated) when scheduled: 0 of 302 children born during a test run
  # ever enrolled. Here every child does, whatever its age at start.
  codec <- edu_codec()
  registry <- edu_registry(codec)
  ages <- data.frame(ID = 1:8,
                     birthDate = c("19990101", "19990701", "19980315",
                                   "19970820", "19960101", "19950601",
                                   "19940303", "19931010"),
                     initState = rep(c("m/no", "f/no"), 4L),
                     stringsAsFactors = FALSE)
  run <- new_engine(codec, registry, ages, horizon = c(20000101, 20101231),
                    max_age = 120,
                    subsystems = list(compiled(
                      school_enrolment(move("edu", "no", "low"), 9L),
                      codec, registry)))$run(seed = 1L)
  enrolled <- run$pop[!is.na(run$pop$To), ]
  expect_setequal(enrolled$ID, ages$ID)
  expect_true(all(enrolled$transitionAge > 6.9 & enrolled$transitionAge < 8))
})


test_that("a different enrolment month moves every enrolment", {
  set.seed(3L)
  january <- edu_engine(month = 1L)$run()$pop
  enrolments <- january[!is.na(january$To) &
                          grepl("low", january$To, fixed = TRUE), ,
                        drop = FALSE]
  expect_gt(nrow(enrolments), 0L)
  dates <- getInDays(enrolments$transitionTime)
  expect_true(all(getMonth(dates) == 1L))
  expect_true(all(getDay(dates) == 1L))
})


test_that("the snap never moves the clock backwards", {
  # MicSim 3.0.0 moved a draw later in the enrolment month back to the
  # 1st, sometimes to before the present. Now a candidate is never placed
  # before now, so the engine's assertion (an error on any reversal)
  # holds and the run completes.
  set.seed(4L)
  d <- edu_engine(month = 9L)$run()$diagnostics
  expect_null(d$n_clock_reversals)  # the counter no longer exists
  expect_gt(d$n_forced, 0L)        # the Inf rate did fire

  # A child already seven, scheduled on 15 September: the rule says
  # 1 September, which is past, so the enrolment happens now.
  late <- data.frame(ID = 1L, birthDate = "19930301", initState = "f/no",
                     stringsAsFactors = FALSE)
  codec <- edu_codec()
  registry <- edu_registry(codec)
  pop <- new_engine(codec, registry, late, horizon = c(20000915, 20051231),
                    max_age = 120,
                    subsystems = list(compiled(
                      school_enrolment(move("edu", "no", "low"), 9L),
                      codec, registry)))$run(seed = 1L)$pop
  expect_identical(pop$transitionTime, "20000915")
})


test_that("other risks compete over the time an enrolment is postponed", {
  # A child whose seventh birthday is in October enrols the next
  # September. With a certain death in between, the death must win; in
  # MicSim 3.0.0 the enrolment had already won when its date was moved.
  codec <- edu_codec()
  doomed <- function(age) ifelse(age >= 7.2, 1e4, 0)
  registry <- new_rate_registry(
    list(list(from = "f/no", to = "f/low", rate = enrol_rate),
         list(from = NA, to = "dead", rate = doomed)),
    codec
  )
  child <- data.frame(ID = 1L, birthDate = "19961015", initState = "f/no",
                      stringsAsFactors = FALSE)
  pop <- new_engine(codec, registry, child, horizon = c(20000101, 20101231),
                    max_age = 120,
                    subsystems = list(compiled(
                      school_enrolment(move("edu", "no", "low"), 9L),
                      codec, registry)))$run(seed = 3L)$pop
  expect_identical(pop$To, "dead")
  expect_lt(pop$transitionAge, 7.3)
})


test_that("the snap changes the date and never the destination", {
  set.seed(5L)
  pop <- edu_engine(month = 9L)$run()$pop
  moved <- pop[!is.na(pop$From), , drop = FALSE]
  expect_gt(nrow(moved), 0L)
  # the only transition in this model is no -> low, within one sex
  expect_true(all(sub("/.*$", "", moved$From) == sub("/.*$", "", moved$To)))
  expect_true(all(grepl("/no$", moved$From)))
  expect_true(all(grepl("/low$", moved$To)))
})


test_that("the run is reproducible with enrolment on", {
  engine <- edu_engine(month = 9L)
  expect_identical(engine$run(seed = 41L)$pop, engine$run(seed = 41L)$pop)
})
