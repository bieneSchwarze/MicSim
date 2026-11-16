# =====================================================================
# test-engine.R
# ---------------------------------------------------------------------
# The first end-to-end tests in the package's history: a model goes in, a
# simulated population comes out. MicSim 3.0.0 had none, because its
# engine could not be called without its input-validation and
# output-assembly code, and none of the three could be called without the
# calling frame it reached into.
#
# Three kinds of test here, in increasing order of what they prove:
#
#   STRUCTURAL   the output has the shape the contract demands
#   BEHAVIOURAL  the loop does what the algorithm says (clock monotone,
#                absorbing states absorb, the horizon closes)
#   STATISTICAL  the simulated population reproduces a distribution we
#                know in closed form
#
# The statistical tests are the point. Everything else in this package can
# be checked in isolation; only here can we ask whether the machine as a
# whole draws from the right distribution.
#
# Tests are grouped by what they establish:
#   1. the smallest possible run
#   2. the output contract
#   3. the loop's invariants
#   4. maxAge and the horizon
#   5. determinism and the random stream
#   6. statistical correctness
#   7. diagnostics
#   8. refusal of the unimplemented stages
#   9. input validation
# =====================================================================


# ---------------------------------------------------------------------
# Fixtures
#
# A one-variable model -- sex, plus an absorbing "dead" -- so states are
# "m", "f" and "dead" with indices 1, 2 and 3. Mortality is the only
# transition, which makes every individual's life course at most one event
# long and therefore easy to reason about exactly.
# ---------------------------------------------------------------------

mortality_codec <- function() {
  new_state_codec(list(sex = c("m", "f")), absorbing = "dead")
}

#' A registry with a constant mortality hazard
#'
#' A constant hazard is the case whose waiting-time distribution is known
#' in closed form, which is what lets section 6 test the whole engine
#' against Exponential(rate) rather than against itself.
#'
#' The probe grid matters for speed, not correctness: without it every
#' edge is conservatively treated as possibly forced and the solver builds
#' a daily grid on every single scheduling call.
constant_registry <- function(rate_value = 0.05, codec = mortality_codec()) {
  constant <- function(age, calTime) rep(rate_value, length(age))
  new_rate_registry(
    transitions = list(list(from = NA, to = "dead", rate = constant)),
    codec = codec,
    probe = micsim_rate_probe(ages = 0:120, cal_years = 1990:2060,
                              durations = 0:1)
  )
}

#' A cohort born on the same day, half male and half female
same_age_cohort <- function(n, birth_ymd = "19900101") {
  data.frame(
    ID        = seq_len(n),
    birthDate = rep(birth_ymd, n),
    initState = rep(c("m", "f"), length.out = n),
    stringsAsFactors = FALSE
  )
}

mortality_engine <- function(n = 50L, rate_value = 0.05,
                             horizon = c(20000101, 20101231),
                             max_age = 120, ...) {
  codec <- mortality_codec()
  new_engine(codec, constant_registry(rate_value, codec),
             same_age_cohort(n), horizon = horizon, max_age = max_age, ...)
}


# ---------------------------------------------------------------------
# 1. The smallest possible run
# ---------------------------------------------------------------------

test_that("a run returns a population, diagnostics and extensions", {
  set.seed(1L)
  result <- mortality_engine(n = 10L)$run()

  expect_named(result, c("pop", "newborns", "diagnostics", "extensions"))
  expect_identical(result$extensions, setNames(list(), character(0)))
  expect_s3_class(result$pop, "data.frame")
  expect_identical(result$diagnostics$n_individuals, 10L)
})


test_that("an engine can be run more than once without carrying state", {
  # The population register, the queue and the transition record are all
  # created inside run(), so two runs cannot influence each other. This is
  # what makes the parallel backend straightforward later, and it is the
  # opposite of MicSim 3.0.0, whose state lived in the calling frame.
  engine <- mortality_engine(n = 40L)

  first  <- engine$run(seed = 7L)
  second <- engine$run(seed = 7L)

  expect_identical(first$pop, second$pop)
  expect_identical(first$diagnostics$n_events,
                   second$diagnostics$n_events)
})


# ---------------------------------------------------------------------
# 2. The output contract
# ---------------------------------------------------------------------

test_that("pop has the contract columns", {
  set.seed(2L)
  pop <- mortality_engine(n = 20L)$run()$pop
  expect_identical(
    names(pop),
    c("ID", "birthDate", "initState", "From", "To",
      "transitionTime", "transitionAge")
  )
})


test_that("every individual appears, and nobody appears twice without cause", {
  # "dead" is absorbing and it is the only transition, so each individual
  # has at most one event. That makes the row count exactly the population
  # size: one row per death, one NA row per survivor. It is the cleanest
  # possible check that no synthetic entry rows have crept back in
  # (deviation M1).
  set.seed(3L)
  n <- 60L
  pop <- mortality_engine(n = n)$run()$pop

  expect_identical(nrow(pop), n)
  expect_setequal(pop$ID, seq_len(n))
})


test_that("survivors carry NA in all four transition columns", {
  set.seed(4L)
  pop <- mortality_engine(n = 60L)$run()$pop
  survivors <- pop[is.na(pop$From), , drop = FALSE]

  expect_gt(nrow(survivors), 0L)
  expect_true(all(is.na(survivors$To)))
  expect_true(all(is.na(survivors$transitionTime)))
  expect_true(all(is.na(survivors$transitionAge)))
  # but they are still identified
  expect_false(anyNA(survivors$ID))
  expect_false(anyNA(survivors$initState))
})


test_that("the birth date is returned exactly as it was supplied", {
  # The contract hands back the user's own column. Re-deriving it from the
  # day count would silently change their data.
  set.seed(5L)
  pop <- mortality_engine(n = 10L)$run()$pop
  expect_true(all(pop$birthDate == "19900101"))
  expect_type(pop$birthDate, "character")
})


test_that("deaths are recorded as a move from the initial state to dead", {
  set.seed(6L)
  pop <- mortality_engine(n = 60L)$run()$pop
  died <- pop[!is.na(pop$From), , drop = FALSE]

  expect_gt(nrow(died), 0L)
  expect_true(all(died$To == "dead"))
  expect_true(all(died$From %in% c("m", "f")))
  expect_identical(died$From, died$initState)
})


# ---------------------------------------------------------------------
# 3. The loop's invariants
# ---------------------------------------------------------------------

test_that("an absorbing state absorbs", {
  # Nobody may have a second event after entering "dead". With one row per
  # individual guaranteed by the previous section, this is equivalent to
  # no ID appearing twice.
  set.seed(7L)
  pop <- mortality_engine(n = 80L)$run()$pop
  expect_false(any(duplicated(pop$ID)))
})


test_that("no event is recorded after the horizon closes", {
  # The loop stops on the first event that falls past the end date, and
  # because the queue is a min-heap that event is the earliest remaining
  # one, so nothing later can be missed.
  set.seed(8L)
  horizon <- c(20000101, 20101231)
  pop <- mortality_engine(n = 200L, horizon = horizon)$run()$pop
  times <- getInDays(pop$transitionTime[!is.na(pop$transitionTime)])

  expect_true(all(times >= getInDays(horizon[1])))
  expect_true(all(times <= getInDays(horizon[2])))
})


test_that("an individual can pass through several states in order", {
  # A three-stage model, so a life course is more than one event long.
  # This is the first test of the RESCHEDULING path: after a transition
  # the engine must put the individual back in the queue with a new
  # destination drawn from their new state.
  codec <- new_state_codec(list(stage = c("0", "1", "2")))
  quick <- function(age, calTime) rep(2, length(age))
  registry <- new_rate_registry(
    list(
      list(from = "0", to = "1", rate = quick),
      list(from = "1", to = "2", rate = quick)
    ),
    codec,
    probe = micsim_rate_probe(0:120, 2000:2040, 0:1)
  )
  init_pop <- data.frame(
    ID = 1:100, birthDate = rep("19900101", 100L),
    initState = rep("0", 100L), stringsAsFactors = FALSE
  )
  engine <- new_engine(codec, registry, init_pop,
                       horizon = c(20000101, 20301231), max_age = 120)

  set.seed(9L)
  pop <- engine$run()$pop

  completed <- names(which(table(pop$ID) == 2L))
  expect_gt(length(completed), 0L)

  for (id in completed) {
    rows <- pop[as.character(pop$ID) == id, , drop = FALSE]
    expect_identical(rows$From, c("0", "1"))
    expect_identical(rows$To,   c("1", "2"))
    expect_gt(rows$transitionAge[2L], rows$transitionAge[1L])
  }
})


test_that("duration is measured from the change of the affected variable", {
  # The end-to-end proof of the M2 fix. The second transition's rate is
  # zero until five years have passed SINCE ENTERING state 1, so nobody
  # may reach state 2 less than five years after reaching state 1.
  #
  # The mock-up measured duration from the last transition of any kind,
  # which happens to agree here; what it could not do is tell two state
  # variables apart. This test at least pins the clock to the right
  # origin, which is the half of M2 an end-to-end run can observe.
  codec <- new_state_codec(list(stage = c("0", "1", "2")))
  quick <- function(age, calTime) rep(2, length(age))
  after_five_years <- function(duration) ifelse(duration < 5, 0, 5)

  registry <- new_rate_registry(
    list(
      list(from = "0", to = "1", rate = quick),
      list(from = "1", to = "2", rate = after_five_years)
    ),
    codec,
    probe = micsim_rate_probe(0:120, 2000:2040, 0:120)
  )
  init_pop <- data.frame(
    ID = 1:100, birthDate = rep("19900101", 100L),
    initState = rep("0", 100L), stringsAsFactors = FALSE
  )
  engine <- new_engine(codec, registry, init_pop,
                       horizon = c(20000101, 20401231), max_age = 120)

  set.seed(16L)
  pop <- engine$run()$pop

  completed <- names(which(table(pop$ID) == 2L))
  expect_gt(length(completed), 0L)

  for (id in completed) {
    ages <- pop$transitionAge[as.character(pop$ID) == id]
    # 0.02 of slack for the two-decimal rounding of transitionAge
    expect_gte(ages[2L] - ages[1L], 5 - 0.02)
  }
})


test_that("an individual with no reachable transition is simply left alone", {
  # A transient state the model gives no way out of is legitimate: think
  # of a residual category. The registry answers NULL for it and the
  # engine schedules nothing, rather than erroring.
  codec <- new_state_codec(list(sex = c("m", "f")), absorbing = "dead")
  only_women_die <- function(age, calTime) rep(0.5, length(age))
  registry <- new_rate_registry(
    list(list(from = "f", to = "dead", rate = only_women_die)),
    codec
  )
  engine <- new_engine(codec, registry, same_age_cohort(20L),
                       horizon = c(20000101, 20101231), max_age = 120)

  set.seed(10L)
  pop <- engine$run()$pop
  men <- pop[pop$initState == "m", , drop = FALSE]

  expect_identical(nrow(men), 10L)
  expect_true(all(is.na(men$From)))
})


# ---------------------------------------------------------------------
# 4. maxAge and the horizon
# ---------------------------------------------------------------------

test_that("nobody transitions after reaching maxAge", {
  # On reaching maxAge an individual stops being at risk; their episode is
  # censored at the end of the simulation rather than the individual being
  # removed, so they still appear in the output.
  codec <- mortality_codec()
  init_pop <- data.frame(
    ID = 1:40, birthDate = rep("19600101", 40L),
    initState = rep(c("m", "f"), 20L), stringsAsFactors = FALSE
  )
  engine <- new_engine(codec, constant_registry(0.3, codec), init_pop,
                       horizon = c(20000101, 20401231), max_age = 45)

  set.seed(11L)
  pop <- engine$run()$pop
  ages <- pop$transitionAge[!is.na(pop$transitionAge)]

  expect_true(all(ages < 45))
  # and everyone is still listed
  expect_identical(nrow(pop), 40L)
})


test_that("a horizon too short for anything to happen yields no events", {
  codec <- mortality_codec()
  engine <- new_engine(codec, constant_registry(1e-8, codec),
                       same_age_cohort(30L),
                       horizon = c(20000101, 20000102), max_age = 120)

  set.seed(12L)
  result <- engine$run()
  expect_identical(result$diagnostics$n_events, 0L)
  expect_identical(nrow(result$pop), 30L)
  expect_true(all(is.na(result$pop$From)))
})


# ---------------------------------------------------------------------
# 5. Determinism and the random stream
# ---------------------------------------------------------------------

test_that("the same seed gives the same population", {
  engine <- mortality_engine(n = 100L)
  expect_identical(engine$run(seed = 99L)$pop, engine$run(seed = 99L)$pop)
})


test_that("different seeds give different populations", {
  engine <- mortality_engine(n = 100L)
  expect_false(identical(engine$run(seed = 1L)$pop,
                         engine$run(seed = 2L)$pop))
})


test_that("with no seed the engine consumes the caller's stream", {
  # This reproduces the MicSim 3.0.0 idiom exactly: micSim() has no seed
  # argument, so `set.seed(234); micSim(...)` is how a run is made
  # reproducible. The compatibility shim depends on this behaviour.
  engine <- mortality_engine(n = 40L)

  set.seed(4242L)
  first <- engine$run()
  set.seed(4242L)
  second <- engine$run()

  expect_identical(first$pop, second$pop)
})


test_that("without a seed a run takes exactly one draw from the caller's stream", {
  # The draw is the run's seed. Everything else comes from the individuals'
  # own streams, so the caller's stream moves by one draw whatever the
  # model does, here for someone with no window left.
  space <- c("m", "f")
  attr(space, "name") <- "sex"
  tm <- buildTransitionMatrix(NULL, cbind("dead", "e4Rate"), space)
  e4Rate <- function(age, calTime) rep(0.01, length(age))
  # 99.995 years old on 2000-01-01, maxAge 100: window = 99.99 - 99.995 < 0
  init_pop <- data.frame(ID = 1L, birthDate = "19000102", initState = "f",
                         stringsAsFactors = FALSE)

  set.seed(99L)
  pop <- micSim(initPop = init_pop, transitionMatrix = tm,
                absStates = "dead", maxAge = 100,
                simHorizon = c(20000101, 20101231), verbose = FALSE)
  after_sim <- .Random.seed
  set.seed(99L); invisible(sample.int(.Machine$integer.max, 1L))
  expect_identical(after_sim, .Random.seed)
  expect_true(is.na(pop$To))
})


test_that("supplying a seed restores the caller's stream afterwards", {
  # A simulation must not quietly reposition the random number generator
  # of the script that called it (checklist item T7.9).
  engine <- mortality_engine(n = 20L)

  set.seed(5L)
  before <- .Random.seed
  engine$run(seed = 12345L)
  after <- .Random.seed

  expect_identical(before, after)
})


# ---------------------------------------------------------------------
# 6. Statistical correctness
# ---------------------------------------------------------------------

test_that("a constant hazard reproduces the exponential survival function", {
  # THE test. With a constant hazard lambda and everyone starting at the
  # same age, the probability of dying within a window of T years is
  # 1 - exp(-lambda * T), exactly. Simulating a large cohort and counting
  # deaths is then a binomial experiment with a known success probability,
  # and the observed proportion must sit within a few standard errors of
  # it.
  #
  # Nothing short of the whole engine can produce this number: it needs
  # the solver, the registry, the queue, the person store and the loop to
  # all be right simultaneously.
  skip_on_cran()

  rate_value <- 0.05
  n <- 4000L
  start <- 20000101
  end   <- 20091231

  codec <- mortality_codec()
  engine <- new_engine(codec, constant_registry(rate_value, codec),
                       same_age_cohort(n),
                       horizon = c(start, end), max_age = 120)

  set.seed(20260915L)
  pop <- engine$run()$pop

  window_years <- (getInDays(end) - getInDays(start)) / 365.25
  expected_p <- 1 - exp(-rate_value * window_years)
  observed_p <- mean(!is.na(pop$From))
  standard_error <- sqrt(expected_p * (1 - expected_p) / n)

  expect_lt(abs(observed_p - expected_p), 4 * standard_error)
})


test_that("ages at death follow the exponential waiting-time distribution", {
  # The same experiment looked at from the other side: not how many died
  # but when. Conditional on dying inside the window, the waiting time is
  # exponential truncated to that window, and a Kolmogorov-Smirnov test
  # against that law should not reject.
  skip_on_cran()

  rate_value <- 0.08
  n <- 3000L
  start <- 20000101
  end   <- 20291231

  codec <- mortality_codec()
  engine <- new_engine(codec, constant_registry(rate_value, codec),
                       same_age_cohort(n, birth_ymd = "20000101"),
                       horizon = c(start, end), max_age = 120)

  set.seed(20260916L)
  pop <- engine$run()$pop

  # Everyone is born on the simulation start date, so age at death is the
  # waiting time.
  waiting <- pop$transitionAge[!is.na(pop$transitionAge)]
  window_years <- (getInDays(end) - getInDays(start)) / 365.25

  truncated_exponential <- function(q) {
    stats::pexp(q, rate = rate_value) /
      stats::pexp(window_years, rate = rate_value)
  }
  ks <- suppressWarnings(stats::ks.test(waiting, truncated_exponential))
  expect_gt(ks$p.value, 0.01)
})


test_that("a higher hazard kills more people", {
  # A sanity check that costs nothing and would catch a rate that was
  # being ignored, mis-scaled, or evaluated on the wrong time scale.
  low  <- mortality_engine(n = 500L, rate_value = 0.01)
  high <- mortality_engine(n = 500L, rate_value = 0.20)

  set.seed(13L); n_low  <- sum(!is.na(low$run()$pop$From))
  set.seed(13L); n_high <- sum(!is.na(high$run()$pop$From))

  expect_gt(n_high, n_low)
})


# ---------------------------------------------------------------------
# 7. Diagnostics
# ---------------------------------------------------------------------

test_that("diagnostics report what the run did", {
  set.seed(14L)
  result <- mortality_engine(n = 50L)$run(seed = 21L)
  d <- result$diagnostics

  expect_identical(d$n_events, nrow(result$pop[!is.na(result$pop$From), ]))
  expect_identical(d$n_individuals, 50L)
  expect_identical(d$n_unreliable, 0L)   # a constant hazard integrates cleanly
  expect_identical(d$n_forced, 0L)       # no Inf rates in this model
  expect_identical(d$seed, 21L)
  expect_gt(d$peak_queue, 0L)
  expect_type(d$seconds, "double")
})


test_that("a seed drawn from the caller's stream is reported and reproduces the run", {
  engine <- mortality_engine(n = 5L)
  set.seed(15L)
  run <- engine$run()
  set.seed(15L)
  expected <- sample.int(.Machine$integer.max, 1L)

  expect_identical(run$diagnostics$seed, expected)
  expect_true(run$diagnostics$seed_drawn)
  expect_identical(engine$run(seed = expected)$pop, run$pop)
  expect_false(engine$run(seed = expected)$diagnostics$seed_drawn)
})


# ---------------------------------------------------------------------
# 8. Optional features are off unless asked for
# ---------------------------------------------------------------------

test_that("a model with no optional features installs no subsystem", {
  # The seam costs nothing when nothing is installed: one is.null() test
  # per hook point per event, and no context is ever built.
  set.seed(16L)
  result <- mortality_engine(n = 20L)$run()

  expect_identical(result$diagnostics$n_immigrants, 0L)
  expect_null(result$diagnostics$n_clock_reversals)
  expect_false("motherID" %in% names(result$pop))
})


test_that("the engine takes only compiled subsystems", {
  # Subsystems reach the engine compiled against the model, as
  # add_subsystem() stores them. A raw specification, or anything else, is
  # refused at the boundary rather than failing inside the loop.
  codec <- mortality_codec()
  registry <- constant_registry(0.05, codec)
  args <- list(codec = codec, registry = registry,
               init_pop = same_age_cohort(5L),
               horizon = c(20000101, 20101231), max_age = 120)

  expect_error(do.call(new_engine, c(args, list(subsystems = list("x")))),
               class = "micsim_engine_error")
  expect_error(
    do.call(new_engine, c(args, list(subsystems = list(
      school_enrolment(move("sex", "m", "f"), 9L))))),
    class = "micsim_engine_error"
  )
})


test_that("a malformed immigrant table is refused", {
  codec <- mortality_codec()
  registry <- constant_registry(0.05, codec)

  expect_error(
    new_engine(codec, registry, same_age_cohort(5L),
               horizon = c(20000101, 20101231), max_age = 120,
               immigrants = data.frame(ID = 1L)),
    class = "micsim_engine_error"
  )
})


test_that("school enrolment for a transition the model lacks is refused", {
  # Accepting it silently would mean the feature never fires.
  codec <- new_state_codec(list(sex = c("m", "f"), edu = c("no", "low")),
                           absorbing = "dead")
  registry <- constant_registry(0.05, codec)

  expect_error(
    compiled(school_enrolment(move("edu", "no", "low"), 9L), codec, registry),
    class = "micsim_selector_error"
  )
})


# ---------------------------------------------------------------------
# 9. Input validation
# ---------------------------------------------------------------------

test_that("a malformed starting population is rejected", {
  codec <- mortality_codec()
  registry <- constant_registry(0.05, codec)
  build <- function(init_pop) {
    new_engine(codec, registry, init_pop,
               horizon = c(20000101, 20101231), max_age = 120)
  }

  expect_error(build(data.frame(id = 1L, birth = "19900101", st = "m")),
               class = "micsim_engine_error")
  expect_error(build(same_age_cohort(0L)), class = "micsim_engine_error")
})


test_that("an impossible horizon or maxAge is rejected", {
  codec <- mortality_codec()
  registry <- constant_registry(0.05, codec)

  expect_error(
    new_engine(codec, registry, same_age_cohort(5L),
               horizon = c(20101231, 20000101), max_age = 120),
    class = "micsim_engine_error"
  )
  expect_error(
    new_engine(codec, registry, same_age_cohort(5L),
               horizon = c(20000101, 20101231), max_age = -1),
    class = "micsim_engine_error"
  )
})


test_that("an individual born after the simulation starts is rejected", {
  # MicSim 3.0.0 stopped on this too, listing the offending IDs. Being
  # born during the run is what fertility is for. Refused when the engine
  # is built, so a model that exists can run.
  codec <- mortality_codec()
  init_pop <- data.frame(
    ID = 1:3,
    birthDate = c("19900101", "20500101", "19950101"),
    initState = c("m", "f", "m"),
    stringsAsFactors = FALSE
  )
  expect_error(new_engine(codec, constant_registry(0.05, codec), init_pop,
                       horizon = c(20000101, 20101231), max_age = 120),
               class = "micsim_engine_error")
})


test_that("an individual already older than maxAge is rejected", {
  codec <- mortality_codec()
  init_pop <- data.frame(
    ID = 1:2, birthDate = c("19900101", "19000101"),
    initState = c("m", "f"), stringsAsFactors = FALSE
  )
  expect_error(new_engine(codec, constant_registry(0.05, codec), init_pop,
                       horizon = c(20000101, 20101231), max_age = 60),
               class = "micsim_engine_error")
})


test_that("an unknown state label in the starting population is rejected", {
  codec <- mortality_codec()
  init_pop <- data.frame(
    ID = 1:2, birthDate = c("19900101", "19910101"),
    initState = c("m", "x"), stringsAsFactors = FALSE
  )
  expect_error(new_engine(codec, constant_registry(0.05, codec), init_pop,
                          horizon = c(20000101, 20101231), max_age = 120),
               "x", class = "micsim_codec_error")
})
