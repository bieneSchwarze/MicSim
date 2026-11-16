# =====================================================================
# test-registry.R
# ---------------------------------------------------------------------
# The registry is where a model's transitions stop being text and become
# something the engine can execute. Everything it decides is decided ONCE,
# before the run, so everything it decides is testable without running a
# simulation: which states an edge connects, which time scales its rate
# uses, which state variable it changes, and whether it can be forced.
#
# Tests are grouped by what they establish:
#   1. construction and indexing
#   2. the from = NA shorthand
#   3. time-scale dependence from formals()
#   4. the affected state variable (what `duration` is measured against)
#   5. rejection of malformed rates
#   6. bind_hazard: the three clocks, and the truncated probe
#   7. may_be_infinite probing
# =====================================================================


toy_codec <- function() {
  new_state_codec(
    state_space = list(sex = c("m", "f"), fert = c("0", "1")),
    absorbing   = "dead"
  )
}

mortality <- function(age, calTime) 0.00003 * exp(0.1 * age)
fertility <- function(age, calTime) ifelse(age > 15 & age < 45, 0.1, 0)


# ---------------------------------------------------------------------
# 1. Construction and indexing
# ---------------------------------------------------------------------

test_that("a registry indexes edges by origin state", {
  codec <- toy_codec()
  reg <- new_rate_registry(
    transitions = list(
      list(from = "f/0", to = "f/1",  rate = fertility),
      list(from = "f/0", to = "dead", rate = mortality)
    ),
    codec = codec
  )

  out <- reg$feasible(codec$encode("f/0"))
  expect_length(out, 2L)
  expect_setequal(vapply(out, function(e) e$to_label, character(1L)),
                  c("f/1", "dead"))

  # A state nobody can leave reports NULL, not an empty list: NULL is the
  # honest answer and the engine tests for it.
  expect_null(reg$feasible(codec$encode("m/0")))
})


test_that("a non-codec first argument is rejected", {
  expect_error(
    new_rate_registry(
      transitions = list(list(from = "a", to = "b", rate = mortality)),
      codec = list()
    ),
    class = "micsim_registry_error"
  )
})


test_that("an unknown state label is rejected", {
  expect_error(
    new_rate_registry(
      transitions = list(list(from = "f/9", to = "dead", rate = mortality)),
      codec = toy_codec()
    ),
    class = "micsim_registry_error"
  )
})


# ---------------------------------------------------------------------
# 2. The from = NA shorthand
# ---------------------------------------------------------------------

test_that("from = NA expands to every transient state", {
  # This is how mortality is normally declared: everyone is always at risk
  # of dying. Four transient states here, so four edges from one
  # specification.
  codec <- toy_codec()
  reg <- new_rate_registry(
    transitions = list(list(from = NA, to = "dead", rate = mortality)),
    codec = codec
  )

  expect_length(reg$edges, length(codec$transient))
  for (state in codec$transient) {
    edges <- reg$feasible(codec$encode(state))
    expect_length(edges, 1L)
    expect_identical(edges[[1L]]$to_label, "dead")
  }
})


test_that("from = NA into a transient state leaves out the state itself", {
  # "From every state to f/1" means from every OTHER state. A self-loop
  # f/1 -> f/1 would be simulated as an event that changes nothing.
  codec <- toy_codec()
  reg <- new_rate_registry(
    transitions = list(list(from = NA, to = "f/1", rate = fertility)),
    codec = codec
  )
  expect_length(reg$edges, length(codec$transient) - 1L)
  expect_null(reg$feasible(codec$encode("f/1")))
  from <- vapply(reg$edges, `[[`, character(1L), "from_label")
  expect_false("f/1" %in% from)
})


test_that("a transition into its own origin is refused", {
  expect_error(
    new_rate_registry(
      transitions = list(list(from = c("f/0", "f/1"), to = "f/1",
                              rate = fertility)),
      codec = toy_codec()
    ),
    class = "micsim_registry_error"
  )
  one_state <- new_state_codec(state_space = list(sex = "f"),
                               absorbing = "dead")
  expect_error(
    new_rate_registry(
      transitions = list(list(from = NA, to = "f", rate = fertility)),
      codec = one_state
    ),
    class = "micsim_registry_error"
  )
})


# ---------------------------------------------------------------------
# 3. Time-scale dependence
# ---------------------------------------------------------------------

test_that("dependence on age, calTime and duration comes from formals()", {
  codec <- toy_codec()
  dur_dependent <- function(age, duration) ifelse(duration < 0.75, 0, 0.1)

  reg <- new_rate_registry(
    transitions = list(
      list(from = "f/0", to = "f/1",  rate = dur_dependent),
      list(from = "f/0", to = "dead", rate = mortality)
    ),
    codec = codec
  )

  fert_edge <- reg$edges[[1L]]
  expect_true(fert_edge$uses_age)
  expect_false(fert_edge$uses_caltime)
  expect_true(fert_edge$uses_duration)

  mort_edge <- reg$edges[[2L]]
  expect_true(mort_edge$uses_age)
  expect_true(mort_edge$uses_caltime)
  expect_false(mort_edge$uses_duration)
})


# ---------------------------------------------------------------------
# 4. The affected state variable
# ---------------------------------------------------------------------

test_that("changed_vars lists every component the transition changes", {
  # affected_var answers "which variable does duration refer to" and is NA
  # unless exactly one changes. changed_vars answers "which clocks must be
  # reset" and lists all of them. The person store needs the second.
  expect_identical(changed_state_variables("f/DE/0", "f/DE/1"), 3L)
  expect_identical(changed_state_variables("f/DE/0", "m/AT/0"), c(1L, 2L))
  expect_identical(changed_state_variables("f/DE/0", "f/DE/0"), integer(0))
  expect_identical(changed_state_variables("f/DE/0", "dead"), integer(0))
})


test_that("every edge carries both affected_var and changed_vars", {
  codec <- new_state_codec(
    state_space = list(sex = c("m", "f"), nat = c("DE", "AT")),
    absorbing   = "dead"
  )
  flat <- function(age, calTime) rep(0.1, length(age))

  reg <- new_rate_registry(
    transitions = list(
      list(from = "f/DE", to = "f/AT", rate = flat),   # one variable
      list(from = "f/DE", to = "m/AT", rate = flat),   # two variables
      list(from = "f/DE", to = "dead", rate = flat)    # absorbing
    ),
    codec = codec
  )

  expect_identical(reg$edges[[1L]]$affected_var, 2L)
  expect_identical(reg$edges[[1L]]$changed_vars, 2L)

  expect_true(is.na(reg$edges[[2L]]$affected_var))
  expect_identical(reg$edges[[2L]]$changed_vars, c(1L, 2L))

  expect_true(is.na(reg$edges[[3L]]$affected_var))
  expect_identical(reg$edges[[3L]]$changed_vars, integer(0))
})


test_that("a duration-dependent multi-variable transition is rejected", {
  # "Time since the last change of the affected variable" is undefined
  # when two variables change at once. MicSim 3.0.0 printed a note with
  # cat() and continued with Inf; here it is a model error.
  codec <- new_state_codec(
    state_space = list(sex = c("m", "f"), nat = c("DE", "AT")),
    absorbing   = "dead"
  )
  dur_dependent <- function(age, duration) rep(0.1, length(age))

  expect_error(
    new_rate_registry(
      transitions = list(list(from = "f/DE", to = "m/AT", rate = dur_dependent)),
      codec = codec
    ),
    class = "micsim_rate_error"
  )
})


test_that("a duration-dependent transition into an absorbing state is allowed", {
  # Absorbing destinations have an atomic label, get duration = Inf, and
  # that is exactly what MicSim 3.0.0 does. Nothing is undefined, so
  # nothing is rejected.
  codec <- toy_codec()
  dur_dependent <- function(age, duration) rep(0.1, length(age))
  expect_no_error(
    new_rate_registry(
      transitions = list(list(from = "f/0", to = "dead", rate = dur_dependent)),
      codec = codec
    )
  )
})


# ---------------------------------------------------------------------
# 5. Rejection of malformed rates
# ---------------------------------------------------------------------

test_that("a rate passed as a string is rejected", {
  # MicSim 3.0.0 stored rate function NAMES and resolved them on every
  # event through the calling scope (defect D12). Names are no longer
  # accepted here; the compatibility shim resolves them once, up front.
  expect_error(
    new_rate_registry(
      transitions = list(list(from = "f/0", to = "dead", rate = "mortality")),
      codec = toy_codec()
    ),
    class = "micsim_rate_error"
  )
})


test_that("an illegal argument name is rejected", {
  bad_args <- function(age, income) rep(0.1, length(age))
  expect_error(
    new_rate_registry(
      transitions = list(list(from = "f/0", to = "dead", rate = bad_args)),
      codec = toy_codec()
    ),
    class = "micsim_rate_error"
  )
})


test_that("a non-vectorised rate is rejected at construction", {
  # Caught here, when the model is built, rather than minutes into a run.
  not_vectorised <- function(age, calTime) 0.1
  expect_error(
    new_rate_registry(
      transitions = list(list(from = "f/0", to = "dead", rate = not_vectorised)),
      codec = toy_codec()
    ),
    class = "micsim_rate_error"
  )
})


# ---------------------------------------------------------------------
# 6. bind_hazard
# ---------------------------------------------------------------------

test_that("bind_hazard advances all three clocks together", {
  # An individual aged 30 in calendar year 2014 who has been in the
  # affected state for 2 years is, x years later, aged 30 + x in year
  # 2014 + x with duration 2 + x. bind_hazard() is the only place that
  # threading happens.
  reveal <- function(age, calTime, duration) age + 1000 * (calTime - 2014) +
    1e6 * duration

  edge <- list(
    rate = reveal, argnames = c("age", "calTime", "duration"),
    uses_age = TRUE, uses_caltime = TRUE, uses_duration = TRUE,
    may_be_infinite = FALSE
  )
  h <- bind_hazard(edge, age0_years = 30, caltime0_years = 2014,
                   duration0_years = 2)

  expect_equal(h$lambda(0), 30 + 0 + 2e6)
  expect_equal(h$lambda(5), 35 + 5000 + 7e6)
  expect_null(h$forced_time)
})


test_that("bind_hazard only supplies the arguments the rate declares", {
  age_only <- function(age) age
  edge <- list(
    rate = age_only, argnames = "age",
    uses_age = TRUE, uses_caltime = FALSE, uses_duration = FALSE,
    may_be_infinite = FALSE
  )
  h <- bind_hazard(edge, age0_years = 30, caltime0_years = 2014,
                   duration0_years = Inf)
  expect_equal(h$lambda(c(0, 1, 2)), c(30, 31, 32))
})


test_that("a forced age is found exactly, whatever the current age", {
  # School-enrolment rates are written as ifelse(age == 7, Inf, 0), an
  # exact equality that only ever holds on whole years. MicSim 3.0.0
  # looked for it at trunc(age0) + k / 365.25, which is a whole number
  # only every 1461 days: a child aged 6.4 was found at x = 4 (age 10.4),
  # and a newborn never. The finder sets age to exactly 7 at its crossing.
  enrol <- function(age) ifelse(age == 7, Inf, 0)
  edge <- list(
    rate = enrol, argnames = "age", label = "no -> low",
    uses_age = TRUE, uses_caltime = FALSE, uses_duration = FALSE,
    may_be_infinite = TRUE
  )
  at <- function(age0, H = 20) {
    bind_hazard(edge, age0_years = age0, caltime0_years = 2014,
                duration0_years = Inf)$forced_time(H)
  }
  expect_equal(at(6.4), 0.6)
  expect_equal(at(0), 7)                  # a newborn
  expect_equal(at(3.5), 3.5)
  expect_identical(at(7.4), 0)            # completed age 7: forced now
  expect_true(is.na(at(8.2)))             # past it
  expect_true(is.na(at(0, H = 6.5)))      # not within the window
})


test_that("an interval of Inf is found at its first day", {
  from_seven <- function(age) ifelse(age >= 7.5, Inf, 0)
  edge <- list(rate = from_seven, argnames = "age", label = "x",
               uses_age = TRUE, uses_caltime = FALSE, uses_duration = FALSE,
               may_be_infinite = TRUE)
  t <- bind_hazard(edge, 6, 2014, Inf)$forced_time(5)
  expect_equal(t, 1.5, tolerance = 1 / 365.25)
})


test_that("forced calendar years and durations are found at their crossing", {
  in_2020 <- function(calTime) ifelse(calTime == 2020, Inf, 0)
  after_5 <- function(duration) ifelse(duration == 5, Inf, 0)
  e1 <- list(rate = in_2020, argnames = "calTime", label = "a",
             uses_age = FALSE, uses_caltime = TRUE, uses_duration = FALSE,
             may_be_infinite = TRUE)
  e2 <- list(rate = after_5, argnames = "duration", label = "b",
             uses_age = FALSE, uses_caltime = FALSE, uses_duration = TRUE,
             may_be_infinite = TRUE)
  expect_equal(bind_hazard(e1, 30, 2014.25, Inf)$forced_time(10), 5.75)
  expect_equal(bind_hazard(e2, 30, 2014, 1.5)$forced_time(10), 3.5)
})


test_that("the registry probes each rate once and skips the search for finite ones", {
  codec <- new_state_codec(list(edu = c("no", "low"), sex = c("m", "f")),
                           absorbing = "dead")
  probe <- forced_probe_grid(100, c(20000101, 20101231))
  registry <- new_rate_registry(
    list(list(from = NA, to = "dead",
              rate = function(age, calTime) rep(0.01, length(age))),
         list(move = move("edu", "no", "low"),
              rate = function(age) ifelse(age == 7, Inf, 0))),
    codec, probe = probe
  )
  infinite <- vapply(registry$edges, `[[`, logical(1L), "may_be_infinite")
  labels <- vapply(registry$edges, `[[`, character(1L), "label")
  expect_false(any(infinite[grepl("dead", labels)]))
  expect_true(all(infinite[!grepl("dead", labels)]))
})


test_that("a rate declaring no arguments is recycled to the input length", {
  # A constant written as function() 0.1 is legal. The solver requires one
  # value per element of x, so bind_hazard() recycles.
  edge <- list(
    rate = function() 0.1, argnames = character(0),
    uses_age = FALSE, uses_caltime = FALSE, uses_duration = FALSE,
    may_be_infinite = FALSE
  )
  h <- bind_hazard(edge, 30, 2014, Inf)
  expect_equal(h$lambda(c(0, 1, 2)), rep(0.1, 3L))
})


# ---------------------------------------------------------------------
# 7. may_be_infinite probing
# ---------------------------------------------------------------------

test_that("without a probe grid every edge is conservatively forceable", {
  # TRUE is the safe default: wrongly ruling Inf out would make forced
  # transitions silently never fire, which is far worse than paying for a
  # daily grid scan.
  reg <- new_rate_registry(
    transitions = list(list(from = "f/0", to = "dead", rate = mortality)),
    codec = toy_codec()
  )
  expect_true(reg$edges[[1L]]$may_be_infinite)
})


test_that("a probe grid distinguishes forceable from ordinary rates", {
  codec <- toy_codec()
  probe <- micsim_rate_probe(ages = 0:20, cal_years = 2000:2005,
                             durations = 0:5)

  enrol <- function(age) ifelse(age == 7, Inf, 0)
  reg <- new_rate_registry(
    transitions = list(
      list(from = "f/0", to = "f/1",  rate = enrol),
      list(from = "f/0", to = "dead", rate = mortality)
    ),
    codec = codec, probe = probe
  )

  expect_true(reg$edges[[1L]]$may_be_infinite)
  expect_false(reg$edges[[2L]]$may_be_infinite)
})


test_that("micsim_rate_probe thins the duration grid under the point cap", {
  probe <- micsim_rate_probe(ages = 0:120, cal_years = 1900:2100,
                             durations = 0:120, max_points = 1000)
  total <- length(probe$ages) * length(probe$cal_years) *
    length(probe$durations)
  expect_lte(length(probe$durations), 121L)
  expect_gte(length(probe$durations), 2L)
  expect_lt(total, 121 * 201 * 121)
})
