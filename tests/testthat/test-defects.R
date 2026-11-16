# =====================================================================
# test-defects.R
# ---------------------------------------------------------------------
# One test per defect fixed for the 4.0.0 release that is not pinned
# elsewhere. Each says what MicSim 3.0.0 (or an earlier 4.0.0 build) did,
# so that a regression names the defect it brings back. Forced transitions
# and the school calendar are pinned in test-enrolment.R and
# test-registry.R; the all-zero shortcut and quadrature reporting in
# test-solver.R.
# =====================================================================


def_mort <- function(age, calTime) 0.00005 * exp(0.09 * age)


test_that("episodes are censored at max_age, not at the end of the run", {
  # A person reaching max_age before the horizon stops being at risk;
  # MicSim 3.0.0 still counted exposure to the end of the simulation.
  model <- micsim_model(
    state_space(s = c("a", "b")),
    transitions(transition("a", "b", function(age) rep(0, length(age)))),
    data.frame(ID = 1:2, birthDate = c("19450101", "19900101"),
               initState = "a"),
    c(20200101, 20501231), max_age = 80
  )
  long <- as_long(simulate(model, seed = 1L))
  stop_1 <- getInDays(long$Tstop[long$ID == 1L])
  expect_equal((stop_1 - getInDays("19450101")) / DAYS_PER_YEAR,
               80 - MAX_AGE_MARGIN, tolerance = 1 / DAYS_PER_YEAR)
  expect_identical(long$Tstop[long$ID == 2L], "20501231")   # horizon binds
  expect_true(all(long$statusExit == 0L))
})


test_that("the same censoring applies to a transition followed by survival", {
  model <- micsim_model(
    state_space(s = c("a", "b")),
    transitions(transition("a", "b", function(age) rep(5, length(age)))),
    data.frame(ID = 1L, birthDate = "19450101", initState = "a"),
    c(20200101, 20501231), max_age = 80
  )
  long <- as_long(simulate(model, seed = 1L))
  expect_identical(nrow(long), 2L)
  expect_lt(getInDays(long$Tstop[2L]), getInDays(20250101))
})


test_that("classic results record maxAge for the converter", {
  ss <- expand.grid(s = c("a", "b"), stringsAsFactors = FALSE)
  none <- function(age) rep(0, length(age))
  tm <- buildTransitionMatrix(cbind("a->b", "none"), cbind("dead", "none"), ss)
  pop <- micSim(data.frame(ID = 1L, birthDate = "19450101", initState = "a"),
                transitionMatrix = tm, maxAge = 80,
                simHorizon = c(20200101, 20501231), verbose = FALSE)
  expect_identical(attr(pop, "micsim_spec")$maxAge, 80)
  long <- convertToLongFormat(pop, stateSpace = ss)
  expect_lt(getInDays(long$Tstop), getInDays(20250101))
  # without the record, the old behaviour
  bare <- pop; attr(bare, "micsim_spec") <- NULL
  long_old <- convertToLongFormat(bare, stateSpace = ss, absStates = "dead",
                                  simHorizon = c(20200101, 20501231))
  expect_identical(long_old$Tstop, "20501231")
})


test_that("validate() checks durations up to the oldest age", {
  # Duration clocks of the starting population start at birth, so a
  # 94-year-old has duration 94 at once. The grid used to stop at the
  # spread of ages (here 5) and let this rate through.
  model <- micsim_model(
    state_space(sex = c("m", "f"), mar = c("s", "w")),
    transitions(transition(move("mar", "s", "w"), function(age, duration)
      ifelse(duration > 80, -1, 0.01))),
    data.frame(ID = 1:3, birthDate = "19050101", initState = "f/s"),
    c(20000101, 20101231), max_age = 100
  )
  expect_identical(max(reachable_grid(model)$durations), 99L)
  expect_error(validate(model), class = "micsim_rate_error")
})


test_that("ambiguous classic state values are warned about", {
  ss <- expand.grid(sex = c("1", "2"), parity = c("0", "1", "2"),
                    stringsAsFactors = FALSE)
  expect_warning(buildTransitionMatrix(cbind("1->2", "r"), cbind("dead", "m"), ss),
                 class = "micsim_ambiguous_values_warning")
  ok <- expand.grid(sex = c("m", "f"), parity = c("0", "1"),
                    stringsAsFactors = FALSE)
  expect_no_warning(buildTransitionMatrix(cbind("0->1", "r"),
                                          cbind("dead", "m"), ok))
})


test_that("identifiers must be whole numbers", {
  # as.integer() used to truncate 1.5 to 1 without a word.
  build <- function(ids) {
    micsim_model(state_space(sex = c("m", "f")),
                 transitions(transition(NA, "dead", def_mort)),
                 data.frame(ID = ids, birthDate = "19800101", initState = "f"),
                 c(20000101, 20101231))
  }
  expect_error(build(c(1, 1.5)), "whole numbers", class = "micsim_engine_error")
  expect_error(build(c(1, 1)), "repeated", class = "micsim_engine_error")
  expect_s3_class(build(c(1, 2)), "micsim_model")
})


test_that("an enrolment month must be a whole number", {
  e <- move("edu", "no", "low")
  expect_error(school_enrolment(e, 9.7), class = "micsim_enrolment_error")
  expect_error(school_enrolment(e, "9"), class = "micsim_enrolment_error")
  expect_s3_class(school_enrolment(e, 9), "micsim_subsystem")
})


test_that("nobody may enter the population already absorbed", {
  expect_error(
    micsim_model(state_space(sex = c("m", "f")),
                 transitions(transition(NA, "dead", def_mort)),
                 data.frame(ID = 1:2, birthDate = "19800101",
                            initState = c("f", "dead")),
                 c(20000101, 20101231)),
    class = "micsim_engine_error"
  )
})


test_that("micSimParallel() records what the converter needs", {
  skip_on_cran()
  ss <- expand.grid(sex = c("m", "f"), stringsAsFactors = FALSE)
  mortRates <- def_mort
  tm <- buildTransitionMatrix(NULL, cbind("dead", "mortRates"), ss)
  pop <- micSimParallel(
    initPop = data.frame(ID = 1:20, birthDate = "19500101", initState = "f"),
    transitionMatrix = tm, maxAge = 100, simHorizon = c(20000101, 20201231),
    cores = 2L, seeds = 3L, verbose = FALSE
  )
  spec <- attr(pop, "micsim_spec")
  expect_identical(spec$absStates, "dead")
  expect_identical(spec$maxAge, 100)
  expect_s3_class(convertToLongFormat(pop, stateSpace = ss), "data.frame")
})
